from collections import defaultdict
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pandas as pd
import numpy as np
import math
import os
import scipy.stats as stats

from minigene_analysis import minigene_utils

class seq_run:
    def __init__(self, data_folder, f_adaptor, r_adaptor, target_seqs, sample_info_file, target_info_file,
                 mdf_file='md5sum_list.txt', threads=4, cutadapt_global_out_dir = 'cutadapt_output',
                 star_index_folder = 'STAR_indices',mapped_folder = 'mapped_reads',
                 assembled_folder = 'assembled_transcripts', quant_folder='transcript_counts'):
        """
        Constructor for seq_run class which will hold information for the overall experiment
        """
        self.data_folder = data_folder
        self.f_adaptor = f_adaptor
        self.r_adaptor = r_adaptor
        self.target_seqs = minigene_utils.convertFastaToDict(target_seqs)
        self.threads = threads
        self.cutadapt_out_dir = cutadapt_global_out_dir
        self.star_index_folder = star_index_folder
        self.mapped_folder = mapped_folder
        self.quant_folder = quant_folder
        self.assembled_folder = assembled_folder

        #get info for all targets
        #target_info maps {reference_name:{transcript_name:[exon_boundaries]}}
        self.target_info_file = target_info_file
        self.target_info = defaultdict(dict)
        self.parse_target_info_file()

        self.reference_gtf_files = {} #{ref_name:ref_gtf_filename}
        self.make_reference_gtfs()

        #get info for all samples
        self.sample_info_file = sample_info_file
        self.reference_names=None
        self.samples = {}  # map sample_name to sample objects
        # map pairs of files to all pairs of barcodes used for them,
        # along with the sample name for that combo
        # {(file_1,file2):{(bar1,bar2):sample_name}}
        self.file_barcode_mapping = defaultdict(dict)
        self.reference_names = set()
        self.expected_MD5s = minigene_utils.get_MD5s_From_file(os.path.join(data_folder, mdf_file))
        self.reference_to_samples = defaultdict(list) #maps reference names to list of samples
        self.merged_stringtie_gtfs = {} #maps reference names stringtie merged GTF file
        self.parse_sample_info_file()

        self.split_and_trim_reads()
        self.map_reads()
        self.assemble_transcripts()
        self.map_reads_for_quant()
        minigene_utils.make_dir('tables')
        self.make_tables('tables/')

    def parse_sample_info_file(self, header_lines=1):
        file_barcode_tuples = []
        with open(self.sample_info_file) as f:
            lines = f.readlines()
            for line in lines[1:]:
                if not line.strip() == '':  # ignore whitespace lines
                    ll = line.strip().split('\t')
                    sample_name, f_reads, r_reads, wt_seq_name, f_bar, r_bar = ll
                    self.reference_names.add(wt_seq_name)
                    f_reads_md5 = minigene_utils.compute_md5(os.path.join(self.data_folder, f_reads))
                    r_reads_md5 = minigene_utils.compute_md5(os.path.join(self.data_folder, r_reads))
                    assert (self.expected_MD5s[f_reads] == f_reads_md5)
                    assert (self.expected_MD5s[r_reads] == r_reads_md5)
                    file_barcode_tuples.append((f_reads, r_reads, f_bar, r_bar))
                    self.file_barcode_mapping[(f_reads, r_reads)][(f_bar, r_bar)] = sample_name
                    self.samples[sample_name] = sample(self, sample_name, f_reads, r_reads, wt_seq_name, f_bar, r_bar)
                    self.reference_to_samples[wt_seq_name].append(self.samples[sample_name])
        # check that no set of files and barcodes appears twice (this would imply that 2 samples
        # with same barcodes are in the same file and thus cannot be seperated)
        assert len(file_barcode_tuples) == len(set(file_barcode_tuples))
        # check that all sample names are unique
        assert len(file_barcode_tuples) == len(set(file_barcode_tuples))

    def parse_target_info_file(self, header_lines=1):
        #for each target sequence, compile all given transcript info
        with open(self.target_info_file) as f:
            lines = f.readlines()
            for line in lines[1:]:
                if not line.strip() == '':  # ignore whitespace lines
                    ll = line.strip().split('\t')
                    reference_name, transcript_name, exons = ll
                    self.target_info[reference_name][transcript_name] = minigene_utils.get_exon_boundaries(exons)

    def split_and_trim_reads(self):
        # for use with cutadapt 4.0
        # python3 -m pip install --user --upgrade cutadapt
        """
        cutadapt \
        -e 0.15 --no-indels \
        -g file:barcodes_fwd.fasta \
        -G file:barcodes_rev.fasta \
        -o {name1}-{name2}.1.fastq.gz -p {name1}-{name2}.2.fastq.gz \
        input.1.fastq.gz input.2.fastq.gz
            """
        minigene_utils.make_dir(self.cutadapt_out_dir)
        for file_pair in self.file_barcode_mapping:
            read1_file, read2_file = file_pair
            read1_path = os.path.join(self.data_folder, read1_file)
            read2_path = os.path.join(self.data_folder, read2_file)
            pair_name = read1_file.split('.')[0]
            cutadapt_file_out_dir = os.path.join(self.cutadapt_out_dir, pair_name)
            if not minigene_utils.file_exists(cutadapt_file_out_dir):
                minigene_utils.make_dir(cutadapt_file_out_dir)
                barcodes = self.file_barcode_mapping[file_pair].keys()
                barfile_f_name = f'{pair_name}_forward_bar.fa'
                barfile_r_name = f'{pair_name}_reverse_bar.fa'
                barfile_f_name = os.path.join(cutadapt_file_out_dir, barfile_f_name)
                barfile_r_name = os.path.join(cutadapt_file_out_dir, barfile_r_name)
                with open(barfile_f_name, 'w') as barfile_f, open(barfile_r_name, 'w') as barfile_r:
                    used_barf = []
                    used_barr = []
                    for barcode_f, barcode_r in barcodes:
                        if barcode_f not in used_barf:
                            barfile_f.write('>%s\n^%s%s\n' % (barcode_f, barcode_f, self.f_adaptor))
                            used_barf.append(barcode_f)
                        if barcode_r not in used_barr:
                            barfile_r.write('>%s\n^%s%s\n' % (barcode_r, barcode_r, self.r_adaptor))
                            used_barr.append(barcode_r)

                f_out = os.path.join(cutadapt_file_out_dir, '%s_{name1}-{name2}.1.fastq' % (pair_name))
                r_out = os.path.join(cutadapt_file_out_dir, '%s_{name1}-{name2}.2.fastq' % (pair_name))

                log_file = os.path.join(cutadapt_file_out_dir, '%s.log' % (pair_name))

                command_parts = ['cutadapt', '--no-indels', '-g', 'file:%s' % (barfile_f_name), '-q', '30', '-m', '50',
                                 '-G', 'file:%s' % (barfile_r_name), '-e', '0', '-o', f_out, '-p', r_out, '--cores=4',
                                 '--discard-untrimmed', '--max-n', '0', '--pair-filter=any', read1_path, read2_path,
                                 '>', log_file]
                command = ' '.join(command_parts)
                print(command)
                os.system(command)

            # let sample objects know where to find fastq files
            for barcode_pair in self.file_barcode_mapping[file_pair]:
                sample_name = self.file_barcode_mapping[file_pair][barcode_pair]
                sample = self.samples[sample_name]
                sample.demuxed_f_reads = os.path.join(cutadapt_file_out_dir, '%s_%s-%s.1.fastq' % (
                    pair_name, barcode_pair[0], barcode_pair[1]))
                sample.demuxed_r_reads = os.path.join(cutadapt_file_out_dir, '%s_%s-%s.2.fastq' % (
                    pair_name, barcode_pair[0], barcode_pair[1]))

    def map_reads(self):
        minigene_utils.make_dir(self.star_index_folder)
        minigene_utils.make_dir(self.mapped_folder)
        for sample in self.samples.values():
            sample.generate_STAR_index(self.star_index_folder)
            sample.map_reads_star(self.mapped_folder)

    def assemble_transcripts(self):

        for sample in self.samples.values():
            sample.assemble_transcripts_stringtie(self.assembled_folder)
        for ref_name in self.reference_names:
            merged_output = os.path.join(self.assembled_folder, ref_name + '_merged.gtf')
            self.merged_stringtie_gtfs[ref_name] = merged_output
            #for each reference, assemble the list of all of the stringtie output GTFs
            samples = self.reference_to_samples[ref_name]
            stringtie_gtfs = [sample.stringtie_out for sample in samples]
            with open('stringtie_merge_temp.txt', 'w') as f:
                for gtf in stringtie_gtfs:
                    f.write(f'{gtf}\n')
            reference_gtf = self.reference_gtf_files[ref_name]
            #merge all gtfs
            command_parts = ['stringtie', '--merge', '-p 4' ,f'-G {reference_gtf}',f'-o {merged_output}',
                             'stringtie_merge_temp.txt']
            command = ' '.join(command_parts)
            print(command)
            os.system(command)
            '''
            for sample in samples:
                sample.quantify_transcripts_stringtie(self.assembled_folder, merged_output)
            #now assemble the read counts for all of the samples for this reference into a single file using the
            #prepDE.py script from JHU
            #first assemble the gtf output for each sample into a single text file
            with open('prepDE_samples_temp.txt', 'w') as f:
                for sample in samples:
                    f.write(f'{sample.sample_name}\t{sample.stringtie_quant}\n')
            command_parts = ['python', 'prepDE.py3', '-i', 'prepDE_samples_temp.txt' ,f'-G {reference_gtf}',f'-o {merged_output}',
                             'stringtie_merge_temp.txt']
            command = ' '.join(command_parts)
            print(command)
            os.system(command)
            '''


        #for sample in self.samples.values():
        #    sample.parse_stringtie_output()

    def map_reads_for_quant(self):
        minigene_utils.make_dir(self.quant_folder)
        for reference in self.merged_stringtie_gtfs:
            ref_gtf = self.merged_stringtie_gtfs[reference]
            for sample in self.reference_to_samples[reference]:
                sample.quant_reads_star(self.quant_folder, ref_gtf)

    def make_reference_gtfs(self):
        for reference_seq_name in self.target_info:
            # make a gtf file for the expected transcripts to focus on
            minigene_utils.make_dir(self.assembled_folder)
            ref_gtf = os.path.join(self.assembled_folder, reference_seq_name + '_ref.gtf')
            self.reference_gtf_files[reference_seq_name] = ref_gtf
            with open(ref_gtf, 'w') as f:
                for transcript in self.target_info[reference_seq_name]:
                    exons = self.target_info[reference_seq_name][transcript]
                    tx_start = min([exon[0] for exon in exons])
                    tx_end = max([exon[1] for exon in exons])
                    tx_line = '\t'.join(
                        [reference_seq_name, 'expected', 'transcript', str(tx_start), str(tx_end), '.', '+', '.',
                         f'gene_id "{reference_seq_name}"; transcript_id "{transcript}";']) + '\n'
                    f.write(tx_line)
                    for exon in exons:
                        exon_line = '\t'.join([reference_seq_name, 'expected', 'exon', str(exon[0]), str(exon[1]), '.',
                                               '+', '.', f'gene_id "{reference_seq_name}"; transcript_id "{transcript}";']) + \
                                    '\n'
                        f.write(exon_line)

    def make_tables(self, out_prefix):
        count_table_dict = defaultdict(lambda: defaultdict(list))
        percent_table_dict = defaultdict(lambda: defaultdict(list))
        for reference in self.reference_to_samples:
            tx_lists = [sample.tx_counts.keys() for sample in self.reference_to_samples[reference]]
            all_tx =  sorted(set([item for sublist in tx_lists for item in sublist]))
            transcript_order=None
            for sample in self.reference_to_samples[reference]:
                count_table_dict[sample.ref_seq_name]['sample name'].append(sample.sample_name)
                percent_table_dict[sample.ref_seq_name]['sample name'].append(sample.sample_name)
                total_counts = sum(sample.tx_counts.values())
                if transcript_order is None:
                    # i'm doing this convoluted thing to make sure iget a count  of zero for any tx not in the first
                    # sample instead of it just being excluded
                    all_tx_counts = [sample.tx_counts[tx] for tx in all_tx]
                    transcript_order = [key for key, value in sorted(zip(all_tx, all_tx_counts), key=lambda x:x[1],
                                                                    reverse=True)]
                for transcript in transcript_order:
                    count_table_dict[sample.ref_seq_name][transcript].append(sample.tx_counts[transcript])
                    if total_counts>0:
                        percent_table_dict[sample.ref_seq_name][transcript].append(100*sample.tx_counts[transcript] /
                                                                               total_counts)
                    else:
                        percent_table_dict[sample.ref_seq_name][transcript].append(0)

        for target in count_table_dict.keys():
            count_table = pd.DataFrame(count_table_dict[target])
            count_table.to_csv('%s%s_splicing_analysis_counts.tsv' % (out_prefix, target), sep='\t')
            percent_table = pd.DataFrame(percent_table_dict[target])
            percent_table.to_csv('%s%s_splicing_analysis_percent.tsv' % (out_prefix, target), sep='\t')
            self.stacked_bar(percent_table, f'{out_prefix}bar.pdf', transcript_order)
            gtf_file = self.merged_stringtie_gtfs[target]
            self.draw_isoforms(gtf_file, f'{out_prefix}transcripts.pdf')

    def stacked_bar(self, summary_df, out_name, order):
        summary_df = summary_df.sort_values(['sample name'])
        ind = np.arange(len(summary_df))  # the x locations for the groups
        width = 0.8  # the width of the bars: can also be len(x) sequence
        plotLayers = []
        bottoms = [0] * len(summary_df)
        bottoms = np.array(bottoms)

        fig = plt.figure(figsize=(len(summary_df)/2, 4))
        num_plots_wide = 1
        num_plots_high = 1
        plot = fig.add_subplot(num_plots_high, num_plots_wide, 1)

        color_index = 0
        for anno_type in order:
            # print(summary_df[anno_type].values)
            plotLayers.append(
                plot.bar(ind, summary_df[anno_type].values, width, bottom=bottoms,
                         color=minigene_utils.colors[color_index],
                         label=anno_type, hatch=None))
            color_index += 1
            bottoms = bottoms + np.array(summary_df[anno_type].values)

        # plt.title('summary of read loss in pipeline')
        plot.set_xticks(ind)
        plot.set_xticklabels(summary_df['sample name'].values, rotation=85)
        plt.tight_layout()
        # Shrink current axis by 40%
        box = plot.get_position()
        plot.set_position([box.x0, box.y0, box.width * 0.6, box.height])
        # Put a legend to the right of the current axis
        plot.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plot.set_ylim(0, 100)
        plot.spines['right'].set_visible(False)
        plot.spines['top'].set_visible(False)
        # plt.subplots_adjust(bottom=0.38, right=0.8)
        plt.savefig(out_name, transparent=True)

    def draw_isoforms(self, gtf_file, out_file):
        plt.rcParams['pdf.fonttype'] = 42  # leaves most text as actual text in PDFs, not outlines
        plt.rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})

        transcript_info = minigene_utils.gtf_to_dict(gtf_file)

        num_plots_wide = 1
        num_plots_high = len(transcript_info)
        fig, plots = plt.subplots(num_plots_high, figsize=(10 * num_plots_wide, 1 * num_plots_high), sharex=True)

        positions = np.arange(1447)
        plot_index = 0
        for transcript in transcript_info:
            transcript_bottom = 35
            CDS_width = 30
            plot=plots[plot_index]
            plot.plot(transcript_info[transcript]['span'],
                      [transcript_bottom + CDS_width/2, transcript_bottom + CDS_width/2],
                      color=minigene_utils.black)
            exons = transcript_info[transcript]['exons']
            for exon in exons:
                plot.add_patch(patches.Rectangle((exon[0], transcript_bottom), 1+exon[1]-exon[0], CDS_width,
                                                 facecolor="black", edgecolor="black", zorder=10))

            plot.spines["bottom"].set_visible(False)
            plot.spines["left"].set_visible(False)
            plot.spines["top"].set_visible(False)
            plot.spines["right"].set_visible(False)
            plot.set_title(transcript)
            plot.set_yticks([])
            plot.set_ylim(0, 100)
            plot_index += 1
        #for plot in plots[:-1]:
        #    plot.set_xticks([])
        plt.savefig(out_file, transparent=True)

class sample:
    def __init__(self, seq_run, sample_name, f_reads_original, r_reads_original, ref_seq_name, f_bar, r_bar):
        """
        Constructor for sample class which will hold information for a single sample after read demultiplexing
        """
        self.seq_run = seq_run
        self.sample_name = sample_name
        self.f_reads_original = f_reads_original
        self.r_reads_original = r_reads_original
        self.ref_seq_name = ref_seq_name
        self.ref_seq = self.seq_run.target_seqs[self.ref_seq_name]
        self.ref_gtf = self.seq_run.reference_gtf_files[self.ref_seq_name]
        self.f_bar = f_bar
        self.r_bar = r_bar
        self.demuxed_f_reads = None
        self.demuxed_r_reads = None

        # self.wt_alleles = [self.wt_seq[i] for i in self.edit_positions]
        # self.edit_alleles = [self.edited_seq[i] for i in self.edit_positions]

        #self.deduplicated_reads = defaultdict(int)  # map tuple of read pairs to counts, not currently used

        self.window_read_seqs = defaultdict(str)
        self.window_read_seq_counts = defaultdict(int)

        self.categorized_reads = defaultdict(dict)
        self.categorized_read_counts = defaultdict(int)


    '''
    def collapse_identical_reads(self):
        line_counter = 0
        for read_f_line, read_r_line in zip(gzip.open(self.demuxed_f_reads, mode='rt'),
                                            gzip.open(self.demuxed_r_reads, mode='rt')):
            if line_counter % 4 == 1:
                # print(read_f_line.strip(), read_r_line.strip())
                # may want to trim for read quality!!!!
                fwd_read = read_f_line.strip()
                rev_read = read_r_line.strip()
                self.deduplicated_reads[(fwd_read, rev_read)] += 1
            line_counter += 1
    '''
    def generate_STAR_index(self, index_folder):
        # build STAR reference index for WT seq
        sa_index = math.log2(len(self.ref_seq)) // 2 - 1
        self.genomeDir = os.path.join(index_folder, self.ref_seq_name)
        if not minigene_utils.file_exists(self.genomeDir):
            with open('temp.fa', 'w') as temp_fasta:
                temp_fasta.write(f'>{self.ref_seq_name}\n')
                temp_fasta.write(f'{self.ref_seq}\n')
            os.mkdir(self.genomeDir)
            command_parts = ['STAR', '--runThreadN', '%d' % self.seq_run.threads, '--runMode genomeGenerate',
                             f'--genomeDir {self.genomeDir}', '--genomeFastaFiles temp.fa',
                             f'--genomeSAindexNbases {sa_index}']
            command = ' '.join(command_parts)
            print(command)
            os.system(command)

    def map_reads_star(self, mapped_folder, peOverlapNbasesMin=20, peOverlapMMp = .01):
        '''
        STAR --genomeDir /home/sanderson/largevolume/zcpan/Boris/GRCh37_star_gencode --alignEndsType EndToEnd  --sjdbGTFfile /home/sanderson/largevolume/zcpan/Boris/GRCh37_star_gencode/gencode.v19.annotation.gtf --outSAMtype BAM SortedByCoordinate --twopassMode Basic --readFilesCommand zcat --readFilesIn reads_1.fq.gz reads_2.fq.gz --runThreadN 4  --outFileNamePrefix file_name_prefix
        but zcat doesn't right work on osx
        '''

        self.mapping_log = os.path.join(mapped_folder, self.sample_name + '.log')
        self.mapped_reads_prefix = os.path.join(mapped_folder, self.sample_name)
        self.sorted_bam = os.path.join(mapped_folder, self.sample_name + 'Aligned.sortedByCoord.out.bam')
        self.umapped_mate1 = os.path.join(mapped_folder, self.sample_name + 'Unmapped.out.mate1')
        self.umapped_mate2 = os.path.join(mapped_folder, self.sample_name + 'Unmapped.out.mate2')
        self.umapped_summary = os.path.join(mapped_folder, self.sample_name + '_abundant_unmapped.tsv')
        self.splice_junctions = os.path.join(mapped_folder, self.sample_name + 'SJ.out.tab')
        #f'--sjdbGTFfile {self.ref_gtf}',
        if not minigene_utils.file_exists(self.sorted_bam):
            command_parts = ['STAR', '--genomeDir', self.genomeDir, '--alignEndsType EndToEnd',
                             '--outSAMtype BAM SortedByCoordinate',
                             '--twopassMode Basic', '--outReadsUnmapped Fastx', '--limitBAMsortRAM 10000000000',
                             f'--peOverlapNbasesMin {peOverlapNbasesMin}'
                             f'--peOverlapMMp {peOverlapMMp}'
                             f'--alignEndsProtrude 30 ConcordantPair',
                             '--runThreadN %d' % (self.seq_run.threads),
                             '--readFilesIn %s %s' % (self.demuxed_f_reads, self.demuxed_r_reads),
                             '--outFileNamePrefix', self.mapped_reads_prefix,
                             '1>>', self.mapping_log, '2>>', self.mapping_log]
            command = ' '.join(command_parts)
            print(command)
            os.system(command)
            command_parts = ['samtools', 'index', self.sorted_bam]
            command = ' '.join(command_parts)
            print(command)
            os.system(command)
        minigene_utils.summarize_abundant_read_pairs(self.umapped_mate1, self.umapped_mate2, self.umapped_summary,
                                                     output_number=100)
        # self.mapped_pairs = functools.reduce(lambda x, y: x + y, [ int(l.rstrip('\n').split('\t')[2]) for l in pysam.idxstats(pysam.AlignmentFile(self.sorted_bam, "rb")) ])
    def quant_reads_star(self, quant_folder, merged_GTF, peOverlapNbasesMin=20, peOverlapMMp = .01, ):
        '''
        STAR --genomeDir /home/sanderson/largevolume/zcpan/Boris/GRCh37_star_gencode --alignEndsType EndToEnd  --sjdbGTFfile /home/sanderson/largevolume/zcpan/Boris/GRCh37_star_gencode/gencode.v19.annotation.gtf --outSAMtype BAM SortedByCoordinate --twopassMode Basic --readFilesCommand zcat --readFilesIn reads_1.fq.gz reads_2.fq.gz --runThreadN 4  --outFileNamePrefix file_name_prefix
        but zcat doesn't right work on osx
        '''

        self.quant_mapping_log = os.path.join(quant_folder, self.sample_name + '.log')
        self.quant_mapped_reads_prefix = os.path.join(quant_folder, self.sample_name)
        self.quant_bam = os.path.join(quant_folder, self.sample_name + 'Aligned.toTranscriptome.out.bam')
        self.quant_sorted_bam = os.path.join(quant_folder, self.sample_name + 'Aligned.toTranscriptome.out.sorted.bam')
        self.quant_umapped_mate1 = os.path.join(quant_folder, self.sample_name + 'Unmapped.out.mate1')
        self.quant_umapped_mate2 = os.path.join(quant_folder, self.sample_name + 'Unmapped.out.mate2')
        self.quant_umapped_summary = os.path.join(quant_folder, self.sample_name + '_abundant_unmapped.tsv')
        self.quant_splice_junctions = os.path.join(quant_folder, self.sample_name + 'SJ.out.tab')
        if not minigene_utils.file_exists(self.quant_sorted_bam):
            command_parts = ['STAR', '--genomeDir', self.genomeDir, '--alignEndsType EndToEnd',
                             '--outSAMtype BAM SortedByCoordinate',
                             '--quantMode TranscriptomeSAM GeneCounts',
                             f'--sjdbGTFfile {merged_GTF}',
                             f'--alignEndsProtrude 30 ConcordantPair',
                             '--outReadsUnmapped Fastx', '--limitBAMsortRAM 10000000000',
                             f'--peOverlapNbasesMin {peOverlapNbasesMin}'
                             f'--peOverlapMMp {peOverlapMMp}'
                             '--runThreadN %d' % (self.seq_run.threads),
                             '--readFilesIn %s %s' % (self.demuxed_f_reads, self.demuxed_r_reads),
                             '--outFileNamePrefix', self.quant_mapped_reads_prefix,
                             '1>>', self.quant_mapping_log, '2>>', self.quant_mapping_log]
            command = ' '.join(command_parts)
            print(command)
            os.system(command)
            command_parts = ['samtools', 'sort', f'-o {self.quant_sorted_bam}', '-m 1G', f'{self.quant_bam}']
            command = ' '.join(command_parts)
            print(command)
            os.system(command)
            command_parts = ['samtools', 'index', self.quant_sorted_bam]
            command = ' '.join(command_parts)
            print(command)
            os.system(command)
        minigene_utils.summarize_abundant_read_pairs(self.umapped_mate1, self.umapped_mate2, self.umapped_summary,
                                                     output_number=100)
        self.tx_counts = minigene_utils.count_transcript_reads_from_bam(self.quant_sorted_bam)
        # self.mapped_pairs = functools.reduce(lambda x, y: x + y, [ int(l.rstrip('\n').split('\t')[2]) for l in pysam.idxstats(pysam.AlignmentFile(self.sorted_bam, "rb")) ])

    def assemble_transcripts_stringtie(self, assembled_folder):
        """
        stringtie [-o <output.gtf>] [other_options] <read_alignments.bam>
        """
        self.assembly_log = os.path.join(assembled_folder, self.sample_name + '.log')
        self.stringtie_out = os.path.join(assembled_folder, self.sample_name + '.gtf')
        if not minigene_utils.file_exists(self.stringtie_out):
            command_parts = ['stringtie', '-o', self.stringtie_out, '--fr', '-m 50', f'-G {self.ref_gtf}',
                             self.sorted_bam]
            command = ' '.join(command_parts)
            print(command)
            os.system(command)

    def quantify_transcripts_stringtie(self, assembled_folder, merged_gtf):
        """
        stringtie [-o <output.gtf>] [other_options] <read_alignments.bam>
        """
        self.stringtie_quant = os.path.join(assembled_folder, self.sample_name + '_quant.gtf')
        self.ctab_folder = os.path.join(assembled_folder, self.sample_name)
        minigene_utils.make_dir(self.ctab_folder)
        if not minigene_utils.file_exists(self.stringtie_quant):
            command_parts = ['stringtie', '-o', self.stringtie_quant, f'-b {self.ctab_folder}', '-e', '--fr', '-m 50', f'-G {merged_gtf}',
                             self.sorted_bam]
            command = ' '.join(command_parts)
            print(command)
            os.system(command)



    def parse_stringtie_output(self):
        self.assembled_transcript_FPKMs = {}
        self.assembled_transcript_lengths = defaultdict(int)
        f = open(self.stringtie_out, 'r')
        for line in f:
            if line.startswith('#'):
                pass
            else:
                seqname, source, feature, start, end, score, strand, frame, attributes = line.strip().split('\t')
                att_dict = {}
                split_att = attributes.strip(';').split(';')
                split_att = [att.strip() for att in split_att]  # remove whitespace
                for att_pair in split_att:
                    att, value = att_pair.split(' ')
                    value = value.strip('"')  # remove quotes
                    att_dict[att] = value
                if feature == 'transcript':
                    if 'reference_id' in att_dict:
                        self.assembled_transcript_FPKMs[att_dict['reference_id']] = float(att_dict['FPKM'])
                    else:
                        self.assembled_transcript_FPKMs[att_dict['transcript_id']] = float(att_dict['FPKM'])
                if feature == 'exon':
                    exon_length = 1 + int(end) - int(start)
                    if 'reference_id' in att_dict:
                        self.assembled_transcript_lengths[att_dict['reference_id']] += exon_length
                    else:
                        self.assembled_transcript_lengths[att_dict['transcript_id']] += exon_length
        f.close()
        self.assembled_transcript_rpms = {}
        for tx in self.assembled_transcript_FPKMs:
            self.assembled_transcript_rpms[tx] = self.assembled_transcript_FPKMs[tx] * (
                        self.assembled_transcript_lengths[tx] / 1000.)