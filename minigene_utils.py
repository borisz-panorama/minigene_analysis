import os.path
import os
import hashlib
import gzip
import matplotlib.pyplot as plt

plt.rcParams['pdf.fonttype'] = 42 #leaves most text as actual text in PDFs, not outlines

black = (0,0,0)
gray = (127/255.0,127/255.0,127/255.0)
orange = (230/255.0,159/255.0,0)
skyBlue = (86/255.0,180/255.0,233/255.0)
bluishGreen = (0,158/255.0,115/255.0)
yellow = (240/255.0,228/255.0,66/255.0)
blue = (0,114/255.0,178/255.0)
vermillion = (213/255.0,94/255.0,0)
reddishPurple = (204/255.0,121/255.0,167/255.0)
colors = [gray, orange, skyBlue, bluishGreen, blue, reddishPurple, vermillion, yellow]

def reverse_complement(seq, isRNA = False):
    seq = seq.upper()
    compDict = {'A':'T', 'T':'A', 'U':'A', 'C':'G', 'G':'C', 'N':'N', '-':'-', '.':'.', '*':'*'}
    revComp = ''.join([compDict[c] for c in seq[::-1]])
    if isRNA:
        return revComp.replace('T', 'U')
    return revComp

def convertFastaToDict(fastaFile):
        '''
        converts a fasta file to a dict of {sequenceName:sequence}
        can take extra files in * args
        '''
        if isinstance(fastaFile, list):
            files = fastaFile
        else:
            files = [fastaFile]
        currentName = None
        currentSequence = None
        seqDict = {}
        for currentFile in files:
            if currentFile.endswith('.gz'):
                f = gzip.open(currentFile)
            else:
                f = open(currentFile)
            for line in f:
                if not line.strip() == '' and not line.startswith('#'):  # ignore empty lines and commented out lines
                    if line.startswith('>'):  # > marks the start of a new sequence
                        if not currentName == None:  # after we've reached the firtst > line, we know what the sequence corresponds to
                            seqDict[currentName] = currentSequence.upper()
                        currentName = line.strip()[1:].split()[0]  # i've noticed the gencode names have extraneous numbering after some whitespace. This doens't match the GTF files, so I'm removing it.
                        currentSequence = ''
                    else:
                        currentSequence += line.strip()
            f.close()
        seqDict[currentName] = currentSequence.upper()
        return seqDict

def get_MD5s_From_file(md5_file):
    file_to_md5 = {}
    f = open(md5_file)
    for line in f:
        md5, file = line.strip().split()
        file_to_md5[file[2:]] = md5
    f.close()
    #print(file_to_md5)
    return file_to_md5

def compute_md5(file_name):
    md5_hash = hashlib.md5()
    #log_file.write('computing md5sum for %s\n' %(file_name))
    #print('computing md5sum for ', file_name)
    with open(file_name,"rb") as f:
        # Read and update hash in chunks of 4K
        for byte_block in iter(lambda: f.read(4096),b""):
            md5_hash.update(byte_block)
        md5 = md5_hash.hexdigest()
    f.close()
    return md5

def make_dir(dirname):
    """
    Makes the directory; doesn't throw an error if it exists.
    """
    if not os.path.exists(dirname):
        try:
            os.makedirs(dirname)
        except:
            print('The directory %s exists' % (dirname))

def file_exists(fname):
    """
    makes sure a given file exists
    """
    if not os.path.exists(fname):
        return False
    fstats = os.stat(fname)
    if not fstats[6]:
        return False
    if not os.access(fname, os.R_OK):
        raise ValueError('Input File %s cannot be read' % fname)
    return True

def get_exon_boundaries(exon_positions):
    '''
    boundaries are of format:
    exon1_start,exon1_end;exon2_start,exon2_end; etc
    as in:
    1,153;254,338;1296;1447
    numbers are 1-based and inclusive on both ends
    '''
    exons = exon_positions.strip().strip('"').split(';')
    exon_boundaries = []
    for exon in exons:
        start, end = [int(x) for x in exon.split(',')]
        exon_boundaries.append((start, end))
    return exon_boundaries