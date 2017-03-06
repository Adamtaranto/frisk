#!/usr/bin/env python

from matplotlib import pyplot as plt
from matplotlib.collections import BrokenBarHCollection
import pandas as pd
import pybedtools
import logging
import itertools
from operator import itemgetter

logging.basicConfig(level=logging.INFO, format=("%(asctime)s - %(funcName)s - %(message)s"))

# Here's the function that we'll call for each dataframe (once for chromosome
# ideograms, once for genes).  The rest of this script will be prepping data
# for input to this function
#
def chromosome_collections(df, y_positions, height,  **kwargs):
    """
    Yields BrokenBarHCollection of features that can be added to an Axes
    object.
    Parameters
    ----------
    df : pd.DataFrame
        Must at least have columns ['chrom', 'start', 'end', 'color']. If no
        column 'width', it will be calculated from start/end.
    y_positions : dict
        Keys are chromosomes, values are y-value at which to anchor the
        BrokenBarHCollection
    height : float
        Height of each BrokenBarHCollection
    Additional kwargs are passed to BrokenBarHCollection
    """
    del_width = False
    if 'width' not in df.columns:
        del_width = True
        df['width'] = df['end'] - df['start']
    for chrom, group in df.groupby('chrom'):
        print chrom
        yrange = (y_positions[chrom], height)
        xranges = group[['start', 'width']].values
        yield BrokenBarHCollection(
            xranges, yrange, facecolors=df['colors'], **kwargs)
    if del_width:
        del df['width']

def gffFilter(rec, **kwargs):
    """ Given list of feature types return GFF3 records that match type."""
    if 'feature' in kwargs.keys():
        if rec[2] in kwargs['feature']:
            return rec
    else:
        return rec

def iterFasta(path):
    """Iterates over the sequences of a fasta file"""
    logging.info("Loading fasta files from %s" % path)
    name = None
    seq = []
    if path.endswith('.gz') or path.endswith('.gz"'):
        handle = gzip.open(path)
    else:
        handle = open(path)
    for line in handle:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name: #then first seq has been read and it is time to yield
                yield (name, ''.join(seq))
            # Update the name to current seq and blank the sequence
            name = line.strip('>').split()[0]
            seq = []
        else:
            # Extend the current sequence with uppercase sequence
            # Note: Could break on masked regions to exclude them from global training set.
            seq.append(line.upper())
    if name:
        # Yield the final sequence at end of file.
        yield (name, ''.join(seq))
    handle.close()

def getFasta(fastaPath):
    """Write fasta to dictionary, key by scaffold name."""
    seqDict = dict()
    nRanges = list()
    for name, seq in iterFasta(fastaPath):
        seqDict[name] = seq
        for i in findBaseRanges(seq, 'N', name=name, minlen=10):
            nRanges.append(i)
    nBlocks = pybedtools.BedTool(nRanges)
    return seqDict, nBlocks

def findBaseRanges(s, ch, name=None, minlen=0):
    """For string return intervals of character longer than minimum length."""
    data    = [i for i, ltr in enumerate(s) if ltr == ch]
    ranges  = list()
    for k, g in itertools.groupby(enumerate(data), lambda (i, x): i-x):
        group =  map(itemgetter(1), g)
        if (group[-1] - group[0]) < minlen:
            continue
        else:
            if name:
                ranges.append((name, group[0], group[-1]))
            else:
                ranges.append((group[0], group[-1]))
    return ranges #Format = [('ScaffName,start,stop'),('ScaffName,start,stop')]

# Height of each ideogram
chrom_height = 1

# Spacing between consecutive ideograms
chrom_spacing = 1

# Height of the gene track. Should be smaller than `chrom_spacing` in order to
# fit correctly
gene_height = 0.4

# Padding between the top of a gene track and its corresponding ideogram
gene_padding = 0.1

# Width, height (in inches)
figsize = (6, 8)


selfGenome, nBlocks = getFasta(querySeq)

# Decide which chromosomes to use
chromosome_list = sorted(selfGenome.keys())

chromo_dict = dict()
for chr in chromosome_list:
    chromo_dict[chr] = len(selfGenome[chr])

# Keep track of the y positions for ideograms and genes for each chromosome,
# and the center of each ideogram (which is where we'll put the ytick labels)
ybase = 0
chrom_ybase = {}
gene_ybase = {}
chrom_centers = {}

# Iterate in reverse so that items in the beginning of `chromosome_list` will
# appear at the top of the plot
for chrom in chromosome_list[::-1]:
    chrom_ybase[chrom] = ybase
    chrom_centers[chrom] = ybase + chrom_height / 2.
    gene_ybase[chrom] = ybase - gene_height - gene_padding
    ybase += chrom_height + chrom_spacing

#Blank scaffold
#figure out how to read dict to pandas {'chr1':len}
scaffolds = pd.DataFrame(list(chromo_dict.iteritems()),columns=['chrom','end'])
scaffolds['start'] = 0
# Filter out chromosomes not in our list
scaffolds = scaffolds[scaffolds.chrom.apply(lambda x: x in chromosome_list)]
scaffolds['width'] = scaffolds.end - scaffolds.start
scaffolds['colors'] = '#eeeeee'

##Read in anomalies
testAnom = [['scaffold_001',1000,60000,'text',42],['scaffold_001',80000,99000,'text',42],['scaffold_001',200000,290000,'text',42]]
testBED = pybedtools.BedTool(testAnom)
features = pd.read_table(testBED.fn, names=['chrom', 'start', 'end'],usecols=range(3))
# Filter out chromosomes not in our list
features = features[features.chrom.apply(lambda x: x in chromosome_list)]
features['width'] = features.end - features.start
features['colors'] = '#ff6600'

# Same thing for genes
gffIn = 'data/SN15v6.gff3'
gffFeatures = ['gene']

if type(gffFeatures) is list:
    gffType = gffFeatures[0]
else:
    gffType = gffFeatures

filteredGff = list()

handle = open(gffIn)
for line in handle:
    line = line.strip()
    if not line:
        continue
    if not line.startswith("#"):
        rec = line.split()
        reclist = gffFilter(rec, feature=gffType)
        if not reclist:
            continue
        recInterval = (reclist[0],int(reclist[3]),int(reclist[4]))
        filteredGff.append(recInterval)

handle.close()


genes = pd.DataFrame.from_records(filteredGff, columns=['chrom', 'start', 'end'])
genes = genes[genes.chrom.apply(lambda x: x in chromosome_list)]
genes['width'] = genes.end - genes.start
genes['colors'] = '#2243a8'



fig = plt.figure()
ax = fig.add_subplot(111)

# Now all we have to do is call our function for the ideogram data...

print("adding scaffolds...")
for collection in chromosome_collections(scaffolds, chrom_ybase, chrom_height):
    ax.add_collection(collection)

print("adding anomalies...")
for collection in chromosome_collections(features, chrom_ybase, chrom_height):
    ax.add_collection(collection)

# ...and the gene data
print("adding genes...")
for collection in chromosome_collections(
    genes, gene_ybase, gene_height, alpha=0.5, linewidths=0
):
    ax.add_collection(collection)

# Axes tweaking
ax.set_yticks([chrom_centers[i] for i in chromosome_list])
ax.set_yticklabels(chromosome_list)
ax.axis('tight')
plt.show()