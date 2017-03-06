#!/usr/bin/env python

from collections import Counter
import copy
import gzip
import itertools
import matplotlib.pyplot as plt
import numpy as np
from operator import itemgetter
import os
import os.path
import sys

LETTERS = ('A', 'T', 'G', 'C')

def iterFasta(path):
    """Iterates over the sequences of a fasta file"""
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
            if name: # Then first seq has been read and it is time to yield
                yield (name, ''.join(seq))
            # Update the name to current seq and blank the sequence
            name = line.strip('>').split()[0]
            seq = []
        else:
            # Extend the current sequence with line
            seq.append(line)
    if name:
        # Yield the final sequence at end of file.
        yield (name, ''.join(seq))
    handle.close()

def countN(sequence):
    """Count number of ATGC and non-ATGC bases in sequence string"""
    count  		= 0
    countN 		= 0
    baseTally 	= Counter()
    baseTally += Counter(sequence)
    for i in baseTally:
        if i not in LETTERS:
            countN 	+= baseTally[i]
        else:
            count 	+= baseTally[i]
    # Return tuple with readable base count and non-ATGC base count
    return (count, countN)

''' Make blank dict '''
def prepareMaps(k, maxk, kmers):
    """Prepares the kmer maps for the specified kmer length"""
    if k == maxk:
        kmer2int = {}
        for kmer in kmers:
            kmer2int[kmer] = 0
        return kmer2int
    newKmers = []
    for kmer in kmers:
        for letter in LETTERS:
            newKmers.append(kmer + letter)
    kmers = newKmers
    return prepareMaps(k + 1, maxk, kmers)

''' counts from window / genome '''
def computeKmers(genomefasta=None, genomepickle=None, window=None, genomeMode=False, pcaMode=False, kmerMap=None, getMeta=True, sym=False):
    """	Compute the kmer counts throughout the kmer range for a sequence or collection of sequences,
        optionally write the output to file."""
    # args					'arguments object'
    # genomepickle=None		'If run in genome mode dump kmer map to this file'
    # window=None			'tuple(windowID,Seq)'
    # genomeMode=None		'If genome mode: load genome fasta from args.hostSeq, count kmers symmetrically'
    # kmerMap=None			'Blank kmer map - list of dictionaries'
    # getMeta=True			'Store metadata in final output: SeqLen, excluded maxmers, N's in seq'
    # sym=False				'Calculate kmers symmetrically - as double stranded seq'

    if sum(kmerMap[0].itervalues()) != 0:
        logging.info('kmer template is not blank!')
        sys.exit(1)

    if genomeMode:

        targetSeq = iterFasta(genomefasta)
    else:
        targetSeq = window

    excludedMaxMer = dict()
    excludedMaxMer['exMax'] = 0
    surveySeqLen = dict()
    surveySeqLen['totalLen'] = 0
    nCounter 	= dict()
    nCounter['nnTotal'] = 0

    # Prepare all maps
    if pcaMode:
        kMin    = args.pcaMin
        kMax    = args.pcaMax
    else:
        kMin    = args.minWordSize
        kMax    = args.maxWordSize

    # Redefine the blank kmer dict being to avoid edits to the original
    maps    = copy.deepcopy(kmerMap)
    # Iterate over sequences
    nSeqs   = 0
    # Read (name, sequence) tuple
    for name, seq in targetSeq:
        size   = len(seq)
        surveySeqLen['totalLen'] += size
        baseCount, nnCount 	 = countN(seq)
        nCounter['nnTotal']	+=	nnCount
        # For each kmer length
        for i in xrange(kMin, kMax + 1):
            # Crawl sequence one base at a time (stop k bases from end)
            for j in xrange(size - i + 1):
                # Extract word of len i (kmer) starting at current base
                word = seq[j:j + i]
                # Do not mask query windows
                # Optionally include lowercase masking from host genome so that kmer returns idx = None.
                if not genomeMode:
                    word = word.upper()
                elif not args.maskHost:
                    word = word.upper()
                # This checks that the string is a legit kmer in the map i.e. does not contain Ns
                # .get will return 'None' if kmer key does not exist, or value if it does.
                # Updates idx original 0 count of word
                idx  = maps[i - kMin].get(word, None)
                # If the word contains characters other than ATGC, skip to next
                if idx is None:
                    if i == kMax:
                        excludedMaxMer['exMax'] += 1
                    continue
                # Increment instances of current word
                maps[i - kMin][word] += 1
                # Increment count for reverse complement of current word
                if genomeMode or sym:
                    maps[i - kMin][revComplement(word)] += 1

        # Iterate to next sequence in input fasta
        nSeqs += 1

    if getMeta:
        maps.append(surveySeqLen)
        maps.append(excludedMaxMer)
        maps.append(nCounter)

    if genomeMode:
        # Save genome kmers to pickle
        pickle.dump(maps, open(genomepickle, "wb"))
        # Trace
        logging.info('Processed %d sequences' % nSeqs)
    # Return a dictionary of kmer counts
    return maps

''' IVOM from window '''
def IvomBuild(querySeq, maxKlen=8, minKlen=1, GenomeKmers, isGenomeIVOM):
    ''' Input is a dictionary of all kmers in the current window.
        Should include a final dict with N count for each kmer length.
        Calculates weights Wi(=Counts*deg_freedom) and obs_freqs (Pi)'''
    
    windowKmers = computeKmers(querySeq)
    windowlen = len(querySeq)
    kRange      = maxKlen - minKlen # Equal to index position of maxmer

    # Note: Alternative define total genome space as total readable number of maxmers.
    genomeSpace = GenomeKmers[kRange+1]['totalLen'] - GenomeKmers[kRange+3]['nnTotal']
    windowSpace = windowKmers[kRange+1]['totalLen'] - windowKmers[kRange+3]['nnTotal']

    sumWindowIVOM = 0
    subKmers = dict()
    storeMaxMerIVOM = dict()
    # For each max len kmer
    for k in windowKmers[kRange]:
        # Skip kmers that do not occur in window - this could be more efficient.
        if windowKmers[kRange][k] == 0:
            continue
        count = windowKmers[kRange][k]
        subKmers['w'] = dict()
        subKmers['p'] = dict()

        if not isGenomeIVOM: # Calculating IVOM weights for window kmers against window counts
            for x in range(minKlen, maxKlen+1):
                # Process maxmer
                if x == maxKlen:
                    # Note:Weighting could probably be achieved by observations multiplied by kmerlen squared i.e. count * x**2
                    subKmers['w'][x] = count * 4**x # kmer_count * (DNA bases^seqlength)
                    # Note: Should probability of observing x count of k be = x /total available non-N kmers in window?
                    subKmers['p'][x] = float(count) / ((windowSpace-(x-1)) * 2)
                elif x >= 2:
                    subK = k[0:x] # Grab the first x bases in maxmer
                    subKmers['w'][x] = windowKmers[x-minKlen][subK] * 4**x
                    subKmers['p'][x] = float(windowKmers[x-minKlen][subK]) / ((windowSpace - (x-1)) * 2)
                else:
                    subK = k[0] # Grab the first base in maxmer
                    subKmers['w'][x] = windowKmers[x-minKlen][subK] * 4**x
                    subKmers['p'][x] = float(windowKmers[x-minKlen][subK]) / ((windowSpace) * 2)

        else: # Calculating IVOM weights for current window kmer against genome counts
            for x in range(minKlen, maxKlen+1):
                # Process maxmer
                if x == maxKlen:
                    subKmers['w'][x] = GenomeKmers[x-minKlen][k] * 4**x
                    subKmers['p'][x] = float(GenomeKmers[x-minKlen][k]) / ((genomeSpace - (x-1)) * 2)
                elif x >= 2:
                    subK = k[0:x] # Grab the first x bases in maxmer
                    subKmers['w'][x] = GenomeKmers[x-minKlen][subK] * 4**x
                    subKmers['p'][x] = float(GenomeKmers[x-minKlen][subK]) / ((genomeSpace - (x-1)) * 2)
                else:
                    subK = k[0] # Grab the first base in maxmer
                    subKmers['w'][x] = GenomeKmers[x-minKlen][subK] * 4**x
                    subKmers['p'][x] = float(GenomeKmers[x-minKlen][subK]) / (genomeSpace * 2)

        w_total = 0
        w_running_totals = dict()
        for x in subKmers['w']:
            # Running total of sub-kmer raw weights
            w_total += subKmers['w'][x]
            # Add total at eact position to dict
            w_running_totals[x] = w_total

        scaledW = dict()
        for x in range(minKlen, maxKlen+1):
            # Note: Weighting of sub-kmers not explictily addressed in Vernikos and Parkhill.
            scaledW[x] = float(subKmers['w'][x]) / w_running_totals[x]

        kmerIVOM = dict()
        for x in range(minKlen, maxKlen+1):
            if x == minKlen:
                kmerIVOM[x] = scaledW[x] * subKmers['p'][x]
            elif x < maxKlen:
                kmerIVOM[x] = scaledW[x] * subKmers['p'][x] + ((1-scaledW[x]) * kmerIVOM[x-1])
            else: # IVOM for max len kmer
                kmerIVOM[x] = scaledW[x] * subKmers['p'][x] + ((1-scaledW[x]) * kmerIVOM[x-1])

        storeMaxMerIVOM[k] = kmerIVOM[maxKlen]
        # Note: In Genome mode 'window' is eqivalent to whole readable genome space.
        sumWindowIVOM += kmerIVOM[maxKlen]

    # Rescale each kmer from 0 to 1 (because relative entropy)
    for k in storeMaxMerIVOM:
        storeMaxMerIVOM[k] = float(storeMaxMerIVOM[k]) / sumWindowIVOM

    # Return dictionary of scaled IVOM scores keyed by kmer
    return storeMaxMerIVOM

''' KLD '''
def KLD(GenomeIVOM, windowIVOM):
    """	Kullback-Leibler Divergence: Measure of relative entropy between sets of
        probability distributions."""
    # IVOM structure: dict[kmer][IVOM Score]
    # Negative number indicates kmer depleted in window relative to genome
    # Positive number indicates enriched in window relative to genome
    windowKLD = 0
    for k in windowIVOM:
        w = float(windowIVOM[k])
        G = float(GenomeIVOM[k])
        if G != 0:
            windowKLD += (w*math.log((w/G), 2))
        # Note: Using log2 instead of log10 to accentuate variance within small range.
    return windowKLD



#Example usage

genomefasta = 'path/to/test.fa.gz'
GenomeKmers = computeKmers(genomefasta,gmode=True)


#calc genome kmers


