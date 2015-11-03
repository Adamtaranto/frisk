#!/usr/bin/env python
# Copyright (C) 2015 Adam Taranto <adam.p.taranto@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

#########################################################################################
# Detects regions of unusual sequence composition by comparison of local kmer           #
# frequencies to whole genome abundance.                                                #
# For use in detection of genomic islands and segmental lateral gene transfer events.   #
#########################################################################################

from __future__ import print_function, division, absolute_import
import argparse
import array
from collections import Counter
import copy
import datetime
import gzip
from hmmlearn import hmm
import itertools
import logging
import math
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.collections import BrokenBarHCollection
from matplotlib.ticker import FuncFormatter
import numpy as np
from operator import itemgetter
import os
import os.path
import pandas as pd
import pickle
import pybedtools
import re
from sklearn.decomposition import PCA
from sklearn.cluster import DBSCAN
from scipy import stats
import seaborn as sns
import sys
from . import tsne

#######################
#######################
### Global contants ###
#######################
#######################

from ._version import get_versions
__version__ = get_versions()['version']
FRISK_VERSION = get_versions()['version']
del get_versions

LETTERS = ('A', 'T', 'G', 'C')

#######################
#######################
###### Functions ######
#######################
#######################

def tempPathCheck(args):
    """Check if the temporary folder exists, else creates it"""
    logging.info('Checking for temporary folder')
    tempFolder = os.path.abspath(args.tempDir)
    if not os.path.isdir(tempFolder):
        os.makedirs(tempFolder)

def natural_sort(l,key=str): 
    """Sorting of weird scaffold names for chromosome painting."""
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda item: [ convert(c) for c in re.split('([0-9]+)', key(item)) ] 
    return sorted(l, key = alphanum_key)

def findBaseRanges(s, ch, name=None, minlen=0):
    """For string return intervals of character longer than minimum length."""
    data 	= [i for i, ltr in enumerate(s) if ltr == ch]
    ranges 	= list()
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
    # return tuple with readable base count and non-ATGC base count
    return (count, countN)

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
    nRanges = [ (x,str(y),str(z)) for x,y,z in nRanges ]
    nBlocks = pybedtools.BedTool(nRanges)
    return seqDict, nBlocks

def getBEDSeq(fastaDict, BEDintervals):
    """Given BED object with scaffold coordinates,
        fetch sequence from dictionary of sequences, keyed by scaffold name. """
    for rec in BEDintervals:
        if rec[0] in fastaDict:
            logging.info('From %s try get range %s to %s.' % (rec[0], str(int(rec[1])), str(int(rec[2]))))
            seq = fastaDict[rec[0]][int(rec[1])-1:int(rec[2])]
            name = ":".join([rec[0], str(rec[1]), str(rec[2])])
            if len(seq) > 0:
                yield (name, seq)
            else:
                logging.error('Error: Retrieved zero len sequence for %s' % name)
        else:
            logging.error('Scaffold %s not found in reference.' % str(rec[0]))
            continue

def crawlGenome(args, querySeq):
    """Take genome fasta extract scaffolds, extract incremented windows,
        yield (seq, scaffold name, window start, window end)"""
    w 			= args.windowlen
    i 			= args.increment
    saveSmalls 	= args.scaffoldsAll
    nSeqs 		= 0

    # Iterate over scaffolds
    for name, seq in iterFasta(querySeq):
        # Reset stats and counters
        windowCount		= 0
        winExclude 		= 0
        bases, nonbase 	= countN(seq)
        size 			= len(seq)
        jumpback 		= False
        # Process scaffolds that fall below window size + minimum acceptable first increment, as a single window.
        if size <= w + ((w * 0.75) - i) and saveSmalls:
            scaffB, scaffN = countN(seq)
            if scaffN >= 0.3 * size:
                logging.info('%s excluded as > 30 percent unresolved sequence.' % name)
                winExclude += 1
                continue
            else:
                logging.info('Rescuing small scaffold %s' % name)
                yield (seq, name, 1, size)
                windowCount += 1
                nSeqs += 1
        else:
            # Crawl sequence in multiples of increment length until one increment from end of seq
            for j in xrange(0, size - i + 1, i):
                # If window overshoots sequence, jump back one w len from last base.
                if j+w > size:
                    winSeq 		= seq[size - w: size]
                    jumpback 	= True
                else:
                    # Extract window of len w starting at current base
                    winSeq = seq[j: j+w]
                # Screen N content = (baseCount, Ncount)
                winB, winN = countN(winSeq)
                if winN >= 0.3 * len(winSeq):
                    logging.info('Window from %s excluded as > 30 percent unresolved sequence.' % name)
                    winExclude += 1
                    continue
                elif jumpback:
                    yield (winSeq, name, size - w, size)
                else:
                    yield (winSeq, name, j+1, j + w)
                windowCount += 1
        # Iterate to next sequence in input fasta
        logging.info('Extracted %s windows from %s bases in %s' % (windowCount, str(size), name))
        logging.info('Excluded %s windows from %s' % (winExclude, name))
        nSeqs += 1
    logging.info('Successfully processed sequences: %s' % nSeqs)

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

def rangeMaps(kMin, kMax):
    """Calls prepareMaps to write list of kmer dictionaries for range kmin to kmax"""
    maps    = []
    logging.info('Preparing kmer maps')
    # Prepare possible kmer combinations for each len k
    for i in xrange(kMin, kMax + 1):
        maps.append(prepareMaps(0, i, ['']))
    return maps

def revComplement(kmer):
    revcompl = lambda x: ''.join([{'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}[B] for B in x][::-1])
    return revcompl(kmer)

def computeKmers(args, genomepickle=None, window=None, genomeMode=False, kmerMap=None, getMeta=True, sym=False):
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
        logging.info('Computing kmers for %s' % args.hostSeq)
        targetSeq = iterFasta(args.hostSeq)
    else:
        targetSeq = window

    excludedMaxMer = dict()
    excludedMaxMer['exMax'] = 0
    surveySeqLen = dict()
    surveySeqLen['totalLen'] = 0
    nCounter 	= dict()
    nCounter['nnTotal'] = 0

    # Prepare all maps
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

def IvomBuild(windowKmers, args, GenomeKmers, isGenomeIVOM):
    ''' Input is a dictionary of all kmers in the current window.
        Should include a final dict with N count for each kmer length.
        Calculates weights Wi(=Counts*deg_freedom) and obs_freqs (Pi)'''
    windowlen   = args.windowlen
    klen        = args.maxWordSize

    # Note: Alternative define total genome space as total readable number of maxmers.
    genomeSpace = GenomeKmers[klen]['totalLen'] - GenomeKmers[klen+2]['nnTotal']
    windowSpace = windowKmers[klen]['totalLen'] - windowKmers[klen+2]['nnTotal']

    sumWindowIVOM = 0
    subKmers = dict()
    storeMaxMerIVOM = dict()
    # For each max len kmer
    for k in windowKmers[klen-1]:
        # Skip kmers that do not occur in window - this could be more efficient.
        if windowKmers[klen-1][k] == 0:
            continue
        count = windowKmers[klen-1][k]
        subKmers['w'] = dict()
        subKmers['p'] = dict()

        if not isGenomeIVOM: #Calculating IVOM weights for window kmers against window counts
            for x in range(1, klen+1):
                # Process maxmer
                if x == klen:
                    # Note:Weighting could probably be achieved by observations multiplied by kmerlen squared i.e. count * x**2
                    subKmers['w'][x] = count * 4**x # kmer_count * (DNA bases^seqlength)
                    # Note: Should probability of observing x count of k be = x /total available non-N kmers in window?
                    subKmers['p'][x] = float(count) / ((windowSpace-(x-1)) * 2)
                elif x >= 2:
                    subK = k[0:x] #Grab the first x bases in maxmer
                    subKmers['w'][x] = windowKmers[x-1][subK] * 4**x
                    subKmers['p'][x] = float(windowKmers[x-1][subK]) / ((windowSpace - (x-1)) * 2)
                else:
                    subK = k[0] #Grab the first base in maxmer
                    subKmers['w'][x] = windowKmers[x-1][subK] * 4**x
                    subKmers['p'][x] = float(windowKmers[x-1][subK]) / ((windowSpace) * 2)

        else: #Calculating IVOM weights for current window kmer against genome counts
            for x in range(1, klen+1):
                # Process maxmer
                if x == klen:
                    subKmers['w'][x] = GenomeKmers[x-1][k] * 4**x
                    subKmers['p'][x] = float(GenomeKmers[x-1][k]) / ((genomeSpace - (x-1)) * 2)
                elif x >= 2:
                    subK = k[0:x] #Grab the first x bases in maxmer
                    subKmers['w'][x] = GenomeKmers[x-1][subK] * 4**x
                    subKmers['p'][x] = float(GenomeKmers[x-1][subK]) / ((genomeSpace - (x-1)) * 2)
                else:
                    subK = k[0] #Grab the first base in maxmer
                    subKmers['w'][x] = GenomeKmers[x-1][subK] * 4**x
                    subKmers['p'][x] = float(GenomeKmers[x-1][subK]) / (genomeSpace * 2)

        w_total = 0
        w_running_totals = dict()
        for x in subKmers['w']:
            # Running total of sub-kmer raw weights
            w_total += subKmers['w'][x]
            # Add total at eact position to dict
            w_running_totals[x] = w_total

        scaledW = dict()
        for x in range(1, klen+1):
            # Note: Weighting of sub-kmers not explictily addressed in Vernikos and Parkhill.
            scaledW[x] = float(subKmers['w'][x]) / w_running_totals[x]

        kmerIVOM = dict()
        for x in range(1, klen+1):
            if x == 1:
                kmerIVOM[x] = scaledW[x] * subKmers['p'][x]
            elif x < klen:
                kmerIVOM[x] = scaledW[x] * subKmers['p'][x] + ((1-scaledW[x]) * kmerIVOM[x-1])
            else: #IVOM for max len kmer
                kmerIVOM[x] = scaledW[x] * subKmers['p'][x] + ((1-scaledW[x]) * kmerIVOM[x-1])

        storeMaxMerIVOM[k] = kmerIVOM[klen]
        # Note: In Genome mode 'window' is eqivalent to whole readable genome space.
        sumWindowIVOM += kmerIVOM[klen]

    # Rescale each kmer from 0 to 1 (because relative entropy)
    for k in storeMaxMerIVOM:
        storeMaxMerIVOM[k] = float(storeMaxMerIVOM[k]) / sumWindowIVOM

    # Return dictionary of scaled IVOM scores keyed by kmer
    return storeMaxMerIVOM

def KLI(GenomeIVOM, windowIVOM, args):
    """	Kullback-Leibler Index: Measure of relative entropy between sets of
        probability distributions."""
    # IVOM structure: dict[kmer][IVOM Score]
    # Negative number indicates kmer depleted in window relative to genome
    # Positive number indicates enriched in window relative to genome
    windowKLI = 0
    for k in windowIVOM:
        w = float(windowIVOM[k])
        G = float(GenomeIVOM[k])
        if G != 0:
            windowKLI += (w*math.log((w/G), 2))
        # Note: Using log2 instead of log10 to accentuate variance within small range.
    return windowKLI

def calcRIP(windowKmers):
    """	Calculate measures of Repeat Induced Point mutation (RIP);
        a process in fungal genomes which affects CpA --> TpA
        transition mutations. """
    # Product Index
    if windowKmers[1]['AT'] > 0:
        PI = windowKmers[1]['TA'] / float(windowKmers[1]['AT'])
    else:
        PI = None
    # Substrate index
    AC_GT = (windowKmers[1]['AC'] + windowKmers[1]['GT'])
    if AC_GT > 0:
        SI = (windowKmers[1]['CA'] + windowKmers[1]['TG']) / float(AC_GT)
    else:
        SI = None
    # Composite RIP Index
    if PI and SI:
        CRI = PI - SI
    else:
        CRI = None
    return (PI, SI, CRI)

def makePicklePath(args, **kwargs):
    """ Generate pathname for export of genome kmer counts or window KLI scores. """
    pathbase = os.path.basename(args.hostSeq)
    if kwargs['space'] == 'genome':
        pickleOut = os.path.join(args.tempDir, pathbase + "_kmers_" + str(args.minWordSize) + "_" + str(args.maxWordSize) + "_" + kwargs['space'] + '.p')
    elif kwargs['space'] == 'window':
        if args.querySeq:
            pathbase = os.path.basename(args.querySeq)
        pickleOut = os.path.join(args.tempDir, pathbase + "_kmers_" + str(args.minWordSize) + "_" + str(args.maxWordSize) + "_KLI_" + kwargs['space'] + "_" + str(args.windowlen) + "_increment_" + str(args.increment) + '.p')
    return pickleOut

def FDBins(data):
    """	Freedman-Diaconis method for calculating optimal number of bins for dataset."""
    IQR 	= np.subtract(*np.percentile(data, [75, 25])) #Magic '*' unpacks tuple and passes two values to subtract function.
    bins 	= 2 * IQR * math.pow(len(data), (1.0/3.0)) #Number of bins
    bins 	= int(round(bins)) #Otsu needs rounded integer
    return bins

def otsu(data, optBins):
    # data is array of log10(KLI)
    raw 			= data
    data 			= np.atleast_1d(data)
    data 			= data[~ np.isnan(data)]
    data 			= data/(max(abs(data))*-1.0) #Scale to 0-1 and make positive
    hist, binEdges 	= np.histogram(data, bins=optBins)
    hist 			= hist * 1.0 #Covert to floats
    hist_norm 		= hist.ravel()/hist.max() #Normalise hist to largest bin
    Q 				= hist_norm.cumsum()
    bins 			= np.arange(optBins)
    fn_min 			= np.inf
    thresh 			= -1
    for i in xrange(1, optBins): #Start from second position (1) to compare all to bin 0
        p1, p2 = np.hsplit(hist_norm, [i]) # Split normalised hist values into two brackets at bin i (bin 1 < i, bin 2 >=i)
        q1, q2 = Q[i-1], Q[optBins-1] - Q[i-1] # cum sum of bin values #!! yields zero on final bin
        b1, b2 = np.hsplit(bins, [i]) # Split bins into 2 brackets at bin i
        # Finding means and variances
        m1, m2 = q1/len(p1), q2/len(p2)
        v1, v2 = np.sum(np.square(p1-m1))/len(p1), np.sum(np.square(p2-m2))/len(p2)
        # Calculates the minimization function
        fn = (v1*q1) + (v2*q2)
        if fn < fn_min:
            fn_min = fn
            thresh = i
    logging.info('OTSU selected bin %s as threshold position' % str(thresh))
    # Convert bact to log10(KLI)
    otsuNum = binEdges[thresh] * max(abs(raw)) * -1.0
    return otsuNum

def gffFilter(rec, **kwargs):
    """	Given list of feature types return GFF3 records that match type."""
    if 'feature' in kwargs.keys():
        if rec[2] in kwargs['feature']:
            return rec
    else:
        return rec

def anomaly2GFF(anomBED, **kwargs):
    n = 1
    if 'category' in kwargs.keys():
        annotType = kwargs['category']
    else:
        annotType = 'Kmer-anomaly'
    for i in anomBED:
        content = [str(i[0]), 'frisk_' + FRISK_VERSION, annotType, str(i[1]), str(i[2]), '.', '+', '.', ';'.join(['ID=Anomaly_' + str(n).zfill(len(str(len(anomBED)))), 'maxKLI=' + str(i[3]), 'minKLI=' + str(i[4]), 'meanKLI=' + str(i[5])])]
        if n == 1:
            yield '##gff-version 3' + '\n'
        yield '\t'.join(content) + '\n'
        n += 1

def RIP2GFF(ripBED):
    # name, start, stop, KLI max, PI min, SI max, CRI min, CRI max
    n = 1
    annotType = 'RIP'
    ripBED = natural_sort(ripBED, key=itemgetter(0))
    for i in ripBED:
        content = [str(i[0]), 'frisk_' + FRISK_VERSION, annotType, str(i[1]), str(i[2]), '.', '+', '.', ';'.join(['ID=Anomaly_' + str(n).zfill(len(str(len(ripBED)))), 'maxKLI=' + str(i[3]), 'minPI=' + str(i[4]), 'maxSI=' + str(i[5]), 'minCRI=' + str(i[6]), 'maxCRI=' + str(i[7])])]
        if n == 1:
            yield '##gff-version 3' + '\n'
        yield '\t'.join(content) + '\n'
        n += 1

def hmmBED2GFF(hmmBED):
    n = 1
    for rec in hmmBED:
        outstring = '\t'.join([rec[0], 'frisk_' + FRISK_VERSION, str(rec[3]), str(rec[1]), str(rec[2]), '.', '+', '.', 'ID=' + rec[3] + '_' + str(n).zfill(len(str(len(hmmBED))))]) + '\n'
        if n == 1:
            yield '##gff-version 3' + '\n'
        yield outstring
        n += 1

def thresholdList(intervalList, threshold, args, threshCol=3, merge=True):
    if args.findSelf:
        tItems = [t for t in intervalList if np.log10(t[3]) <= threshold]
    else:
        tItems = [t for t in intervalList if np.log10(t[3]) >= threshold]
    sortItems = sorted(tItems, key=itemgetter(0, 1, 2))
    if args.RIP:
        strItems = [ (w,x,y,str(z),str(c),str(p),str(s)) for w,x,y,z,c,p,s in sortItems ]
    else:
        strItems = [ (w,x,y,str(z)) for w,x,y,z in sortItems ]
    anomaliesBED  = pybedtools.BedTool(strItems)
    if merge:
        anomalies = anomaliesBED.merge(d=args.mergeDist, c='4,4,4', o='max,min,mean')
    else:
        anomalies = anomaliesBED
    return anomalies

def thresholdRIP(intervalList, args):
    logging.info('Extracting RIP features: CRImin = %s, PImin = %s, SImax = %s, CRIpeak = %s \
    ' % (str(args.minCRI), str(args.minPI), str(args.maxSI), str(args.peakCRI)))
    # Threshold windows on CRI, PI, and SI
    tCRImin	= [t for t in intervalList if t[6] >= args.minCRI] #Indexed from 0
    tPImin	= [t for t in tCRImin if t[4] >= args.minPI]
    tSImax	= [t for t in tPImin if t[5] <= args.maxSI]
    sBasic = sorted(tSImax, key=itemgetter(0, 1, 2))
    sBasic =[ (w,x,y,str(z),str(c),str(p),str(s)) for w,x,y,z,c,p,s in sBasic ]
    # Extract windows scoreing > peak CRI value
    tCRIpeak = [t for t in intervalList if t[6] >= args.peakCRI]
    sPeaks = sorted(tCRIpeak, key=itemgetter(0, 1, 2))
    sPeaks =[ (w,x,y,str(z),str(c),str(p),str(s)) for w,x,y,z,c,p,s in sPeaks ]
    # name, start, stop, KLI max, PI min, SI max, CRI min, CRI max
    if len(sBasic) > 0 and len(sPeaks) > 0:
        RIPbasic	= pybedtools.BedTool(sBasic).merge(d=0, c='4,5,6,7,7', o='max,min,max,min,max') #Indexed from 1
        RIPpeaks	= pybedtools.BedTool(sPeaks)
        #	list_idx	Bed_Idx	Field_Label
        #	0 			[1] 	name
        #	1 			[2] 	start
        #	2 			[3] 	stop
        #	3 			[4] 	windowKLI
        #	4 			[5] 	PI
        #	5 			[6] 	SI
        #	6 			[7] 	CRI
        # Return merged RIP annotations containing at least one peak CRI window.
        return RIPbasic.window(b=RIPpeaks, w=0, u=True)
    else:
        logging.info('No RIP features detected.')
        return None

def	wIndex(index, windows):
    try:
        indexCol = np.asarray([x[index] for x in windows])
    except IndexError:
        logging.error('Index field does not exist in window data.\n')
        sys.exit(1)
    return indexCol

def meanRangeLen(l):
    lens = list()
    for i in l:
        lens.append(i[1]-i[0])
    mean = reduce(lambda x, y: x + y, lens) / len(lens)
    return mean

def updateHMM(smallHMM, bigThresh): #BED interval objects. A = Fine scale guide, B = Starting annotations
    '''Update annotation boundaries in B using nearest interval boundary in A.
        Note: Not currently in use.'''
    hmmBounds = dict()
    for i in smallHMM:
        if i[0] not in hmmBounds:
            hmmBounds[i[0]] = list()
            hmmBounds[i[0]].append(int(i[1]))
            hmmBounds[i[0]].append(int(i[2]))
        else:
            hmmBounds[i[0]].append(int(i[1]))
            hmmBounds[i[0]].append(int(i[2]))
    updatedBoundaries = list()
    for y in bigThresh:
        newLeft  	= min(hmmBounds[y[0]], key=lambda x: abs(x-int(y[1])))
        newRight 	= min(hmmBounds[y[0]], key=lambda x: abs(x-int(y[2])))
        updatedBoundaries.append((y[0], str(newLeft), str(newRight)))
    newAnnotations 	= pybedtools.BedTool(updatedBoundaries)
    return newAnnotations

def	hmm2BED(allWindows, model, dataCol=3):
    allIntervals	=	list()
    # Get unique scaffold names from tuples in allWindows
    scaffoldList	= 	[i for i in Counter(elem[0] for elem in allWindows)]
    for name in scaffoldList:
        scaffoldWin = 	sorted([rec for rec in allWindows if rec[0] == name], key=itemgetter(0, 1, 2))
        dataSeq		= 	wIndex(dataCol, scaffoldWin)
        dataSeq		=	dataSeq[:, np.newaxis]
        stateSeq	= 	model.predict(dataSeq)
        stateRange1 = 	findBaseRanges(stateSeq, 0)
        stateRange2 = 	findBaseRanges(stateSeq, 1)
        for i in range2interval(stateRange1, scaffoldWin, 'State1'):
            allIntervals.append(i)
        for i in range2interval(stateRange2, scaffoldWin, 'State2'):
            allIntervals.append(i)
    allIntervals 	= 	sorted(allIntervals, key=itemgetter(0, 1, 2))
    allBED 			= 	pybedtools.BedTool(allIntervals)
    return allBED

def range2interval(rangeList, scaffoldWindows, state):
    for block in rangeList:
        yield (str(scaffoldWindows[0][0]), str(scaffoldWindows[block[0]][1]), str(scaffoldWindows[block[1]][2]), str(state))

def flattenKmerMap(kMap, window=1, seqLen=1, kmin=1, kmax=5, prop=False):
    '''Takes kmer count map (nested dict), flattens dictionaries to an array
        and scales counts to a standard window length.'''
    kCounts = list()
    for d in kMap:
        # Check kmer length for first key in current dictionary is in range.
        if len(d.keys()[0]) in range(kmin, kmax+1):
            if prop:
                # Option to convert counts in d to proportions
                kProps = {k: float(v)/sum(d.values()) for k, v in d.items()}
                for x in itertools.chain(kProps.values()):
                    kCounts.append(x)
            else:
                for x in itertools.chain(d.values()):
                    kCounts.append(x)
    if prop:
        return	np.array(kCounts)
    else:
        return	(float(window)/seqLen) * np.array(kCounts)

def pcaAxisLabels(pca_X, kmin=1, kmax=6):
    'Generate axis labels for PCA. Format: "PC1: AT (90.01%)"'
    blankmap 	= rangeMaps(kmin, kmax)
    keylist 	= []
    axisLabels 	= []
    idx 		= 0 #PC index
    for x in blankmap:
        for y in x.keys():
            keylist.append(y)
    for x in pca_X.explained_variance_ratio_:
        labelString = 'PC%s: %s (%s%%)' % (str(idx+1), str(keylist[np.argmax(pca_X.components_[idx])]), str(round(x*100,2)))
        axisLabels.append(labelString)
        idx += 1
    return axisLabels

def makeScatter(Y,args,pca_X=None,y_pred=None):
    if not args.culster:
        if args.runProjection == 'PCA':
            axisLabels = pcaAxisLabels(pca_X, kmin=args.pcaMin, kmax=args.pcaMax)
            # Generate graphic
            plt.figure()
            plt.title('PCA: kmer counts in anomalous regions')
            plt.xlabel(axisLabels[0])
            plt.ylabel(axisLabels[1])
            plt.scatter(Y[:, 0], Y[:, 1], color='blue', alpha=0.4, label='kmer_Anomaly')

        elif args.runProjection == 'TSNE':
            plt.figure()
            plt.title('t-SNE: kmer counts in anomalous regions')
            plt.xlabel('X')
            plt.ylabel('Y')
            plt.scatter(Y[:, 0], Y[:, 1], color='blue', alpha=0.4, label='kmer_Anomaly')

    elif args.culster:
        colors = np.array([x for x in 'bgrcmykbgrcmykbgrcmykbgrcmyk'])
        colors = np.hstack([colors] * 20)
        markers = np.array([x for x in '+o^8s^Do+o^8s^Do+o^8s^Do+o^8s^Do*'])
        markers = np.hstack([markers] * 20)

        # List dbscan classes
        y_class = np.unique(y_pred)
        # List to hold scatterplots for each class
        scatter_class = np.unique(y_pred).tolist()

        # List to hold generated class names
        names_class = np.unique(y_pred).tolist()
        count_class = 0
        for i in range(0, len(names_class)):
            count_class += 1
            if names_class[i] == -1:
                names_class[i] = 'Unclassified'
            else:
                names_class[i] = 'Class_' + str(i)

        fig = plt.figure() #figsize=(8,8)
        ax = fig.add_subplot(111)
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

        # Generate scatterplots for each class using pca coords
        for i in y_class:
            X_class = Y[np.where(y_pred == i)]
            scatter_idx = np.where(scatter_class == i)[0][0]
            scatter_class[scatter_idx] = plt.scatter(X_class[:, 0], X_class[:, 1], color=colors[i], marker=markers[i], s=20)

        ax.legend(scatter_class, names_class, loc='center left', bbox_to_anchor=(1, 0.5))

        if args.runProjection == 'PCA':
            axisLabels = pcaAxisLabels(pca_X, kmin=args.pcaMin, kmax=args.pcaMax)
            plt.title('PCA: kmer counts in anomalous regions')
            plt.xlabel(axisLabels[0])
            plt.ylabel(axisLabels[1])

        elif args.runProjection == 'TSNE':
            plt.title('t-SNE: kmer counts in anomalous regions')

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
        yrange = (y_positions[chrom], height)
        xranges = group[['start', 'width']].values
        yield BrokenBarHCollection(
            xranges, yrange, facecolors=df['colors'], **kwargs)
    if del_width:
        del df['width']

def makeChrPainting(selfGenome, args, anomBED, showGfffeatures=False):
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
    # Decide which chromosomes to use
    chromosome_list = natural_sort(selfGenome.keys())
    chromosome_list = [str(i) for i in chromosome_list]
    #Get chr lengths
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

    #Scaffold background
    scaffolds = pd.DataFrame(list(chromo_dict.iteritems()),columns=['chrom','end'])
    scaffolds['start'] = 0
    scaffolds['width'] = scaffolds.end - scaffolds.start
    scaffolds['colors'] = '#eeeeee'

    ##Read in anomalies
    features = pd.read_table(anomBED.fn, names=['chrom', 'start', 'end'],usecols=range(3))
    # Filter out chromosomes not in our list
    features = features[features.chrom.apply(lambda x: x in chromosome_list)]
    features['width'] = features.end - features.start
    features['colors'] = '#ff6600'

    if showGfffeatures and args.gffIn and args.gffFeatures:
        if type(args.gffFeatures) is list:
            gffType = args.gffFeatures[0]
        else:
            gffType = args.gffFeatures

        filteredGff = list()

        handle = open(args.gffIn)
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                break
            if not line.startswith("#"):
                rec = line.split()
                reclist = gffFilter(rec, feature=gffType)
                if not reclist:
                    continue
                recInterval = (str(reclist[0]),int(reclist[3]),int(reclist[4]))
                filteredGff.append(recInterval)

        handle.close()

        genes = pd.DataFrame.from_records(filteredGff, columns=['chrom', 'start', 'end'])
        genes = genes[genes.chrom.apply(lambda x: x in chromosome_list)]
        genes['width'] = genes.end - genes.start
        genes['colors'] = '#2243a8'

    fig = plt.figure()
    ax = fig.add_subplot(111)
    for collection in chromosome_collections(scaffolds, chrom_ybase, chrom_height):
        ax.add_collection(collection)

    for collection in chromosome_collections(features, chrom_ybase, chrom_height):
        ax.add_collection(collection)

    if showGfffeatures and args.gffIn and args.gffFeatures:
        for collection in chromosome_collections(
            genes, gene_ybase, gene_height, alpha=0.5, linewidths=0
        ):
            ax.add_collection(collection)

    # Axes tweaking
    ax.set_yticks([chrom_centers[i] for i in chromosome_list])
    ax.set_yticklabels(chromosome_list)
    ax.get_xaxis().set_major_formatter(FuncFormatter(lambda x, p: format(int(x), ',')))
    ax.axis('tight')

def mainArgs():
    """Process command-line arguments"""

    parser = argparse.ArgumentParser(
        description='Calculate all kmers in a given sequence',
        prog='frisk',
    )

    parser.add_argument('--version', 
                        action='version', 
                        version='frisk --' + str(FRISK_VERSION))
    # Inputs
    parser.add_argument('-H',
                        '--hostSeq',
                        type=str,
                        required=True,
                        help='The input host sequences (single species)')
    parser.add_argument('-Q',
                        '--querySeq',
                        type=str,
                        default=None,
                        help='Detect anomalous regions in this sequence by comparison \
                    	to hostSeq. Defaults to hostSeq.')
    parser.add_argument('--gffIn',
                        type=str,
                        default=None,
                        help='Path to GFF file with annotations for genome being frisked.')

    # Outputs
    parser.add_argument('-O',
                        '--outfile',
                        type=str,
                        help='Write KLI-IVOM bed track to this file')
    parser.add_argument('-t',
                        '--tempDir',
                        type=str,
                        default='temp',
                        help='Name of temporary directory')
    parser.add_argument('--gffOutfile',
                        type=str,
                        default=None,
                        help='Filename: Write merged anomaly annotations to gff3.')
    parser.add_argument('--hmmOutfile',
                        type=str,
                        default='2StateHmm.gff3',
                        help='Filename: Write hmm determined state features to gff3.')
    parser.add_argument('--graphics',
                        type=str,
                        default=None,
                        help='Name of file to print summary graphics to.')

    # Output options
    parser.add_argument('--mergeDist',
                        type=int,
                        default=0,
                        help='Merge anomalies annotations within x bases of each other.')
    parser.add_argument('--gffFeatures',
                        type=str,
                        default=None,
                        nargs='+',
                        help='Space delimited list of feature types to report from gff')
    parser.add_argument('--gffRange',
                        type=int,
                        default=0,
                        help='Report gff annotations within window of anomalous features')

    # Core settings
    parser.add_argument('-m',
                        '--minWordSize',
                        type=int,
                        default='1',
                        help='Minimum value of DNA word length')
    parser.add_argument('-k',
                        '--maxWordSize',
                        type=int,
                        default='8',
                        help='Maxmimum value of DNA word length')
    parser.add_argument('-w',
                        '--windowlen',
                        type=int,
                        default='5000',
                        help='Lenght of survey window')
    parser.add_argument('-i',
                        '--increment',
                        type=int,
                        default='2500',
                        help='Slide survey window by this increment')

    # Optional run settings
    parser.add_argument('--exitAfter',
                        default=None,
                        choices=[None, 'GenomeKmers', 'WindowKLI'],
                        help='Exit after completing task.')
    parser.add_argument('--recalc',
                        action='store_false',
                        default=True,
                        help='Force recalculation of reference sequence kmer counts if set. \
                    	Default uses counts from previous run.')
    parser.add_argument('--recalcWin',
                        action='store_false',
                        default=True,
                        help='Force recalculation of KLI score for specified window length and icrement if set. \
                    	Default uses window KLIs from previous run.')
    parser.add_argument('--scaffoldsAll',
                        action='store_true',
                        default=False,
                        help='If genomic scaffold is below minimum window size, process \
                    	whole scaffold as single window.')

    # Thresholding options for KLI
    parser.add_argument('--threshTypeKLI',
                        default= None,
                        choices=[None, 'percentile', 'otsu'],
                        help='Options for defining non-self threshold on window KLI scores: Otsu binarisation, 2-state HMM, percentile.')
    parser.add_argument('--percentileKLI',
                        type=float,
                        default=99.0,
                        help='Percentile at which to threshold window KLI scores. \
                    	By default, lower 99 percent windows are excluded.')
    parser.add_argument('--hmmKLI',
                        action='store_true',
                        default=False,
                        help='Run 2-state HMM to classify windows as self or non-self based on KLI score.')
    parser.add_argument('-F',
                        '--forceThresholdKLI',
                        type=float,
                        default=None,
                        help='KLI score above which a window will be considered anomalous. \
                    	Note: Given as raw KLI, not log10(KLI).')

    # RIP options
    parser.add_argument('--RIP',
                        action='store_true',
                        default=False,
                        help='Calculate and report RIP indicies for all windows + report GFF3 location \
                    	of RIP features, using same thresholding method options as for kmer anomalies.')
    parser.add_argument('--RIPgff',
                        type=str,
                        default='RIP_annotation.gff3',
                        help='Export RIP annotation track as gff3.')
    parser.add_argument('--minCRI',
                        type=float,
                        default=0.0,
                        help='Return windows where CRI score > than this value.'
                        )
    parser.add_argument('--peakCRI',
                        type=float,
                        default=1.0,
                        help='Require that merged annotations contain at least one window with a CRI score above this score.'
                        )
    parser.add_argument('--minPI',
                        type=float,
                        default=1.0,
                        help='Return windows where Product Index score > than this value.'
                        )
    parser.add_argument('--maxSI',
                        type=float,
                        default=1.0,
                        help='Return windows where Substrate Index score < than this value.'
                        )

    # Dimentionality Reduction Stuff
    parser.add_argument('--runProjection',
                        default=None,
                        choices=[None, 'PCA', 'TSNE' ],
                        help='Project anomalous windows into multidimensional space.')
    parser.add_argument('--dimReduce',
                        default='windows',
                        choices=['features', 'windows'],
                        help='Run dimensionality reduction on all anomalous features, or on merged features.')
    parser.add_argument('--culster',
                        default=None,
                        choices=[None, 'DBSCAN'], #add support for kmeans later
                        help='Attempt clustering on projection.')
    parser.add_argument('--dumpPCAdata',
                        action='store_true',
                        default=False,
                        help='If set export array of kmer proportional counts, and array of labels to temp directory.')
    parser.add_argument('--spikeNormal',
                        action='store_true',
                        default=False,
                        help='Include a sampling of windows from the centre of the self population.')
    parser.add_argument('--pcaMin',
                        type=int,
                        default='1',
                        help='Minimum order kmer set to include in dimensionality reduction data.')
    parser.add_argument('--pcaMax',
                        type=int,
                        default='6',
                        help='Maxmimum order kmer set to include in dimensionality reduction data.')
    parser.add_argument('--perplexity',
                        type=float,
                        default=20.0,
                        help='Used for TSNE. Typical values for the perplexity range between 5 and 50. Test larger values for larger/denser datasets.')
    parser.add_argument('--epsDBSCAN',
                        type=float,
                        default=10,
                        help='Set eps value for DBSCAN, will attempt to cluster points < eps apart. Note: Typically smaller for PCA (0.01 - 3.0), than for t-SNE (5-20).')

    # Revise feature boundaries using HMM split track with reduced window and increment size
    # Note: Not yet implemented
    parser.add_argument('--updateHMM',
                        action='store_true',
                        default=False,
                        help='Revise ')
    parser.add_argument('--updateWin',
                        type=int,
                        default=1000,
                        help='')
    parser.add_argument('--updateInc',
                        type=int,
                        default=500,
                        help='')

    # Reverse mode: Finds windows in query genome with KLI BELOW threshold. i.e. MORE like the host genome.
    parser.add_argument('--findSelf',
                        action='store_true',
                        default=False,
                        help='Finds windows in query genome with KLI BELOW threshold. i.e. MORE like the host genome.')

    args = parser.parse_args()
    if args.minWordSize > args.maxWordSize:
        logging.error('[ERROR] Minimum kmer size (-c/--minKmerSize) must be less than Maximum kmer size (-k/--maxKmerSize)\n')
        sys.exit(1)
    return args

def main():
    ##########################################
    ########## Initial housekeeping ##########
    ##########################################
    
    #Note: Move a lot of this off to class object
    logging.basicConfig(level=logging.INFO, format=("%(asctime)s - %(funcName)s - %(message)s"))
    args = mainArgs()
    print('frisk --', FRISK_VERSION)
    
    # Make path to store genome kmer calculations
    genomepickle 	= makePicklePath(args, space='genome')
    # Make path to store window KLI calculations
    windowsPickle 	= makePicklePath(args, space='window')

    # Set query sequence as self if none provided
    if not args.querySeq:
        querySeq = args.hostSeq
    else:
        querySeq = args.querySeq

    # Check temp file exists
    tempPathCheck(args)

    # Generate blank kmer dictionary
    blankMap = rangeMaps(args.minWordSize, args.maxWordSize)

    # Read in query genome sequences as dict keyed by seq name
    selfGenome, nBlocks = getFasta(querySeq)

    ##########################################
    #########  Import or Calculate   #########
    #########   Genome kmer Counts   #########
    ##########################################

    #Check if genome kmer dict previously generated
    if os.path.isfile(genomepickle) and args.recalc:
        logging.info('Importing previously calculated genome kmers from %s' % genomepickle)
        genomeKmers = pickle.load( open( genomepickle, "rb" ) )
    else:
        logging.info('Calculating kmers for host sequence: %s' % args.hostSeq)
        genomeKmers = computeKmers(args, genomepickle=genomepickle, window=None, genomeMode=True, kmerMap=blankMap, getMeta=True)
        if args.exitAfter == 'GenomeKmers':
            logging.info('Finished counting kmers. Exiting.')
            sys.exit(0)
        else:
            logging.info('Finished counting kmers.')

    ##########################################
    ###########   Extract Windows   ##########
    ############  Calc KLI Score  ############
    ##########################################

    if os.path.isfile(windowsPickle) and args.recalcWin:
        logging.info('Importing previously calculated window KLI scores from: %s' % windowsPickle)
        allWindows = pickle.load(open(windowsPickle, "rb"))

    else:
        # Initialise textfile output
        out     = args.outfile
        outPath = os.path.join(args.tempDir, out)
        handle  = open(outPath, "w")

        # List to store KLI-by-window
        allWindows = list()

        # Loop over genome to get all window KLI scores
        for seq, name, start, stop in crawlGenome(args, querySeq): #Note: coords are inclusive i.e. 1:10 = 0:9 in string
            target 			= [(name, seq)]
            windowKmers 	= computeKmers(args, genomepickle=None, window=target, genomeMode=False, kmerMap=blankMap, getMeta=True)
            GenomeIVOM  	= IvomBuild(windowKmers, args, genomeKmers, True)
            windowIVOM  	= IvomBuild(windowKmers, args, genomeKmers, False)
            windowKLI		= KLI(GenomeIVOM, windowIVOM, args)
            if args.RIP:
                PI, SI, CRI 	= calcRIP(windowKmers)
                outString	= [name, str(start), str(stop), str(windowKLI), str(PI), str(SI), str(CRI)]
                allWindows.append((name, start, stop, windowKLI, PI, SI, CRI))
            else:
                outString	= [name, str(start), str(stop), str(windowKLI)]
                allWindows.append((name, start, stop, windowKLI)) #Mind that KLI does not get stored as scientific notation

            handle.write('\t'.join(outString) + '\n')
            print('\t'.join(outString))

        # Close text output
        handle.close()

        # Write KLI-by-window annotations to file
        logging.info('Saving calculated window KLI scores as: %s' % windowsPickle)
        pickle.dump(allWindows, open(windowsPickle, "wb"))

        if args.exitAfter == 'WindowKLI':
            logging.info('Finished calculating window KLI scores. Exiting.')
            sys.exit(0)
        else:
            logging.info('Finished calculating window KLI scores.')

    ##########################################
    ###########  Calculate and Set  ##########
    #############  RIP threshold  ############
    ##########################################
    #Note: Index error if importing windows that do not include RIP data, add index check + sysExit to wIndex
    if args.RIP:
        PI		=	wIndex(4, allWindows)
        SI		=	wIndex(5, allWindows)
        CRI		=	wIndex(6, allWindows)
        sqrCRI 	= 	CRI**2
        logCRI 	= 	np.log10(CRI)
        # Scrub NANs
        PI		=	PI[~np.isnan(PI)]
        SI		=	SI[~np.isnan(SI)]
        CRI		=	CRI[~np.isnan(CRI)]
        sqrCRI	=	sqrCRI[~np.isnan(sqrCRI)]
        logCRI 	= 	logCRI[np.logical_not(np.isnan(logCRI))]

    ##########################################
    ###########  Calculate and Set  ##########
    #############  KLI threshold  ############
    ##########################################

    # Get KLI data for all surveyed windows
    allKLI 			=	wIndex(3, allWindows)
    logKLI 			=	np.log10(allKLI)
    KLIthreshold	=	None

    # Calc optimal bins for KLI data
    if FDBins(logKLI) < 30:
        optBins = 30
    else:
        optBins = FDBins(logKLI)

    if args.threshTypeKLI or args.forceThresholdKLI:
        # Check for user specified KLI threshold
        if args.forceThresholdKLI:
            KLIthreshold = np.log10(float(args.forceThresholdKLI))
            logging.info('Forcing log10(KLI) threshold = %s' % str(KLIthreshold))

        # Use Otsu method to set threshold
        elif args.threshTypeKLI == 'otsu':
            logging.info('Calculating optimal KLI threshold by Otsu binarization.')
            if FDBins(logKLI) < 10: #Calculate optimal number of bins for logKLI set using Freedman-Diaconnis method.
                logging.warning('[WARNING] Low variance in log10(KLI) data: Review data distribution, \
                    			consider percentile or manual thresholding.')
                KLIthreshold = otsu(logKLI, optBins) #Attempt Otsu, though probably not a good idea.
                logging.info('Optimal log10(KLI) threshold = %s' % str(KLIthreshold))
            else:
                KLIthreshold = otsu(logKLI, optBins) #Use Otsu binarization to calc optimal log10 threshold for weird KLI scores
                logging.info('Optimal log10(KLI) threshold = %s' % str(KLIthreshold))

        # Optional set threshold at percentile
        elif args.threshTypeKLI == 'percentile':
            KLIthreshold = np.percentile(logKLI, args.percentileKLI)
            logging.info('Setting threshold at %s percentile of log10(KLI)= %s' % (str(args.percentileKLI), str(KLIthreshold)))

    ##########################################
    ##########  Set anomaly feature  #########
    ####### boundaries using 2-state HMM  ####
    ##########################################

    if args.hmmKLI:
        # Create 2-state HMM model
        model	=	hmm.GaussianHMM(n_components=2, covariance_type="full")
        # Train on KLI strings from all scaffolds
        dataKLI	=	np.asarray(allKLI)
        dataKLI	=	dataKLI[:, np.newaxis]
        model.fit(dataKLI) #Note: Ideally pass all scaffold data sequences in separately
        # Use model to prdict HMM state for all scaffolds
        hmmBED	= hmm2BED(allWindows, model, dataCol=3)
        # Write hmm states as GFF3 outfile
        handle = open(os.path.join(args.tempDir, args.hmmOutfile), "w")
        for i in hmmBED2GFF(hmmBED):
            handle.write(i)
        handle.close()

    #
    # Dimensionality reduction ############
    # on Anomaly kmer counts #############
    #

    # Generate cordinate list for anomalous regions
    if args.runProjection:
        if args.dimReduce == 'features':
            anomWin  = thresholdList(allWindows, KLIthreshold, args, threshCol=3, merge=True)
        elif args.dimReduce == 'windows':
            # Anomalous windows by KLI threshold without merging
            anomWin  = thresholdList(allWindows, KLIthreshold, args, threshCol=3, merge=False)
        # Extract sequences
        logging.info('Recovering sequences %s for anomalous windows.' % str(len(anomWin)))
        anomSeqs = getBEDSeq(selfGenome, anomWin) #generator object that yields (name,seq) tuples.

        # Counter to manage empty array on first pass
        counter = 0

        for name, target in anomSeqs:
            logging.info('Computing kmers in %s.' % str(name))
            # Symetrical counting of kmers
            countMap 	= computeKmers(args, genomepickle=None, window=[(name, target)], genomeMode=False, kmerMap=blankMap, getMeta=False, sym=True)
            # Flatten kmer count dictionary to 1D array, and make counts proportional within each kmer length class
            sclCounts	= flattenKmerMap(countMap, window=args.windowlen, seqLen=len(target), kmin=args.pcaMin, kmax=args.pcaMax, prop=True)

            if counter >= 1:
                anomLabels = np.vstack([anomLabels, name])
                anomCounts = np.vstack([anomCounts, sclCounts])
            else:
                anomLabels = np.array(name)
                anomCounts = sclCounts
            counter += 1

        # Option to output arrays of counts and their labels for use in alternative dim reduce methods.
        if args.dumpPCAdata:
            pickle.dump(anomLabels, open(os.path.join(args.tempDir, 'anomLabels'), "wb"))
            pickle.dump(anomCounts, open(os.path.join(args.tempDir, 'anomCounts'), "wb"))

        # Run dimension reduction method of choice on anomaly kmer counts
        if args.runProjection == 'PCA':
            pca_X = PCA(n_components=3)
            Y = pca_X.fit(anomCounts).transform(anomCounts)
        elif args.runProjection == 'TSNE':
            Y = tsne.tsne(X=anomCounts, no_dims=3, initial_dims=50, perplexity=args.perplexity)
            pca_X = None

        if args.culster == 'DBSCAN':
            # Run clustering. Optional.
            dbscan = DBSCAN(eps=args.epsDBSCAN, min_samples=50).fit(Y)

            if hasattr(dbscan, 'labels_'):
                y_pred = dbscan.labels_.astype(np.int)
            else:
                y_pred = dbscan.predict(Y)
        else:
            y_pred = None

    ##########################################
    ###########  Threshold windows  ##########
    ######### and/or refine boundaries #######
    ##########################################

    # Threshold and merge anomalous features.
    if args.runProjection:
        anomalies = anomWin
        logging.info('Detected %s features above KLI threshold.' % str(len(anomalies)))
    elif KLIthreshold:
        anomalies = thresholdList(allWindows,KLIthreshold,args,threshCol=3,merge=True)
        logging.info('Detected %s features above KLI threshold.' % str(len(anomalies)))

    ##########################################
    ######## Write windows to GFF3 out #######
    ##########################################

    # Write anomalies as GFF3 outfile
    # Note: Need to incorporate cluster identity if calculated.
    if args.gffOutfile:
        handle  = open(os.path.join(args.tempDir, args.gffOutfile), "w")
        for i in anomaly2GFF(anomalies): #Use category='Class_Label' to pass class
            handle.write(i)
        handle.close()

    if args.RIP:
        RIPbed	= thresholdRIP(allWindows, args) #name, start, stop, KLI max, PI min, SI max, CRI min, CRI max
        if RIPbed:
            handle  = open(os.path.join(args.tempDir, args.RIPgff), "w")
            for i in RIP2GFF(RIPbed):
                handle.write(i)
            handle.close()

    # Find genes in anomalies if gff annotation file provided
    # Note: Separate output for each class, if runProjection and cluster are set
    if args.gffIn:

        if args.threshTypeKLI or args.forceThresholdKLI:
            # Change default name to be query centric, ok.
            gffOutname			= os.path.join(args.tempDir, "featuresIn_thresholded_Anomalies_" + os.path.basename(args.gffIn))
            features			= pybedtools.BedTool(args.gffIn).each(gffFilter, feature=args.gffFeatures)
            extractedFeatures 	= features.window(b=anomalies, w=args.gffRange, u=True)

            if len(extractedFeatures) > 0:
                extractedFeatures.saveas(gffOutname)
                logging.info('Successfully extracted %s features from within %sbp of anomaly annotations.' % (str(len(extractedFeatures)), str(args.gffRange)))

            else:
                logging.info('No features from %s detected within %s bases of anomalies.' % (args.gffIn, str(args.gffRange)))

        if args.hmmKLI:

            # State1 and State2 live here
            hmmOutputFile = os.path.join(args.tempDir, args.hmmOutfile)

            gffOutname1			= os.path.join(args.tempDir, "featuresIn_hmm_State1_" + os.path.basename(args.gffIn))
            bedState1			= pybedtools.BedTool(hmmOutputFile).each(gffFilter, feature='State1')
            features1			= pybedtools.BedTool(args.gffIn).each(gffFilter, feature=args.gffFeatures)
            extractedFeatures1 	= features1.window(b=bedState1, w=args.gffRange, u=True)

            gffOutname2			= os.path.join(args.tempDir, "featuresIn_hmm_State2_" + os.path.basename(args.gffIn))
            bedState2			= pybedtools.BedTool(hmmOutputFile).each(gffFilter, feature='State2')
            features2			= pybedtools.BedTool(args.gffIn).each(gffFilter, feature=args.gffFeatures)
            extractedFeatures2 	= features2.window(b=bedState2, w=args.gffRange, u=True)

            if len(extractedFeatures1) > 0:
                extractedFeatures1.saveas(gffOutname1)
                logging.info('Successfully extracted %s features from within %sbp of State1 hmm annotations.' % (str(len(extractedFeatures1)), str(args.gffRange)))

            else:
                logging.info('No features from %s detected within %s bases of State1 hmm features.' % (args.gffIn, str(args.gffRange)))

            if len(extractedFeatures2) > 0:
                extractedFeatures2.saveas(gffOutname2)
                logging.info('Successfully extracted %s features from within %sbp of State2 hmm annotations.' % (str(len(extractedFeatures2)), str(args.gffRange)))

            else:
                logging.info('No features from %s detected within %s bases of State2 hmm features.' % (args.gffIn, str(args.gffRange)))

    ##########################################
    ######### Write Graphics to File #########
    ##########################################
    ##########################################
    if args.graphics:
        logging.info('Writing graphics to %s' % args.graphics)

        with PdfPages(os.path.join(args.tempDir, args.graphics)) as pdf:
            plt.figure()
            plt.title('Raw KLI Distribution')
            sns.set(color_codes=True)
            sns.distplot(allKLI, hist=True, bins=optBins, kde=False, rug=False, color="b")
            pdf.savefig()  # saves the current figure into a pdf page
            plt.close()

            plt.figure()
            plt.title('log10(KLI) Distribution Optimized Bins')
            sns.set(color_codes=True)
            sns.distplot(logKLI, hist=True, bins=optBins, kde=False, rug=False, color="b")
            if KLIthreshold:
                plt.axvline(KLIthreshold, color='r', linestyle='dashed', linewidth=2)
            pdf.savefig()
            plt.close()

            plt.figure()
            plt.title('log10(KLI) Distribution Fine Bins')
            sns.set(color_codes=True)
            sns.distplot(logKLI, hist=True, bins=100, kde=False, rug=False, color="b")
            if KLIthreshold:
                plt.axvline(KLIthreshold, color='r', linestyle='dashed', linewidth=2)
            pdf.savefig()
            plt.close()

            #Make chromosome painting
            makeChrPainting(selfGenome, args, anomalies, showGfffeatures=True)
            plt.title('Chromosome painting: Anomalies on Query scaffolds')
            pdf.savefig(transparent=True)
            plt.close()
            
            #Compose Dim reduction scatterplot
            if args.runProjection:
                makeScatter(Y, args, pca_X=pca_X, y_pred=y_pred)
                pdf.savefig()
                plt.close()
                
            if args.RIP:
                plt.figure()
                plt.title('RIP Product Index')
                sns.set(color_codes=True)
                sns.distplot(PI, hist=True, bins=100, kde=False, rug=False, color="b")
                pdf.savefig()
                plt.close()

                plt.figure()
                plt.title('RIP Substrate Index')
                sns.set(color_codes=True)
                sns.distplot(SI, hist=True, bins=100, kde=False, rug=False, color="b")
                pdf.savefig()
                plt.close()

                plt.figure()
                plt.title('Composite RIP Index')
                sns.set(color_codes=True)
                sns.distplot(CRI, hist=True, bins=100, kde=False, rug=False, color="b")
                pdf.savefig()
                plt.close()

                plt.figure()
                plt.title('Square Composite RIP Index')
                sns.set(color_codes=True)
                sns.distplot(sqrCRI, hist=True, bins=100, kde=False, rug=False, color="b")
                pdf.savefig()
                plt.close()

                plt.figure()
                plt.title('log10 Composite RIP Index')
                sns.set(color_codes=True)
                sns.distplot(logCRI, hist=True, bins=100, kde=False, rug=False, color="b")
                pdf.savefig()
                plt.close()

            # Set the file's metadata via the PdfPages object:
            d = pdf.infodict()
            d['Title'] = 'Frisk: KLI Distribution Summary'
            d['Author'] = u'Frisk v0.0.1'
            d['Subject'] = 'Summary graphics'
            d['Keywords'] = 'histogram'
            d['CreationDate'] = datetime.datetime.today()
            d['ModDate'] = datetime.datetime.today()

    # Trace
    logging.info('Finished!')
