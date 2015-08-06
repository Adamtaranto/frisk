#!/usr/bin/env python
#frisk.py
#Version 1. Adam Taranto, July 2015
#Contact, Adam Taranto, adam.taranto@anu.edu.au

######################################################################################################################
# Detects regions of unusual sequence composition by comparison of local kmer frequencies to whole genome abundance. #
# For using in detecting genomic islands and potential segmental lateral gene transfer events.                       #
######################################################################################################################

##Completed functions
# Calculate all kmer frequencies (1-8mer) in query genome. [DONE]
# Calculate KLI 'selfiness score' for sliding window 5000bp, increment 2500. [DONE]
# Determine non-selifiness KLI threshold for genome. [DONE]

##Pending functions
# Collapse windows: merge overlapping non-self windows, export as gff3 (store KLI range from constituent windows)
# Cluster extracted non-self windows: 
	# Store kmer table per window, 
	# Normalise kmer counts by window size, 
	# D2/PCA clustering + output graphics.
	# K-means clustering??
#Classify
	# Determine kmers that are driving group separation (Positive control is RIP class, CpA --> TpA conversion = TA dinucleotide enrichment against background + CA depletion)
# Identifiy bayesian-changepoint boundaries to non-self clusters using Rpy2 link to changepoint (Alternative to collapsing windows by threshold).
# Mask N-blocks from annotation.
# Find genes (from gff) that are in or overlap each non-self group.

import argparse
import array
from collections import Counter
import copy
import gzip
import logging
import math
import numpy as np
from operator import itemgetter
import os
import os.path
import pickle
import pybedtools
import matplotlib.pyplot as plt
import rpy2
#import subprocess
import sys
#import shutil
#import threading
#import traceback

#######################
#######################
### Global contants ###
#######################
#######################

FRISK_VERSION = '0.0.1'

#SYMB_TESTING  = SYMB_NAME + '_testing.fasta'
#BINARIES_DIR  = 'binaries'

LETTERS = ('A', 'T', 'G', 'C')

def tempPathCheck(args):
    """Check if the temporary folder exists, else creates it"""

    logging.info('Checking for temporary folder')
    tempFolder = os.path.abspath(args.tempDir)
    if not os.path.isdir(tempFolder):
        os.makedirs(tempFolder)

def countN(sequence):
	#Count number of ATGC and non-ATGC bases in sequence string
	count  = 0
	countN = 0
	baseTally = Counter()

	baseTally += Counter(sequence)

	for i in baseTally:
		if i not in LETTERS:
			countN += baseTally[i]
		else:
			count += baseTally[i]
	#return tuple with readable base count and non-ATGC base count
	return (count,countN)

def iterFasta(path):
	#path = location of fasta to be processed
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
			#If this is the >= 2nd title then first seq has been read and it is time to yield 
			if name:
				yield (name, ''.join(seq))
			#Update the name to current seq and blank the list
			name = line[1:]
			seq = []
		else:
			#Extend the current sequence with uppercase sequence
			#Perhaps it should be an option to ignore masked sequences (may want to mask isochores or transposons?)
			seq.append(line.upper())
	if name:
		#Yield the final sequence at end of file.
		yield (name, ''.join(seq))
	handle.close()


def crawlGenome(args, querySeq):
	#Take genome fasta extract scaffolds
	#Extract incremented windows
	#yield (seq, scaffold name, window start, window end)

	w = args.windowlen
	i = args.increment
	saveSmalls = args.scaffoldsAll
	nSeqs = 0

	#Iterate over scaffolds
	for name,seq in iterFasta(querySeq):
		#Reset stats and counters
		windowCount = 0
		winExclude = 0
		bases,nonbase = countN(seq)
		size = len(seq)
		jumpback = False
		#Process scaffolds that below window size + minimum acceptable first increment, as a single window.
		if size <= w + ((w * 0.75) - i) and saveSmalls:
			scaffB,scaffN = countN(seq)
			if scaffN >= 0.3 * size:
				logging.info('%s excluded as > 30 percent unresolved sequence.' % name)
				winExclude += 1
				continue
			else:
				logging.info('Rescuing small scaffold %s' % name)
				yield (seq, name, 0, size)
				windowCount += 1
				nSeqs += 1
		else:
			#Crawl sequence in multiples of increment length until one increment from end of seq
			for j in xrange(0, size - i + 1, i):
				#If window overshoots sequence, jump back one w len from last base.
				if j+w > size:
					winSeq = seq[size - w: size]
					jumpback = True
				else:
					#Extract window of len w starting at current base
					winSeq = seq[j: j+w]
				#Screen N content = (baseCount, Ncount)
				winB,winN = countN(winSeq)
				if winN >= 0.3 * len(winSeq):
					logging.info('%s excluded as > 30 percent unresolved sequence.' % name)
					winExclude += 1
					continue
				elif jumpback:
					yield (winSeq, name, size - w, size)
				else:
					yield (winSeq, name, j+1, j + w)
				windowCount += 1
		#Iterate to next sequence in input fasta
		logging.info('Extracted %s windows from %s bases in %s' % (windowCount,str(size),name))
		logging.info('Excluded %s windows from %s' % (winExclude,name))
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

def rangeMaps(args):
	#Calls prepareMaps to write list of kmer dictionaries for range kmin to kmax
	kMin    = args.minWordSize
	kMax    = args.maxWordSize
	maps    = []
	logging.info('Preparing kmer maps')
	#Prepare possible kmer combinations for each len k
	for i in xrange(kMin, kMax + 1):
		maps.append(prepareMaps(0, i, ['']))

	return maps

def revComplement(kmer):
	revcompl = lambda x: ''.join([{'A':'T','C':'G','G':'C','T':'A'}[B] for B in x][::-1])
	return revcompl(kmer)

def computeKmers(args, path, genomepickle, window, genomeMode, kmerMap):
	#args, path to genome, tuple(windowID,Seq), fastabase.kmers, genomeMode = true/false, precalculated kmerMap)
	#windowKmers = computeKmers(args, None, None, target, False, blankMap)
	#genomeKmers = computeKmers(args, path, genomepickle, None , True, blankMap)
	"""Compute the kmer counts throughout the kmer range for the complete genome, and
	write the output to a file."""

	if sum(kmerMap[0].itervalues()) != 0:
		logging.info('kmer template is not blank!')
		sys.exit(1)

	if genomeMode:
		logging.info('Computing kmers for %s' % path)
		targetSeq = iterFasta(path)
	else:
		targetSeq = window

	excludedMaxMer = dict()
	excludedMaxMer['exMax'] = 0
	surveySeqLen = dict()
	surveySeqLen['totalLen'] = 0

	# Prepare all maps
	kMin    = args.minWordSize
	kMax    = args.maxWordSize
	#Turns out you need to redefine the dict being copied else alias edits orignal
	maps    = copy.deepcopy(kmerMap)

	# Iterate over sequences
	nSeqs   = 0
	#Read (name, sequence) tuple
	for name,seq in targetSeq:
		size   = len(seq)
		surveySeqLen['totalLen'] += size
		# For each kmer length
		for i in xrange(kMin, kMax + 1):
			# Crawl sequence one base at a time (stop k bases from end)
			for j in xrange(size - i + 1):
				#Extract word of len i (kmer) starting at current base
				word = seq[j:j + i]
				#This checks that the string is a legit kmer in the map i.e. does not contain Ns
				#.get will return 'None' if kmer key does not exist, or value if it does.
				#Updates idx original 0 count of word
				idx  = maps[i - kMin].get(word, None)
				# If the word contains characters other than ATGC, skip to next
				if idx is None:
					if i == kMax:
						excludedMaxMer['exMax'] += 1
					continue
				#Increment instances of current word
				maps[i - kMin][word] += 1
				#Increment count for reverse complement of current word
				if genomeMode:
					maps[i - kMin][revComplement(word)] += 1

		#Iterate to next sequence in input fasta
		nSeqs += 1

	#Add excluded 
	maps.append(surveySeqLen)
	maps.append(excludedMaxMer)

	if genomeMode:
		#Save genome kmers to pickle
		pickle.dump(maps, open(genomepickle, "wb"))
		# Trace
		logging.info('Processed %d sequences' % nSeqs)
	#return a dictionary of kmer counts
	return maps

def IvomBuild(windowKmers, args, GenomeKmers, isGenomeIVOM):
	# Genomelen = genome length excluding unresolveable bases (or should this be total maxmer space - MaxMers with Ns?)
	# Input is dictionary of all kmers in current window. Should include a final dict with N count for each kmer length.
	# Calculates weights Wi(=Counts*deg_freedom) and obs_freqs (Pi)

	#print "windowKmers: "
	#print windowKmers[0]
	#print "GenomeKmers: "
	#print GenomeKmers[0]

	windowlen   = args.windowlen
	klen        = args.maxWordSize

	genomeSpace = GenomeKmers[klen]['totalLen'] - (GenomeKmers[klen+1]['exMax'] * klen)
	windowSpace = windowKmers[klen]['totalLen'] - (windowKmers[klen+1]['exMax'] * klen)

	sumWindowIVOM = 0
	subKmers = dict()
	storeMaxMerIVOM = dict()

	for k in windowKmers[klen-1]:
		#Skip kmers that do not occur in window
		if windowKmers[klen-1][k] == 0:
			continue
		#print k
		count = windowKmers[klen-1][k]
		#print count
		subKmers['w'] = dict()
		subKmers['p'] = dict()

		if not isGenomeIVOM: #Calculating IVOM weights for window kmers against window counts
			for x in range(1,klen+1):
				#Process maxmer
				if x == klen:
					subKmers['w'][x] = count * 4^x
					subKmers['p'][x] = float(count) / ((windowSpace-(x-1)) * 2)
				elif x >= 2:
					subK = k[0:x] #Be sure to grab the first x bases in maxmer
					subKmers['w'][x] = windowKmers[x-1][subK] * 4^x
					subKmers['p'][x] = float(windowKmers[x-1][subK]) / ((windowSpace - (x-1)) * 2)
				else:
					subK = k[0] #Be sure to grab the first x bases in maxmer
					subKmers['w'][x] = windowKmers[x-1][subK] * 4^x
					subKmers['p'][x] = float(windowKmers[x-1][subK]) / ((windowSpace) * 2)

		else: #Calculating IVOM weights for current window kmer against genome counts
			for x in range(1,klen+1):
				#Process maxmer
				if x == klen:
					subKmers['w'][x] = GenomeKmers[x-1][k] * 4^x
					subKmers['p'][x] = float(GenomeKmers[x-1][k]) / ((genomeSpace - (x-1)) * 2)
				elif x >= 2:
					subK = k[0:x] #Be sure to grab the first x bases in maxmer
					subKmers['w'][x] = GenomeKmers[x-1][subK] * 4^x
					subKmers['p'][x] = float(GenomeKmers[x-1][subK]) / ((genomeSpace - (x-1)) * 2)
				else:
					subK = k[0] #Be sure to grab the first x bases in maxmer
					subKmers['w'][x] = GenomeKmers[x-1][subK] * 4^x
					subKmers['p'][x] = float(GenomeKmers[x-1][subK]) / (genomeSpace * 2)

		w_totals = 0
		for w in subKmers['w']:
			w_totals += subKmers['w'][w]

		scaledW = dict()
		for x in range(1,klen+1):
			scaledW[x]= float(subKmers['w'][x]) / w_totals

		kmerIVOM = dict()
		for x in range(1,klen+1):
			if x == 1:
				kmerIVOM[x] = scaledW[x] * subKmers['p'][x]
			elif x < klen:
				kmerIVOM[x] = scaledW[x] * subKmers['p'][x] + ((1-scaledW[x]) * kmerIVOM[x-1])
			else:
				kmerIVOM[x] = scaledW[x] * subKmers['p'][x] + ((1-scaledW[x]) * kmerIVOM[x-1])

		storeMaxMerIVOM[k] = kmerIVOM[klen]
		#Note: In Genome mode 'window' is equal to whole genome length (minus N-containing maxmers).
		sumWindowIVOM += kmerIVOM[klen]

	#Rescale each kmer from 0 to 1 (cause relative entropy)
	for k in storeMaxMerIVOM:
		storeMaxMerIVOM[k] = float(storeMaxMerIVOM[k]) / sumWindowIVOM

	#Return dictionary of scaled IVOM scores keyed by kmer
	return storeMaxMerIVOM

#Kullback-Leiber
def KLI(GenomeIVOM, windowIVOM,args):
#calculates relative entropy (Kullback-Leibler)
#IVOM structure: dict[kmer][IVOM Score]
	windowKLI = 0
	for k in windowIVOM:
		w = float(windowIVOM[k])
		G = float(GenomeIVOM[k])
		#print("Window IVOM: %s" % str(w))
		#print("Genome IVOM: %s" % str(G))
		#Absolute value removes direction of divergence
		#Negative number indicates kmer depleted in window relative to genome
		#Positive number indicates enriched in window relative to genome
		#Magnitude of KLI reflects degree of difference
		#Note: Using log2 instead of log10 to accentuate variance within small range.
		if args.absoluteKLI:
			windowKLI += (w*abs(math.log((w/G),2)))
		else:
			windowKLI += (w*math.log((w/G),2))
	return windowKLI

def calcRIP(windowKmers):
	#Product Index
	if windowKmers[1]['AT'] > 0:
		PI = windowKmers[1]['TA'] / float(windowKmers[1]['AT'])
	else:
		PI = None
	#Substrate index
	AC_GT = (windowKmers[1]['AC'] + windowKmers[1]['GT'])
	if AC_GT > 0:
		SI = (windowKmers[1]['CA'] + windowKmers[1]['TG']) / float(AC_GT)
	else:
		SI = None
	#Composite RIP Index
	if PI and SI:
		CRI = PI - SI
	else:
		CRI = None
	return (PI,SI,CRI)

def makePicklePath(args,space):
	#Note: need to add kmer range to filename
	pathbase = os.path.basename(args.hostSeq)
	pickleOut = os.path.join(args.tempDir,pathbase + "_" + str(args.minWordSize) + "_" + str(args.maxWordSize) + "_" + space + '.p')
	return pickleOut

def FDBins(data):
	#Freedman-Diaconis method for calculating optimal number of bins for dataset
	IQR = np.subtract(*np.percentile(data, [75, 25])) #Magic '*' unpacks tuple and passes two values to subtract function.
	bins = 2 * IQR * math.pow(len(data),(1.0/3.0)) #Number of bins
	bins = int(round(bins)) #Otsu needs rounded integer
	return bins

def otsu(data,fd):
	#data is array of log10(KLI) 
	raw = data
	data = np.atleast_1d(data)
	data = data[~ np.isnan(data)]
	data = data/(max(abs(data))*-1.0) #Scale to 0-1 and make positive
	hist,binEdges = np.histogram(data,bins=fd)
	hist = hist * 1.0 #Covert to floats
	hist_norm = hist.ravel()/hist.max() #Normalise hist to largest bin
	Q = hist_norm.cumsum()
	bins = np.arange(fd)
	fn_min = np.inf
	thresh = -1
	for i in xrange(1,fd): #Start from second position (1) to compare all to bin 0
		p1,p2 = np.hsplit(hist_norm,[i]) # Split normalised hist values into two brackets at bin i (bin 1 < i, bin 2 >=i)
		q1,q2 = Q[i-1], Q[fd-1] - Q[i-1] # cum sum of bin values #!! yields zero on final bin
		b1,b2 = np.hsplit(bins,[i]) # Split bins into 2 brackets at bin i
		# Finding means and variances
		m1,m2 = q1/len(p1),q2/len(p2)
		v1,v2 = np.sum(np.square(p1-m1))/len(p1),np.sum(np.square(p2-m2))/len(p2)
		#Calculates the minimization function
		fn = (v1*q1) + (v2*q2)
		if fn < fn_min:
			fn_min = fn
			thresh = i
	logging.info('OTSU selected bin %s as threshold position' % str(thresh))
	#Convert bact to log10(KLI)
	otsuNum = binEdges[thresh] * max(abs(raw)) * -1.0
	return otsuNum

def otsu_CP(data, bins):
	"""Compute a threshold using Otsu's method
	data           - an array of intensity values between zero and one
	bins           - we bin the data into this many equally-spaced bins, then pick
					 the bin index that optimizes the metric
	"""
	data = np.atleast_1d(data)
	data = data[~ np.isnan(data)]
	data.sort()

	var = running_variance(data)
	rvar = np.flipud(running_variance(np.flipud(data)))
	thresholds = data[1:len(data):len(data)/bins]
	score_low = (var[0:len(data)-1:len(data)/bins] * np.arange(0,len(data)-1,len(data)/bins))
	score_high = (rvar[1:len(data):len(data)/bins] * (len(data) - np.arange(1,len(data),len(data)/bins)))
	scores = score_low + score_high

	if len(scores) == 0:
		return thresholds[0]
	index = np.argwhere(scores == scores.min()).flatten()
	if len(index)==0:
		return thresholds[0]
	index = index[0]

	# Take the average of the thresholds to either side of
	# the chosen value to get an intermediate in cases where there is
	# a steep step between the background and foreground
	
	if index == 0:
		index_low = 0
	else:
		index_low = index-1
	if index == len(thresholds)-1:
		index_high = len(thresholds)-1
	else:
		index_high = index+1 
	#Note: Need to restore scaled data to negative log data 
	return (thresholds[index_low]+thresholds[index_high]) / 2

def running_variance(x):
	n = len(x)
	m = x.cumsum() / np.arange(1,n+1)
	x_minus_mprev = x[1:]-m[:-1]
	x_minus_m = x[1:]-m[1:]
	s = (x_minus_mprev*x_minus_m).cumsum()
	var = s / np.arange(2,n+1)
	return np.hstack(([0],var))

def gffFilter(rec, **kwargs):
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
		content = [str(i[0]),'frisk_' + FRISK_VERSION, annotType, str(i[1]),str(i[2]), '.', '+', '.', ';'.join(['ID=Anomaly_'+ str(n).zfill(8), 'maxKLI='+ str(i[3]),'minKLI='+ str(i[4]),'meanKLI='+ str(i[5])])]
		if n == 1:
			yield '##gff-version 3' + '\n'
		yield '\t'.join(content) + '\n'
		n += 1

def thresholdList(intervalList,threshold):
	tItems = [t for t in intervalList if np.log10(t[3]) >= threshold]
	sItems = sorted(tItems, key=itemgetter(0,1,2))
	anomaliesBED = pybedtools.BedTool(sItems)
	anomalies = anomaliesBED.merge(d=0, c='4,4,4', o='max,min,mean')
	return anomalies

def mainArgs():
	"""Process command-line arguments"""

	parser = argparse.ArgumentParser(description='Calculate all kmers in a given sequence')

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
	parser.add_argument('-O',
						'--outfile',
						type=str,
						help='Write KLI-IVOM bed track to this file')
	parser.add_argument('-G',
						'--gffOutfile',
						type=str,
						help='Write merged anomaly annotations to gff3.')
	parser.add_argument('-A',
						'--absoluteKLI',
						action='store_true',
						default=False,
						help='Calculate absolute value of KLI for each window. \
						Default False, will report net KLI movement.')
	parser.add_argument('-R',
						'--recalc',
						action='store_false',
						default=True,
						help='Force recalculation of reference sequence kmer counts if set. \
						Default uses counts from previous run.')
	parser.add_argument('-s',
						'--scaffoldsAll',
						action='store_true',
						default=False,
						help='If genomic scaffold is below minimum window size, process \
						whole scaffold as single window.')
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
	parser.add_argument('-t',
						'--tempDir',
						type=str,
						default='temp',
						help='Name of temporary directory')
	parser.add_argument('-g',
						'--gffPath',
						type=str,
						default=None,
						help='Path to GFF file with annotations for genome being frisked.')
	parser.add_argument('-r',
						'--gffRange',
						type=int,
						default=0,
						help='Report gff annotations within window of anomalous k-chores')
	parser.add_argument('-f',
						'--gffFeatures',
						type=str,
						default=None,
						nargs='+',
						help='Space delimited list of feature types to report from gff')
	parser.add_argument('-p',
						'--runProjection',
						default=None,
						choices=[None, 'PCA2', 'PCA3', 'MDS2', 'MDS3', 'D2' ],
						help='Project anomalous windows into multidimensional space.')
	parser.add_argument('-c',
						'--culster',
						default=None,
						choices=[None, 'DBSCAN', 'KMEANS', ],
						help='Attempt clustering of windows.')
	parser.add_argument('-n',
						'--spikeNormal',
						action='store_true',
						default=False,
						help='Include a sampling of windows from the centre of the self population.')
	parser.add_argument('-F',
						'--forceThresholdKLI',
						default=None,
						help='KLI score above which a window will be considered anomalous. \
						Note: Give as raw KLI, not log10(KLI).')
	parser.add_argument('-E',
						'--exitAfter',
						default=None,
						choices=[None, 'GenomeKmers'],
						help='Exit after completing task.')

	args = parser.parse_args()
	if args.minWordSize > args.maxWordSize:
		logging.error('[ERROR] Minimum kmer size (-c/--minKmerSize) must be less than Maximum kmer size (-k/--maxKmerSize)\n')
		sys.exit(1)
	return args

def main():
	#Test params
	#freakSeq.py --hostSeq zymoChr_1_2.fasta --outfile test_zymo.bed --minWordSize 1 --maxWordSize 8 --tempDir temp
	#number of seq defaulted to 0 so that whole seq is used for kmer calculation
	logging.basicConfig(level=logging.INFO, format=("%(asctime)s - %(funcName)s - %(message)s"))
	args = mainArgs()
	path = args.hostSeq
	genomepickle = makePicklePath(args, 'genome')

	if not args.querySeq:
		querySeq = args.hostSeq
	else:
		querySeq = args.querySeq

	tempPathCheck(args)
	blankMap = rangeMaps(args)

	#Check if genome kmer dict previously generated
	if os.path.isfile(genomepickle) and args.recalc:
		logging.info('Importing previously calculated genome kmers from %s' % genomepickle)
		genomeKmers = pickle.load( open( genomepickle, "rb" ) )
	else:
		logging.info('Calculating kmers for host sequence: %s' % args.hostSeq)
		genomeKmers = computeKmers(args, path, genomepickle, None , True, blankMap)
		if args.exitAfter == 'GenomeKmers':
			logging.info('Finished counting kmers. Exiting.')
			sys.exit(0)
		else:
			logging.info('Finished counting kmers.')

	#Initialise output
	out     = args.outfile
	outPath = os.path.join(args.tempDir,out)
	handle  = open(outPath, "w")

	allWindows = list()

	#First pass to get all window KLI scores
	for seq,name,start,stop in crawlGenome(args, querySeq): #note coords are inclusive i.e. 1:10 = 0:9 in string
		target = [(name,seq)]
		windowKmers = computeKmers(args, None, None, target, False, blankMap)
		GenomeIVOM = IvomBuild(windowKmers, args, genomeKmers, True)
		windowIVOM = IvomBuild(windowKmers, args, genomeKmers, False)
		windowKLI = KLI(GenomeIVOM, windowIVOM,args)
		outString = [name, str(start), str(stop), str(windowKLI)]
		handle.write('\t'.join(outString) + '\n')
		print('\t'.join(outString))
		allWindows.append((name,start,stop,windowKLI))#Mind that KLI does not get stored as scientific notation

	handle.close()

	windowsPickle = makePicklePath(args, 'windows')
	pickle.dump(allWindows, open(windowsPickle, "wb"))

	#Find optimal KLI threshold
	logging.info('Finding optimal threshold.')
	allKLI = [x[3] for x in allWindows] #Get KLI list
	logKLI = np.log10(allKLI) #log10 transform
	fd = FDBins(logKLI) #Calculate optimal number of bins for logKLI set using Freedman-Diaconnis method.

	#Set KLI threshold
	if args.forceThresholdKLI:
		KLIthreshold = np.log10(float(args.forceThresholdKLI))
		logging.info('Forcing log10(KLI) threshold = %s' % str(KLIthreshold))
	else:
		#Use Otsu binarization to calc optimal log10 threshold for weird KLI scores
		KLIthreshold = otsu(logKLI,fd)
		logging.info('Optimal log10(KLI) threshold = %s' % str(KLIthreshold))

	#Threshold and merge anomalous features.
	anomalies = thresholdList(allWindows,KLIthreshold)

	handle  = open(os.path.join(args.tempDir,args.gffOutfile), "w")
	for i in anomaly2GFF(anomalies):
		handle.write(i)
	handle.close()
	
	#If a gff3 annotation is provided 
	#make function def findOverlaps(args, anomBED):
	if args.gffPath:
		#Change default name to be query centric, ok.
		gffOutname = os.path.join(args.tempDir, "featuresInAnomalies_" + os.path.basename(args.gffPath))
		features = pybedtools.BedTool(args.gffPath).each(gffFilter, feature=args.gffFeatures)
		extractedFeatures = features.window(b=anomalies, w=args.gffRange, u=True)
		if len(extractedFeatures) > 0:
			extractedFeatures.saveas(gffOutname)
			logging.info('Successfully extracted %s features from within %s of anomalies' % (str(len(extractedFeatures)), str(args.gffRange)))
		else:
			logging.info('No features from %s detected within %s bases of anomalies' % (args.gffPath, str(args.gffRange)))

	##Next:
	'''	1) 	Save image to file [log10(KLI),PCA,PCA with clusters]. 

		2) 	Calculate and mask N-blocks from final merged anoms. (find continuous strings)
		
		3) 	Def getSeq(coords): #Feed in list of tuples / bedTool object
				for chr,start,stop in coords:

		4)	For given sequence, get kmer counts 
			and normalise to window len. (rowname= ) Pickle object.

		5)	Driving Kmer report. 
			For subset of windows rank kmers by mean IKW-KLI, print as report.

		6)	Run projection (PCA, MDS). Save image. (figure out how to save images!)
	
		7)	Run clustering on Projection. (DBSCAN, other methods)

		8)	Function to extract super normal windows to spike PCA.

		9)	Report coords for anomalies with category lable in type.

		10) Option to recycle KLI track previously calculated per scaffold >> continue to gff reports and PCA

		11)	Build RIP into output
	'''

	#Options to add:
		#...

	# Trace
	logging.info('Finished!')
	

if __name__ == '__main__':
	main()