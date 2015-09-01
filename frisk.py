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
#Classify
	# Determine kmers that are driving group separation (Positive control is RIP class, CpA --> TpA conversion = TA dinucleotide enrichment against background + CA depletion)
# Identifiy bayesian-changepoint boundaries to non-self clusters
# Mask N-blocks from annotation.
# Find genes (from gff) that are in or overlap each non-self group.

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
import numpy as np
from operator import itemgetter
import os
import os.path
#import pandas as pd
import pickle
import pybedtools
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy import stats
import seaborn as sns
import sys
import versioneer

#######################
#######################
### Global contants ###
#######################
#######################


FRISK_VERSION = '0.0.2'

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

def findBaseRanges(s, ch, name=None, minlen=0):
	#For string return intervals of character longer than minimum length.
	data = [i for i, ltr in enumerate(s) if ltr == ch]
	ranges = list()
	for k, g in itertools.groupby(enumerate(data), lambda (i,x):i-x):
		group =  map(itemgetter(1), g)
		if (group[-1] - group[0]) < minlen:
			continue
		else:
			if name:#Note: Might need to +1 to positions for bedtools coords that are not zero indexed
				ranges.append((name, group[0], group[-1]))
			else:
				ranges.append((group[0], group[-1]))
	return ranges #Format = [('ScaffName,start,stop'),('ScaffName,start,stop')]

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
			name = line.strip('>').split()[0]
			seq = []
		else:
			#Extend the current sequence with uppercase sequence
			#Perhaps it should be an option to ignore masked sequences (may want to mask isochores or transposons?)
			seq.append(line.upper())
	if name:
		#Yield the final sequence at end of file.
		yield (name, ''.join(seq))
	handle.close()

def getFasta(fastaPath):
	"""Write fasta to dictionary, key by scaffold name."""
	seqDict = dict()
	nRanges = list()
	for name,seq in iterFasta(fastaPath):
		seqDict[name] = seq
		nRanges.append(findBaseRanges(seq, 'N', name=name, minlen=10))
	nBlocks = pybedtools.BedTool(nRanges)
	return seqDict,nBlocks

def getBEDSeq(fastaDict,BEDintervals):
	"""Given BED object with scaffold coordinates, 
		fetch sequence from dictionary of sequences, keyed by scaffold name. """
	for rec in BEDintervals:
		if rec[0] in fastaDict:
			seq = fastaDict[rec[0]][int(rec[1])-1:int(rec[2])-1]
			name = ":".join([rec[0],str(rec[1]),str(rec[2])])
			if len(seq) > 0:
				yield (name,seq)
			else: 
				print('Retrieved zero len sequence for %s' % name)
		else:
			continue

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
		#Process scaffolds that fall below window size + minimum acceptable first increment, as a single window.
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
					logging.info('Window from %s excluded as > 30 percent unresolved sequence.' % name)
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
	# Genomelen = genome length excluding unresolveable max-kmer space (i.e. Total bases in max len kmers that contain Ns)
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
	#For each max len kmer
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
					#Note:Weighting could probably be achieved by observations multiplied by kmerlen squared i.e. count * x**2
					subKmers['w'][x] = count * 4**x # kmer_count * (DNA bases^seqlength)
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
			for x in range(1,klen+1):
				#Process maxmer
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
			#Running total of sub-kmer raw weights
			w_total += subKmers['w'][x]
			#Add total at eact position to dict
			w_running_totals[x] = w_total

		scaledW = dict()
		for x in range(1,klen+1):
			#Note: Weighting of sub-kmers not explictily addressed in Vernikos and Parkhill.
			scaledW[x]= float(subKmers['w'][x]) / w_running_totals[x]

		kmerIVOM = dict()
		for x in range(1,klen+1):
			if x == 1:
				kmerIVOM[x] = scaledW[x] * subKmers['p'][x]
			elif x < klen:
				kmerIVOM[x] = scaledW[x] * subKmers['p'][x] + ((1-scaledW[x]) * kmerIVOM[x-1])
			else: #IVOM for max len kmer
				kmerIVOM[x] = scaledW[x] * subKmers['p'][x] + ((1-scaledW[x]) * kmerIVOM[x-1])

		storeMaxMerIVOM[k] = kmerIVOM[klen]
		#Note: In Genome mode 'window' is equal to whole genome length (minus length of 'N'-containing maxmers).
		sumWindowIVOM += kmerIVOM[klen]

	#Rescale each kmer from 0 to 1 (cause relative entropy)
	for k in storeMaxMerIVOM:
		storeMaxMerIVOM[k] = float(storeMaxMerIVOM[k]) / sumWindowIVOM

	#Return dictionary of scaled IVOM scores keyed by kmer
	return storeMaxMerIVOM

def KLI(GenomeIVOM, windowIVOM, args):
#Kullback-Leiber Index: Measure of relative entropy.
#IVOM structure: dict[kmer][IVOM Score]
	windowKLI = 0
	for k in windowIVOM:
		w = float(windowIVOM[k])
		G = float(GenomeIVOM[k])
		#print("Window IVOM: %s" % str(w))
		#print("Genome IVOM: %s" % str(G))
		if args.absoluteKLI:
			#Absolute value removes direction of divergence
			windowKLI += (w*abs(math.log((w/G),2)))
		else:
			#Negative number indicates kmer depleted in window relative to genome
			#Positive number indicates enriched in window relative to genome
			#Magnitude of KLI reflects degree of difference
			windowKLI += (w*math.log((w/G),2))
			altWindowKLI += (w*math.log((w/G),10))
			print('log2KLI' + str(windowKLI))
			print('log10KLI' + str(altWindowKLI))
			#Note: Using log2 instead of log10 to accentuate variance within small range.
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

def makePicklePath(args,**kwargs):
	#Note: need to add kmer range to filename
	pathbase = os.path.basename(args.hostSeq)
	if kwargs['space'] == 'genome':
		pickleOut = os.path.join(args.tempDir,pathbase + "_kmers_" +str(args.minWordSize) + "_" + str(args.maxWordSize) + "_" + kwargs['space'] + '.p')
	else:
		pickleOut = os.path.join(args.tempDir,pathbase + "_kmers_" + str(args.minWordSize) + "_" + str(args.maxWordSize) + "_KLI_" + kwargs['space'] + "_" + str(args.windowlen) + "_increment_" + str(args.increment) + '.p')
	return pickleOut

def FDBins(data):
	#Freedman-Diaconis method for calculating optimal number of bins for dataset
	IQR = np.subtract(*np.percentile(data, [75, 25])) #Magic '*' unpacks tuple and passes two values to subtract function.
	bins = 2 * IQR * math.pow(len(data),(1.0/3.0)) #Number of bins
	bins = int(round(bins)) #Otsu needs rounded integer
	return bins

def otsu(data,optBins):
	#data is array of log10(KLI) 
	raw = data
	data = np.atleast_1d(data)
	data = data[~ np.isnan(data)]
	data = data/(max(abs(data))*-1.0) #Scale to 0-1 and make positive
	hist,binEdges = np.histogram(data,bins=optBins)
	hist = hist * 1.0 #Covert to floats
	hist_norm = hist.ravel()/hist.max() #Normalise hist to largest bin
	Q = hist_norm.cumsum()
	bins = np.arange(optBins)
	fn_min = np.inf
	thresh = -1
	for i in xrange(1,optBins): #Start from second position (1) to compare all to bin 0
		p1,p2 = np.hsplit(hist_norm,[i]) # Split normalised hist values into two brackets at bin i (bin 1 < i, bin 2 >=i)
		q1,q2 = Q[i-1], Q[optBins-1] - Q[i-1] # cum sum of bin values #!! yields zero on final bin
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

def gffFilter(rec, **kwargs):
	#Given list of feature types
	if 'feature' in kwargs.keys():
		#return GFF3 records that match type
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

def thresholdList(intervalList,threshold,args,threshCol=3,merge=True):
	tItems = [t for t in intervalList if np.log10(t[3]) >= threshold]
	sItems = sorted(tItems, key=itemgetter(0,1,2))
	anomaliesBED = pybedtools.BedTool(sItems)
	if merge:
		anomalies = anomaliesBED.merge(d=args.mergeDist, c='4,4,4', o='max,min,mean')
	else:
		anomalies = anomaliesBED
	return anomalies

def	wIndex(index, windows):
	return np.asarray([x[index] for x in windows])

def meanRangeLen(l):
	lens = list()
	for i in l:
			lens.append(i[1]-i[0])
	mean = reduce(lambda x, y: x + y, lens) / len(lens)
	return mean	

def updateHMM(smallHMM, bigThresh): #BED interval objects. A = Fine scale guide, B = Starting annotations
	'''Update annotation boundaries in B using nearest interval boundary in A.'''
	hmmBounds = dict()
	for i in smallHMM:
		if i[0] not in hmmBounds:
			hmmBounds[i[0]] = list()
			hmmBounds[i[0]].append(int(i[1]))
			hmmBounds[i[0]].append(int(i[2]))
		else:
			hmmBounds[i[0]].append(int(i[1]))
			hmmBounds[i[0]].append(int(i[2]))
	updatedBoundaries= list()
	for y in bigThresh:
		newLeft  = min(hmmBounds[y[0]], key=lambda x:abs(x-int(y[1])))
		newRight = min(hmmBounds[y[0]], key=lambda x:abs(x-int(y[2])))
		updatedBoundaries.append((y[0],newLeft,newRight))
	newAnnotations = pybedtools.BedTool(updatedBoundaries)
	return newAnnotations

def mainArgs():
	"""Process command-line arguments"""

	parser = argparse.ArgumentParser(description='Calculate all kmers in a given sequence')

	#Inputs
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
	parser.add_argument('-g',
						'--gffPath',
						type=str,
						default=None,
						help='Path to GFF file with annotations for genome being frisked.')
	
	#Outputs
	parser.add_argument('-O',
						'--outfile',
						type=str,
						help='Write KLI-IVOM bed track to this file')
	parser.add_argument('-G',
						'--gffOutfile',
						type=str,
						help='Write merged anomaly annotations to gff3.')
	parser.add_argument('-t',
						'--tempDir',
						type=str,
						default='temp',
						help='Name of temporary directory')
	parser.add_argument('--graphics',
						type=str,
						default='Summary_graphics.pdf',
						help='Name of file to print graphics to.')

	#Output options
	parser.add_argument('--mergeDist',
						type=int,
						default=0,
						help='Merge anomalies annotations within x bases of each other.')
	parser.add_argument('-f',
						'--gffFeatures',
						type=str,
						default=None,
						nargs='+',
						help='Space delimited list of feature types to report from gff')
	parser.add_argument('-r',
						'--gffRange',
						type=int,
						default=0,
						help='Report gff annotations within window of anomalous k-chores')
	
	#Core settings
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

	#Optional run settings
	parser.add_argument('-E',
						'--exitAfter',
						default=None,
						choices=[None, 'GenomeKmers', 'WindowKLI'],
						help='Exit after completing task.')
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
	parser.add_argument('--recalcWin',
						action='store_false',
						default=True,
						help='Force recalculation of KLI score for specified window length and icrement if set. \
						Default uses window KLIs from previous run.')
	parser.add_argument('-s',
						'--scaffoldsAll',
						action='store_true',
						default=False,
						help='If genomic scaffold is below minimum window size, process \
						whole scaffold as single window.')

	#Thresholding options for KLI
	parser.add_argument('--threshTypeKLI',
						default='percentile',
						choices=[None, 'percentile', 'otsu','hmm'],
						help='Options for defining non-self threshold on window KLI scores: Otsu binarisation, 2-state HMM, percentile.')
	parser.add_argument('--percentileKLI',
						type=float,
						default=99.0,
						help='Percentile at which to threshold window KLI scores. \
						By default, lower 99 percent windows are excluded.')
	parser.add_argument('-F',
						'--forceThresholdKLI',
						type=float,
						default=None,
						help='KLI score above which a window will be considered anomalous. \
						Note: Given as raw KLI, not log10(KLI).')
	
	#RIP options
	parser.add_argument('--RIP',
						action='store_true',
						default=False,
						help='Calculate and report RIP indicies for all windows + report GFF3 location \
						of RIP features, using same thresholding method options as for kmer anomalies.')
	parser.add_argument('--threshTypeRIP',
						default='percentile',
						choices=[None, 'percentile', 'otsu','hmm'],
						help='Options for defining non-self threshold on window KLI scores: Otsu binarisation, 2-state HMM, percentile.')
	
	#PCA Stuff
	parser.add_argument('-p',
						'--runProjection',
						default=None,
						choices=[None, 'PCA2', 'PCA3', 'MDS2', 'D2' ],
						help='Project anomalous windows into multidimensional space.')
	parser.add_argument('-c',
						'--culster',
						default=None,
						choices=[None, 'DBSCAN', 'KMEANS'],
						help='Attempt clustering of PCA projection.')
	parser.add_argument('-n',
						'--spikeNormal',
						action='store_true',
						default=False,
						help='Include a sampling of windows from the centre of the self population.')

	#Revise feature boundaries using HMM split track with reduced window and increment size
	parser.add_argument('--updateHMM',action='store_true',
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

	##Add options to threshold fine-scale HMM tracks
	parser.add_argument('--updateForceThreshKLI',
						type=float,
						default=None,
						help='')
	parser.add_argument('--updateCRImin',
						type=float,
						default=None,
						help='')
	parser.add_argument('--updatePImin',
						type=float,
						default=None,
						help='')
	parser.add_argument('--updateSImax',
						type=float,
						default=None,
						help='')
	parser.add_argument('--updateHmmThreshType',
						default=None,
						choices=[None, 'percentile', 'otsu'],
						help='')
	

	args = parser.parse_args()
	if args.minWordSize > args.maxWordSize:
		logging.error('[ERROR] Minimum kmer size (-c/--minKmerSize) must be less than Maximum kmer size (-k/--maxKmerSize)\n')
		sys.exit(1)
	return args

def main():
	##########################################
	########## Initial housekeeping ##########
	##########################################
	##########################################
	
	#Note: Move a lot of this off to class object
	logging.basicConfig(level=logging.INFO, format=("%(asctime)s - %(funcName)s - %(message)s"))
	args = mainArgs()
	path = args.hostSeq
	#Make path to store genome kmer calculations
	genomepickle = makePicklePath(args, space='genome')
	#Make path to store window KLI calculations
	windowsPickle = makePicklePath(args, space='window')

	#Set query sequence as self if none provided
	if not args.querySeq:
		querySeq = args.hostSeq
	else:
		querySeq = args.querySeq

	#Check temp file exists
	tempPathCheck(args)

	#Generate blank kmer dictionary
	blankMap = rangeMaps(args)

	#Read in query genome sequences as dict keyed by seq name
	selfGenome,nBlocks = getFasta(querySeq)

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
		genomeKmers = computeKmers(args, path, genomepickle, None , True, blankMap)
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
		#Initialise textfile output
		out     = args.outfile
		outPath = os.path.join(args.tempDir,out)
		handle  = open(outPath, "w")

		#List to store KLI-by-window
		allWindows = list()

		#Loop over genome to get all window KLI scores
		for seq,name,start,stop in crawlGenome(args, querySeq): #Note: coords are inclusive i.e. 1:10 = 0:9 in string
			target 		= [(name,seq)]
			windowKmers = computeKmers(args, None, None, target, False, blankMap)
			GenomeIVOM  = IvomBuild(windowKmers, args, genomeKmers, True)
			windowIVOM  = IvomBuild(windowKmers, args, genomeKmers, False)
			windowKLI   = KLI(GenomeIVOM, windowIVOM,args)
			if args.RIP:
				PI,SI,CRI 	= calcRIP(windowKmers)
				outString	= [name, str(start), str(stop), str(windowKLI),str(PI),str(SI),str(CRI)]
				allWindows.append((name,start,stop,windowKLI,PI,SI,CRI)) 
			else:
				outString	= [name, str(start), str(stop), str(windowKLI)]
				allWindows.append((name,start,stop,windowKLI)) #Mind that KLI does not get stored as scientific notation
			
			handle.write('\t'.join(outString) + '\n')
			print('\t'.join(outString))
			
		#Close text output
		handle.close()

		#Write KLI-by-window annotations to file
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

	if args.RIP:
		PI	=	wIndex(4,allWindows)
		SI	=	wIndex(5,allWindows)
		CRI	=	wIndex(6,allWindows)
		sqrCRI = CRI**2
		logCRI = np.log10(CRI)

	##########################################
	###########  Calculate and Set  ##########
	#############  KLI threshold  ############
	##########################################

	#Get KLI data for all surveyed windows
	allKLI = wIndex(3,allWindows)
	logKLI = np.log10(allKLI) #log10 transform
	
	#Calc optimal bins for KLI data
	if FDBins(logKLI) < 30:
		optBins = 30
	else:
		optBins = FDBins(logKLI)
	
	if not args.threshTypeKLI == 'hmm':
		#Check for user specified KLI threshold
		if args.forceThresholdKLI:
			KLIthreshold = np.log10(float(args.forceThresholdKLI))
			logging.info('Forcing log10(KLI) threshold = %s' % str(KLIthreshold))
		
		#Optional set threshold at percentile
		elif args.threshTypeKLI == 'percentile':
			KLIthreshold = np.percentile(logKLI, args.percentileKLI)
			logging.info('Setting threshold at %s percentile of log10(KLI)= %s' % (str(args.percentileKLI),str(KLIthreshold)))

		#Use Otsu method to set threshold
		else:
			logging.info('Calculating optimal KLI threshold by Otsu binarization.')
			if FDBins(logKLI) < 10: #Calculate optimal number of bins for logKLI set using Freedman-Diaconnis method.
				logging.warning('[WARNING] Low variance in log10(KLI) data: Review data distribution, \
								consider percentile or manual thresholding.')
				KLIthreshold = otsu(logKLI,optBins) #Attempt Otsu, though probably not a good idea.
				logging.info('Optimal log10(KLI) threshold = %s' % str(KLIthreshold))
			else:
				KLIthreshold = otsu(logKLI,optBins) #Use Otsu binarization to calc optimal log10 threshold for weird KLI scores
				logging.info('Optimal log10(KLI) threshold = %s' % str(KLIthreshold))

	##########################################
	##########  Set anomaly feature  #########
	####### boundaries using 2-state HMM  ####
	##########################################

	if args.threshTypeKLI == 'hmm':
		#Do the HMM thing
		model	=	hmm.GaussianHMM(n_components=2,covariance_type="full")
		a	=	np.asarray(allKLI)
		b	=	a[:,np.newaxis]
		#Train on all data
		model.fit(b) #Note: Ideally pass all scaffold data sequences in separately
		
		
def	hmm2BED(allWindows, dataCol, model):
	
	allIntervals = list()
	#get all scaffold names from allWindows

	#For name in allNames
		#get scaff windows from allWindows
		#sort scaffwindows by start position
		a = wIndex(dataCol,scaffWindows)
		b =	a[:,np.newaxis]
		Z = model.predict(b)
		state1 = findBaseRanges(Z, 0)
		state2 = findBaseRanges(Z, 1)

		allIntervals.append(range2interval(state1,scaffWindows,dataCol,'state1'))
		allIntervals.append(range2interval(state2,scaffWindows,dataCol,'state2'))

	allBED = pybedtools.BedTool(allIntervals)
	
	return allBED

def range2interval(rangeList, scaffoldWindows, dataCol, state):
	#deal with same start/stop
	#scaffoldWindows = (name, start, stop, KLI, othershiz)
	for block in rangeList:
		yield (scaffoldWindows[0],int(scaffoldWindows[block[0]][1]),int(scaffoldWindows[block[1]][2]),state)

def hmmBED2GFF(hmmBED):
	n = 0
	for rec in hmmBED:
		n += 1
		outstring = '\t'.join([rec[0],'frisk', rec[3], rec[1], rec[2],'.', '+', '.','ID='+ rec[3] + str(n).zfill(8)]) + '\n'
		yield outstring




	##########################################
	###########  Threshold windows  ##########
	######### and/or refine boundaries #######
	##########################################

	#Threshold and merge anomalous features.
	anomalies = thresholdList(allWindows,KLIthreshold,args)
	logging.info('Detected %s features above KLI threshold.' % str(len(anomalies)))

	#Mask 'N' blocks from anomalies
	###if args.maskN:
		###anomalies = anomalies.intersect(nBlocks)

	##########################################
	######## Write windows to GFF3 out #######
	##########################################

	#Write anomalies as GFF3 outfile
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
			logging.info('Successfully extracted %s features from within %sbp of anomaly annotations.' % (str(len(extractedFeatures)), str(args.gffRange)))
		else:
			logging.info('No features from %s detected within %s bases of anomalies.' % (args.gffPath, str(args.gffRange)))


	##########################################
	############ PCA on Anomalies ############
	##########################################
	##########################################
	#Anomalous windows by KLI threshold without merging
	anomWin  = thresholdList(allWindows,KLIthreshold,args,threshCol=3,merge=False)
	anomSeqs = getBEDSeq(selfGenome,anomWin) #yields generator object
	
	#for name,target in anomSeqs:
	#	computeKmers(args, None, None, target, False, blankMap)
	#	allKmers_win =
	#	allKmers_gen = 

	#def getKLIbyKmer(winIVOM,genIVOM):

	#A) Run PCA on window kmer counts
	#B) Run PCA on window kmer KLI scores
	
	##########################################
	######### Write Graphics to File #########
	##########################################
	##########################################

	logging.info('Writing graphics to %s' % args.graphics)

	with PdfPages(os.path.join(args.tempDir,args.graphics)) as pdf:
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
		plt.axvline(KLIthreshold, color='r', linestyle='dashed', linewidth=2)
		pdf.savefig()
		plt.close()

		plt.figure()
		plt.title('log10(KLI) Distribution Fine Bins')
		sns.set(color_codes=True)
		sns.distplot(logKLI, hist=True, bins=100, kde=False, rug=False, color="b")
		plt.axvline(KLIthreshold, color='r', linestyle='dashed', linewidth=2)
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

	##Next:
	'''	1)	Def getSeq(coords): #Feed in list of tuples / bedTool object
				for chr,start,stop in coords:

		2) 	For given sequence, get kmer counts 
			and normalise to window len. (rowname= ) Pickle object.
		
		3) 	Driving Kmer report. 
			For subset of windows rank kmers by mean IKW-KLI, print as report.

		4)	Run projection (PCA, MDS). Save image.

		5)	Run clustering on Projection. (DBSCAN, other methods)

		6)	Save image to file [PCA,PCA with clusters].
	
		7)	Report coords for anomalies with category label in type.

		8)	Function to extract super normal windows to spike PCA.

		9)	Option to recycle KLI track previously calculated per scaffold >> continue to gff 
			reports and PCA

		10) Build RIP into output: * have made 'calcRIP(windowKmers)'

		11)	Mask N-blocks from final merged anoms
		
		12)	 Option to return high self-scoring blocks (i.e All windows below KLI-OTSU threshold)
			This feature to be used for mitochondrial ID test case: Train self as Mt Genome, query 
			assembly and annotate Mt-like regions.

		13)	Use hmmlearn to classifiy 1000bp window 4mer counts as self or non-self 
			using MultinomialHMM trained on otsu thresholded 5000bp windows 
			(use pyBEDTools to retrieve 1000bp windows contained within 5000bp windows)

	'''
	# Trace
	logging.info('Finished!')

if __name__ == '__main__':
	main()