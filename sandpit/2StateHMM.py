#git clone git@github.com:hmmlearn/hmmlearn.git
#cd hmmlearn
#python setup.py install
#from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import seaborn as sns
from hmmlearn import hmm
import itertools
from operator import itemgetter
import pybedtools

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
	return ranges

def meanRangeLen(l):
	lens = list()
	for i in l:
			lens.append(i[1]-i[0])
	mean = reduce(lambda x, y: x + y, lens) / len(lens)
	return mean	

def updateHMM(smallHMM, bigThresh):
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

#Getdata
fname = "XAC_rawTrack.bed"
with open(fname) as f:
	lines = [float(line.rstrip('\n').split()[3]) for line in f]

rawData = np.asarray(lines)
logData = np.log10(rawData)

#Do the HMM thing
model = hmm.GaussianHMM(n_components=2,covariance_type="full")
x = lines
a=np.asarray(x)
b=a[:,np.newaxis]

#Train on all data
model.fit(b) #Note: figure out way to add all scaffolds. Concat?
#Calc per scaffold
Z = model.predict(b) #Run one at a time

#CalcRanges
state1 = findBaseRanges(Z, 0)
state2 = findBaseRanges(Z, 1)

if meanRangeLen(state1) > meanRangeLen(state2):
	selfState = state1
	nonSelfState = state2
else:
	selfState = state2
	nonSelfState = state1

