#draft_histograms.py

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import seaborn as sns
import math

fname = "XAC_rawTrack.bed"
with open(fname) as f:
	lines = [float(line.rstrip('\n').split()[3]) for line in f]

rawData = np.asarray(lines)
logData = np.log10(rawData)
FDBins(logData)
otsu(logData,100)


plt.title('KLI Distribution')
sns.set(color_codes=True)
sns.distplot(logData, hist=True, bins=100, kde=False, rug=False, color="b")
plt.axvline(thresh, color='r', linestyle='dashed', linewidth=2)
plt.show()

####

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
	print('OTSU selected bin %s as threshold position' % str(thresh))
	#Convert bact to log10(KLI)
	otsuNum = binEdges[thresh] * max(abs(raw)) * -1.0
	return otsuNum