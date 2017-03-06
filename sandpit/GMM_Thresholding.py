#!/usr/bin/env python
#GMM_Thresholding.py

import numpy as np
from sklearn import mixture
import matplotlib.pyplot as plt

fname = "just_data.tab"
with open(fname) as f:
	lines = [float(line.rstrip('\n')) for line in f]

rawData = np.asarray(lines)
logData = np.log10(rawData)
X = map(list,zip(*[iter(logData)]*1))



def fit_samples(samples):
	gmix = mixture.GMM(n_components=3, covariance_type='full')
	gmix.fit(samples)
	print gmix.means_
	#colors = ['r' if i==0 else 'g' for i in gmix.predict(samples)]
	#ax = plt.gca()
	#????
	#plt.show()

gmix = fit_samples(X)

myData = logData
hist,binEdges = np.histogram(myData,bins=30)
plt.bar(binEdges[:-1], hist, width = 0.1)
plt.axvline(gmix.means_[0][0], color='r', linestyle='dashed', linewidth=2)
plt.axvline(gmix.means_[1][0], color='r', linestyle='dashed', linewidth=2)
plt.axvline(gmix.means_[2][0], color='r', linestyle='dashed', linewidth=2)
plt.xlim(min(binEdges), max(binEdges))
plt.show()

# fit models with 1-10 components
N = np.arange(1, 11)
models = [None for i in range(len(N))]

for i in range(len(N)):
    models[i] = mixture.GMM(N[i]).fit(X)

# compute the AIC and the BIC
AIC = list()
for m in models:
	AIC.append(m.aic(X))

AIC.index(min(AIC))


