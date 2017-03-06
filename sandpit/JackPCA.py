#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from sklearn.decomposition import PCA
from sklearn import datasets
import copy
import pickle
import itertools

X = pickle.load( open( "anomCountsSymProp_fixLenNorm", "rb" ) )

#slice out 2-5mers from pre-flattened counts
#z = copy.deepcopy(X[:,4:1365])
#z_label = copy.deepcopy(Y[:,4:1365])

#Symetrical + proportional
pca_X = PCA(n_components=2)
X_fit = pca_X.fit(X).transform(X)

plt.figure(figsize = (35, 20))
plt.scatter(X_fit[:, 0], X_fit[:, 1], s=200)
plt.show()


##3D!!!
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import proj3d

z = X

pca3 = PCA(n_components=3)
pca3.explained_variance_ratio_

X_3 = pca3.fit(z).transform(z)

fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(111, projection='3d')
plt.rcParams['legend.fontsize'] = 10
ax.plot(X_3[:,0], X_3[:,1], X_3[:,2],'o', markersize=8, color='blue', alpha=0.4, label='class1')

#ax.plot(class2_sample[0,:], class2_sample[1,:], class2_sample[2,:],
#        '^', markersize=8, alpha=0.5, color='red', label='class2')

plt.title('Samples for class 1')
ax.legend(loc='upper right')

ax.view_init(elev=0, azim=-90)

plt.show()


#########################
#Find driving variables
LETTERS = ('A', 'T', 'G', 'C')

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

def rangeMaps(kmin,kmax):
	#Calls prepareMaps to write list of kmer dictionaries for range kmin to kmax
	kMin    = kmin
	kMax    = kmax
	maps    = []
	#Prepare possible kmer combinations for each len k
	for i in xrange(kMin, kMax + 1):
		maps.append(prepareMaps(0, i, ['']))
	return maps

def pcaAxisLabels(pca_X, kmin=1, kmax=6):
	'Generate axis labels for PCA. Format: "PC1: AT (90.01%)"'
	blankmap 	= rangeMaps(kmin,kmax)
	keylist 	= []
	axisLabels 	= []
	idx 		= 0 #PC index
	for x in blankmap:
		for y in x.keys():
			keylist.append(y)
	for x in pca_X.explained_variance_ratio_:
		labelString = 'PC%s: %s (%s)' % (str(idx+1),str(keylist[np.argmax(pca_X.components_[idx])]),str(x*100))
		axisLabels.append(labelString)
		idx += 1
	return axisLabels

blankmap = rangeMaps(1,6)

keylist = []

for x in blankmap:
	for y in x.keys():
		keylist.append(y)


idx = 0 #PC index
cv = 0 #cumulative variance

for x in pca_X.explained_variance_ratio_:
	print 'PC:',idx+1
	print keylist[np.argmax(pca_X.components_[idx])]
	idx += 1
	cv += x
	print(cv)

####################



##3D with plotly
#python -c "import plotly; plotly.tools.set_credentials_file(username='adamtaranto', api_key='4fit77pb6l')"

import plotly.plotly as py
from plotly.graph_objs import *

trace1 = Scatter3d(
	x=X_3[:,0],
	y=X_3[:,1], 
	z=X_3[:,2],
	mode='markers',
	name='kmer Anomalies',
	text=anomLabels,
	marker=Marker(
		size=8,
		line=Line(
			color='rgba(217, 217, 217, 0.14)',
			width=0.00
		),
		opacity=0.6
	)
)
data = Data([trace1])
layout = Layout(
	margin=Margin(
		l=0,
		r=0,
		b=0,
		t=0
	)
)
fig = Figure(data=data, layout=layout)
plot_url = py.plot(fig, filename='Draft-Kmer-Anomaly-PCA')













"""
=========================================================
PCA example with Iris Data-set
=========================================================

"""
#iris = datasets.load_iris()
#X = iris.data
#y = iris.target

#X is my data as a numpy array, each row was a flattened image (1D) and each column was one of the pixel values from the image
#X = np.array(tag_roi_flat)
#y are the labels for the different images (tag 1, 2, and 3
#y = np.array([2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 3, 3, 3])

np.random.seed(5)

centers = [[1, 1], [-1, -1], [1, -1]]
iris = datasets.load_iris()
X = iris.data
y = iris.target

fig = plt.figure(1, figsize=(4, 3))
plt.clf()
ax = Axes3D(fig, rect=[0, 0, .95, 1], elev=48, azim=134)

plt.cla()
pca = PCA(n_components=3)
pca.fit(X)
X = pca.transform(X)

for name, label in [('Setosa', 0), ('Versicolour', 1), ('Virginica', 2)]:
	ax.text3D(X[y == label, 0].mean(),
			  X[y == label, 1].mean() + 1.5,
			  X[y == label, 2].mean(), name,
			  horizontalalignment='center',
			  bbox=dict(alpha=.5, edgecolor='w', facecolor='w'))
# Reorder the labels to have colors matching the cluster results
y = np.choose(y, [1, 2, 0]).astype(np.float)
ax.scatter(X[:, 0], X[:, 1], X[:, 2], c=y, cmap=plt.cm.spectral)

x_surf = [X[:, 0].min(), X[:, 0].max(),
		  X[:, 0].min(), X[:, 0].max()]
y_surf = [X[:, 0].max(), X[:, 0].max(),
		  X[:, 0].min(), X[:, 0].min()]
x_surf = np.array(x_surf)
y_surf = np.array(y_surf)
v0 = pca.transform(pca.components_[0])
v0 /= v0[-1]
v1 = pca.transform(pca.components_[1])
v1 /= v1[-1]

ax.w_xaxis.set_ticklabels([])
ax.w_yaxis.set_ticklabels([])
ax.w_zaxis.set_ticklabels([])

plt.show()