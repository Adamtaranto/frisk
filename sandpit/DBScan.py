#DBScan
import time
import numpy as np
from sklearn.cluster import DBSCAN
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from sklearn.decomposition import PCA
import copy
import pickle
import itertools
from operator import itemgetter 

X = pickle.load( open( "anomCountsSymProp_fixLenNorm", "rb" ) )
pca_X = PCA(n_components=3)
X = pca_X.fit(X).transform(X)

dbscan = DBSCAN(eps=0.02,min_samples=10).fit(X)

if hasattr(dbscan, 'labels_'):
	y_pred = dbscan.labels_.astype(np.int)
else:
	y_pred = dbscan.predict(X)

#Colours and markers
colors = np.array([x for x in 'bgrcmykbgrcmykbgrcmykbgrcmyk'])
colors = np.hstack([colors] * 20)

markers = np.array([x for x in '+o^8sD*+o^8sD*+o^8sD*+o^8sD*'])
markers = np.hstack([markers] * 20)


#List dbscan classes
y_class = np.unique(y_pred)
#List to hold scatterplots for each class
scatter_class = np.unique(y_pred).tolist()

#List to hold generated class names
names_class = np.unique(y_pred).tolist()
count_class = 0
for i in range(0,len(names_class)):
	count_class += 1
	if names_class[i] == -1:
		names_class[i] = 'Noise'
	else:
		names_class[i] = 'Class_' + str(i)

fig = plt.figure() #figsize=(8,8)
ax = fig.add_subplot(111)

#Generate scatterplots for each class using pca coords
for i in y_class:
	X_class = X[np.where(y_pred == i)]
	scatter_idx = np.where(scatter_class == i)[0][0]
	scatter_class[scatter_idx] = plt.scatter(X_class[:,0], X_class[:,1], color=colors[i], marker=markers[i], s=20)
	##scatter_class[scatter_idx] = ax.plot(X_class[:,0], X_class[:,1], X_class[:,2], markers[i], color=colors[i], markersize=8, alpha=0.4,label=names_class[scatter_idx])

#Plot
#axisLabels = pcaAxisLabels(pca_X, kmin=1, kmax=6)

plt.title('PCA with DBSCAN Clusters', size=18)
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(scatter_class, names_class,loc='center left', bbox_to_anchor=(1, 0.5))
plt.xlabel(axisLabels[0])
plt.ylabel(axisLabels[1])
#plt.rcParams['legend.fontsize'] = 10
plt.show()

#####3D in plotly

import plotly.plotly as py
from plotly.graph_objs import *

trace1 = Scatter3d(
	x=X_3[:,0],
	y=X_3[:,1], 
	z=X_3[:,2],
	surfaceaxis = "-1",
	surfacecolor = "red",
	hoverinfo = "name+text" 
	mode='markers',
	name='kmer Anomalies', #Class label
	text=anomLabels,
	marker=Marker(
		symbol="diamond-open",
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
	margin=Margin(l=0,r=0,b=0,t=0),
	title="Draft-Kmer-Anomaly-PCA",
	scene=Scene(
		xaxis=XAxis(title="X"),
		yaxis=YAxis(title="Y"),
		zaxis=ZAxis(title="Z")
		)
	)

fig = Figure(data=data, layout=layout)
plot_url = py.plot(fig, filename='Draft-Kmer-Anomaly-PCA')



"circle"
"circle-open"
"square"
"square-open"
"diamond"
"diamond-open"
"cross"
"x-open" 
"x-open-dot"
"triangle-up-dot"
"triangle-down"
"triangle-down-open"
"triangle-down-open-dot"
"triangle-se-open-dot"
"triangle-sw-open"
"hexagon-open"
"octagon"
"octagon-open"
"star"

#########################

## Subset anomaly list
from operator import itemgetter 

for i in y_class: # [-1,0,1,2,3,4]
	class_idx = np.where(y_pred == i) # Positions in array where value is in class
	names_class[np.where(y_class == i)[0][0]] # Name of class
	anoms_in_class = itemgetter(*class_idx.tolist())(anom_tuple_list) #if anom master set is in list
	# anoms_in_class = anom_tuple_list[class_idx]
	#Call export to GFF, provide class info


#plt.scatter(X[:, 0], X[:, 1], color=colors[y_pred].tolist(),marker=markers[y_pred].tolist(), s=10)

#for i in range(0, X.shape[0]):
#	if y_pred[i] == 0:
#		c1 = plt.scatter(X[i,0],X[i,1],c='r',marker='+',s=100)
#	elif y_pred[i] == 1:
#		c1 = plt.scatter(X[i,0],X[i,1],c='g',marker='o',s=100)
#	elif y_pred[i] == 2:
#		c3 = plt.scatter(X[i,0],X[i,1],c='c',marker='o',s=100)
#	elif y_pred[i] == -1:
#		c4 = plt.scatter(X[i,0],X[i,1],c='b',marker='*',s=100)

