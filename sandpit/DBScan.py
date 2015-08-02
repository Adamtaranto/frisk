#DBScan
import time
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN

dbscan = DBSCAN(eps=.2)

t0 = time.time()
dbscan.fit(X)
t1 = time.time()
if hasattr(dbscan, 'labels_'):
	y_pred = dbscan.labels_.astype(np.int)
else:
	y_pred = dbscan.predict(X)

# plot
colors = np.array([x for x in 'bgrcmykbgrcmykbgrcmykbgrcmyk'])
colors = np.hstack([colors] * 20)

plt.title('dbscan', size=18)
plt.scatter(X[:, 0], X[:, 1], color=colors[y_pred].tolist(), s=10)

if hasattr(dbscan, 'cluster_centers_'):
	centers = dbscan.cluster_centers_
	center_colors = colors[:len(centers)]
	plt.scatter(centers[:, 0], centers[:, 1], s=100, c=center_colors)
plt.xlim(-2, 2)
plt.ylim(-2, 2)
plt.xticks(())
plt.yticks(())
plt.text(.99, .01, ('%.2fs' % (t1 - t0)).lstrip('0'),
		 transform=plt.gca().transAxes, size=15,
		 horizontalalignment='right')
plt.show()
