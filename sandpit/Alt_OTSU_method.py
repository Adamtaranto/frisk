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