#!/usr/bin/env python

def drivingKmers(seqList):
	'''Take sequences from cluster return mean and variance 
	of IVOM KLI for each Kmer. For each k, rank +/- kmers 
	by KLI. Indicates strong drivers of divergence.'''

	#dictionary of all kmers by k, each kmer = list()
	clusterKLIs = copy.deepcopy(blankKmers)
	#set all values in clusterKLI to empty list
	#????????????

	for seq in seqList:
		counts = getKmers(seq)
		seqIVOMbyKmer = getIVOM(counts,KLIbyKmerMode) #By kmer mode to return dictionary of all kmer IVOM
		genomeIVOMbyKmer = getIVOM(counts, genomecounts, genomemode, KLIbyKmerMode)
		for k in seqIVOMbyKmer:
			for kmer in seqIVOMbyKmer[k]:
				kmerKLI = getKLI(seqIVOMbyKmer[k][kmer], genomeIVOMbyKmer[k][kmer])
				clusterKLIs[k][kmer].append(kmerKLI)

	kmerMeanVars = list()

	for k in clusterKLIs:
		kmerMeanVars[k] = list()
		for kmer in clusterKLIs[k]:
			kmerMeanVars[k].append(tuple(len(kmer),kmer,mean(clusterKLIs[k][kmer]),var(clusterKLIs[k][kmer])))

	return kmerMeanVars

def dkReport(clusterID, kmerMeanVars, reportCap=100, outputFile):

	#open outfile for writing by append:
	outfile.write("Report for " + str(clusterID))

	def formatString(x):
		'''tab delimited kmer + mean + var'''
		counter = 0
		for record in x:
			string = str(record[2]),str(record[3]),str(record[4]),"\n"
			if counter < reportCap:
				yield string
				counter += 1
			else:
			continue

	for k in kmerMeanVars:
		posList = list()
		neglist = list()
		for kmer in kmerMeanVars[k]:
			if kmer[3] < 0:
				neglist.append(kmer)
			elif kmer[3] >= 0:
				posList.append(kmer)
		
		outfile.write("poskmers for k=" + str(k+1) + "\n")
		posSorted = sort(poslist) #sort on mean then var ???
		outfile.write(formatString(posSorted))

		outfile.write("negkmers for k=" + str(k+1) + "\n")
		negSorted = sort(neglist) #sort on abs(mean) then var ???
		outfile.write(formatString(negSorted))


def getPCACounts(seq, krange, args):
	'''Take sequence and k-range. Generate counts normalised 
	to number of kmers observed for each k.'''
	if not krange:
		kmin = args.kmin
		kmax = args.kmax
	else:
		kmin = krange[0]
		kmax = krange[1]

	kDict = getKmers(seq)

	#Normalise to observed k
	for k in kDict:
		kTotal = sum(kDict[k])
		for kmer in kDict[k]:
			kDict[k][kmer] = float(kDict[k][kmer])/kTotal

	#Make list for PCA
	pcaData = list()
	for k in xrange(kmin:kmax+1):
		for kmer in kDict[k]:
			pcaData.append(kmer)

	return pcaData



def extractWindow(rangeList,fastaObject):
	'''Give scaffold ID, start and stop positions, fetch 
		sequence from fastaObject'''
	for x in rangeList:
		seq = extract(scaff+start+stop)
		yield seq