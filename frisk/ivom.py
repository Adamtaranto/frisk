import numpy as np

from pymer import ExactKmerCounter

def ivom(str sequence not None, int maxk, int mink, int alphabet_sz=4):
    '''Calculates the Interpolated Variable-order Motif probablities

    Arguments
    ---------

    - `kmer_freqs`: dict of {klen: np.array(4**k), } for each kmer length
    '''
    cdef unsigned long long subkmer = 0
    cdef double count = 0

    # Preallocate the vector for the maxmer ivom probs
    maxmer_ivoms = np.zeros(kmer_freqs[max_k].shape, dtype=np.float)
    # Sum of each kmer counts for each k
    sums = {k: counts.sum() for k, counts in kmer_freqs.items()}

    # w and p vectors for each maxmer in loop below. Preallocated outside loop
    # for efficiency
    weights = np.zeros(max_k - min_k + 1, dtype=np.float)
    probs = np.zeros(max_k - min_k + 1, dtype=np.float)

    # frequencies of maxmers. A vector of 4**k for maximal k
    maxmer_freqs = kmer_freqs[max_k]
    for kmer in range(len(maxmer_freqs)):
        if maxmer_freqs[kmer] == 0:
            continue
        for k in range(min_k, max_k + 1):
            subkmer = kmer >> ((max_k - k) * 2)
            count = kmer_freqs[k][subkmer]
            weights[k-1] = count * alphabet_sz ** k
            probs[k-1] = count / (sums[k] * 2)

        weights /= weights.cumsum()

        last_ivom = 0.0
        for k in range(min_k - 1, max_k):
            last_ivom = weights[k] * probs[k] + (1 - weights[k] * last_ivom)
        maxmer_ivoms[kmer] = last_ivom
    return maxmer_ivoms / maxmer_ivoms.sum()
