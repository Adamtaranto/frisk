# cython: profile=False

import numpy as np
cimport numpy as cnp


def build_kmer_vec(min_k, max_k, dtype=np.uint64, alphabet="ACGT"):
    """Creates a dictionary of k-mer count vectors

    The dictionary returned is of the form {1: np.ndarray(), 2: np.ndarray()}.
    """
    kmer_vecs = {}
    for k in range(min_k, max_k + 1):
        n_kmers = len(alphabet) ** k
        kmer_vecs[k] = np.zeros(n_kmers, dtype=dtype)
    return kmer_vecs


def kmerhash(str kmer not None):
    """Hash a k-mer, in 2-bit encoding."""
    cdef unsigned long long hash = 0
    cdef unsigned long long n
    if len(kmer) > 32:
        raise ValueError("Only K <= 32 supported")
    for nt in kmer:
        n = (ord(nt) & 0x06) >>1
        n = n ^ (n >> 1)
        hash <<= 2
        hash |= n
    return hash


def hash_seq(str seq not None,
             dict kmer_vec not None):
    cdef int min_k = min(kmer_vec.keys())
    cdef int max_k = max(kmer_vec.keys())
    if len(seq) < max_k:
        raise ValueError("Sequence too short: max(k) > len(seq)!")
    for k in range(min_k, max_k + 1):
        for start in range(len(seq) - k + 1):
            kmer = seq[start:start + k]
            hash = kmerhash(kmer)
            kmer_vec[k][hash] += 1


def kmer_freqs(dict kmer_counts not None):
    kmer_freqs = {}
    for k in kmer_counts:
        kc = kmer_counts[k].astype(np.float)
        kmer_freqs[k] = kc / kc.sum()
    return kmer_freqs


def ivom(dict kmer_freqs not None, int alphabet_sz=4):
    cdef int min_k = min(kmer_freqs.keys())
    cdef int max_k = max(kmer_freqs.keys())
    cdef unsigned long long subkmer = 0
    cdef float count = 0

    maxmer_ivoms = np.zeros(kmer_freqs[max_k].shape, dtype=np.float)
    sums = {k: v.sum() for k, v in kmer_freqs.items()}

    weights = np.zeros(max_k - min_k + 1, dtype=np.float)
    probs = np.zeros(max_k - min_k + 1, dtype=np.float)
    maxmer_freqs = kmer_freqs[max_k]
    for kmer in range(len(maxmer_freqs)):
        if maxmer_freqs[kmer] == 0:
            continue
        for x in range(min_k, max_k + 1):
            subkmer = kmer >> (max_k - x) * 2
            count = kmer_freqs[x][subkmer]
            weights[x-1] = count * alphabet_sz ** x
            probs[x-1] = count / (sums[x] * 2)

        weights /= weights.cumsum()

        last_ivom = 0.0
        for x in range(min_k - 1, max_k):
            last_ivom = weights[x] * probs[x] + (1 - weights[x] * last_ivom)
        maxmer_ivoms[kmer] = last_ivom
    return maxmer_ivoms / maxmer_ivoms.sum()


def kli(cnp.ndarray gen_ivom, cnp.ndarray window_ivom):
    w = window_ivom.astype(np.float)
    g = gen_ivom.astype(np.float)
    return np.nansum(w * np.log2(w/g))