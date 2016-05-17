import numpy as np


def kmer_freqs(kmer_counts):
    '''Turn kmer count dict into freqs'''
    kmer_freqs = {}
    for k in kmer_counts:
        kc = kmer_counts[k].astype(np.float)
        kmer_freqs[k] = kc / kc.sum()
    return kmer_freqs


def build_kmer_vec(min_k, max_k, dtype=np.uint64, alphabet="ACGT"):
    """Creates a dictionary of k-mer count vectors

    The dictionary returned is of the form {1: np.ndarray(), 2: np.ndarray()}.
    """
    kmer_vecs = {}
    for k in range(min_k, max_k + 1):
        n_kmers = len(alphabet) ** k
        kmer_vecs[k] = np.zeros(n_kmers, dtype=dtype)
    return kmer_vecs


def kli(genome_ivom, window_ivom):
    w = window_ivom.astype(np.float)
    g = genome_ivom.astype(np.float)

    #kli = [w[i] * log2(w[i] / g[i]) for i in len(w) if g[i] > 0 and w[i] > 0]

    kli = w * np.log2(w / g)
    # nansum sums over non-NAN values. NaN values are introduced by div by zero
    # in the above kli score calculation. This is the same as skipping
    return np.nansum(kli)
