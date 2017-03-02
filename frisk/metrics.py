# Copyright (c) 2017 Kevin Murray <kdmfoss@gmail.com>
# Copyright (c) 2017 Adam Taranto
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

'''Anomaly Detection algorithms that underly Frisk.

This module implements two Anomaly Detection algorithms: CRE and IVOM.
'''

import pymer
import numpy as np

from .funcs import each_window


class BaseCalc(object):
    '''Base class and common utilities for CRE and IVOM methods'''
    counter_type = pymer.BaseCounter

    def __init__(self, k, window_size=5000, window_overlap=0.5):
        self.k = k
        self.genome = self.make_counter()
        self.window_size = window_size
        self.window_overlap = window_overlap

    def make_counter(self):
        return self.counter_type(self.k)

    def load_genome_file(self, filename):
        self.genome.consume_file(filename)

    def load_genome_seq(self, seq):
        self.genome.consume(seq)

    def count_window(self, seq):
        counter = self.make_counter()
        counter.consume(seq)

    def window_score(self, window):
        return 0

    def window_scores(self, seq):
        scores = []
        for window in each_window(seq, self.window_size, self.window_overlap):
            scores.append(self.window_score(window))
        return scores


class CREAnomalyDetector(BaseCalc):
    '''Conditional Relative Entropy of Transition Frequencies'''
    counter_type = pymer.TransitionKmerCounter

    def window_score(self, window):
        counts = self.count_window(window)
        P_g = self.genome.transitions
        P_w = counts.transitions
        s_w = counts.stem_frequencies
        cre = np.sum(s_w * np.nansum(P_w * np.log2(P_w/P_g), 1))
        return cre


class IVOMAnomalyDetector(BaseCalc):

    def ivom(self, sequence):
        '''Calculates the IVOM k-mer frequencies of a sequence'''
        counters = {k: self.make_counter(k) for k in range(1, self.k+1)}

        for counter in counters.values():
            counter.consume(sequence)

        maxmers = counters[self.k]

        for k, counter in sorted(counters.items()):
            pass
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
        '''
