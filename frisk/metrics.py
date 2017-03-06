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

from .funcs import each_window, each_seqeunce

# Ignore floating point errors, as we use nansum
np.seterr(all="ignore")


class BaseCalc(object):
    '''Base class and common utilities for CRE and IVOM methods'''
    counter_type = pymer.BaseCounter

    def __init__(self, k, window_size=5000, window_overlap=0.5):
        self.k = k
        self.genome = self._make_counter()
        self.window_size = window_size
        self.window_overlap = window_overlap

    def _make_counter(self):
        return self.counter_type(self.k)

    def _count_window(self, seq):
        counter = self._make_counter()
        counter.consume(seq)
        return counter

    def load_genome_file(self, filename):
        self.genome.consume_file(filename)

    def load_genome_seq(self, seq):
        self.genome.consume(seq)

    def window_score(self, window):
        return 0

    def window_scores_seq(self, seq):
        for start, stop, window in each_window(seq, size=self.window_size,
                                               offset=self.window_overlap,
                                               coordinates=True):
            yield (start, stop, self.window_score(window))

    def window_scores_file(self, filename):
        for seq in each_seqeunce(filename):
            for s, e, w in self.window_scores_seq(seq.sequence):
                yield (seq.name, s, e, w)


class CREAnomalyDetector(BaseCalc):
    '''Conditional Relative Entropy of Transition Frequencies'''
    counter_type = pymer.TransitionKmerCounter

    def window_score(self, window):
        counts = self._count_window(window)
        P_g = self.genome.transitions
        P_w = counts.transitions
        s_w = counts.stem_frequencies
        cre = np.sum(s_w * np.nansum(P_w * np.log2(P_w/P_g), 1))
        return cre
