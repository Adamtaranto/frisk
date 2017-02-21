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


class CRECalc(BaseCalc):
    '''Conditional Relative Entropy of Transition Frequencies
    '''
    counter_type = pymer.TransitionKmerCounter

    def window_score(self, window):
        counts = self.count_window(window)
        P_g = self.genome.transitions
        P_w = counts.transitions
        s_w = counts.stem_frequencies
        cre = np.sum(s_w * np.nansum(P_w * np.log2(P_w/P_g), 1))
        return cre

