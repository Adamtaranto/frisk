#!/usr/bin/env python
from __future__ import division, print_function, absolute_import
import numpy as np
import sys


def gc2p(gc):
    '''Takes percentage GC and returns base probs for ACGT'''
    gc = float(gc)
    if gc < 0. or gc > 1.:
        raise ValueError("Invalid proportion {}".format(gc))
    a = t = (1. - gc) / 2.
    g = c = gc / 2.
    return [a, c, g, t]


def randseq(length, gc=0.5):
    NTs = np.array(list("ACGT"))
    seq = ''.join(np.random.choice(NTs, length, p=gc2p(gc)))
    return seq


def parse_compositiion_cfg(fpath):
    '''expects a file of:

        name    length  %gc

    e.g.:

        chr1    100000  0.5
        chr1    10000   0.9
        chr1    100000  0.5

    Can be more than one line per seq, in which case they are concatenated.
    '''
    parts = {}
    with open(fpath) as hdl:
        for line in hdl:
            if line.lower().startswith("name\t"):
                # someone left the header in!
                continue
            name, length, gc = line.split()
            length = int(length)
            gc = float(gc)
            if name in parts:
                parts[name].append((length, gc))
            else:
                parts[name] = [(length, gc), ]
    return parts


def gen_seqs(cfg):
    seqs = []
    for seqname, sections in cfg.items():
        seq = []
        for length, gc in sections:
            seq.append(randseq(length, gc=gc))
        seq = ''.join(seq)
        seqs.append((seqname, seq))
    return seqs


def write_seqs(seqs, file=sys.stdout, width=80):
    for seqname, seq in seqs:
        print('>', seqname, sep='', file=file)
        for i in range(0, len(seq), width):
            print(seq[i:i+width], file=file)


if __name__ == '__main__':
    from docopt import docopt

    CLI = '''
    USAGE:
        seqgen.py <cfgfile>

    <cfgfile> is a plain text file, with whitespace-separated columns for a
    sequence name, length and gc content. Sequences can have many sections,
    e.g.:

        chr1    100000  0.5
        chr1    10000   0.9
        chr1    100000  0.5
    '''

    args = docopt(CLI)

    cfg = parse_compositiion_cfg(args['<cfgfile>'])
    seqs = gen_seqs(cfg)
    write_seqs(seqs, file=sys.stdout)
