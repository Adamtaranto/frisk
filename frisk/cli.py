# Copyright (c) 2017 Kevin Murray <kdmfoss@gmail.com>
# Copyright (c) 2017 Adam Taranto
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

from sys import stdin, stdout, stderr, exit
import argparse

import numpy as np

import frisk
from .metrics import (
    CREAnomalyDetector,
)
from .log import CLI_LOG as LOG
from .funcs import write_bed


def kld_args():
    p = argparse.ArgumentParser(
        description="Calculate KLD between genome windows and the whole genome",
        prog="frisk-kld",
    )

    p.add_argument(
        "-V", "--version", action="version",
        version="frisk -- version {}".format(frisk.__version__)
    )
    p.add_argument(
        "-s", "--self-genome", type=str, required=True,
        help="The input self genome sequences (as FASTA)"
    )
    p.add_argument(
        "-n", "--non-self-genome", type=str, required=False,
        help="The input non-self genome sequences (as FASTA)"
    )
    p.add_argument(
        "-o", "--outfile", type=str,
        help="Write KLD bed track to this file",
    )
    p.add_argument(
        "-k", "--kmer-size", type=int, default=4,
        help="Length of k-mers",
    )
    p.add_argument(
        "-w", "--window-size", type=int, default=5000,
        help="Length of genome sequence windows",
    )
    p.add_argument(
        "-m", "--mode", choices=("CRE", "IVOM"), default="CRE",
        help="KLD Calculation mode to use"
    )

    args = p.parse_args()
    return args


def kld_main():
    args = kld_args()
    self_file = args.self_genome
    nonself_file= args.self_genome  # FIXME: args.non_self_genome if supplied

    # Instantiate kld calculator
    calc_classes = {
        "CRE": CREAnomalyDetector,
        "IVOM": NotImplementedError,
    }
    kldcalc = calc_classes[args.mode](args.kmer_size, window_size=args.window_size)

    LOG.info("Loading \"self\" genome")
    kldcalc.load_genome_file(self_file)
    LOG.info("\"Self\" genome loaded")

    LOG.info("Calculating window scores")
    with open(args.outfile, "w") as bedfh:
        for i, window in enumerate(kldcalc.window_scores_file(nonself_file)):
            write_bed(bedfh, *window)
            if i % 100000 == 0:
                LOG.info("Processed {} windows".format(i))
        LOG.info("Done, processed {} windows".format(i))

    LOG.info("Finished")
