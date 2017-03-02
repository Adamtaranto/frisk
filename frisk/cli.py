# Copyright (c) 2017 Kevin Murray <kdmfoss@gmail.com>
# Copyright (c) 2017 Adam Taranto
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import argparse

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions


def kli_args():
    p = argparse.ArgumentParser(
        description="Calculate KLD between genome windows and the whole genome",
        prog="frisk-kli",
    )

    p.add_argument(
        "-V,--version",
        action="version",
        version="frisk -- version {}".format(__version__)
    )
    p.add_argument(
        "-s,--self-genome",
        type=str,
        required=True,
        help="The input self genome sequences (as FASTA)"
    )
    p.add_argument(
        "-n,--non-self-genome",
        type=str,
        required=False,
        help="The input non-self genome sequences (as FASTA)"
    )
    p.add_argument(
        "-o,--outfile",
        type=str,
        help="Write KLD bed track to this file",
    )
    p.add_argument(
        "-k,--kmer-size",
        type=int,
        default=4,
        help="Length of k-mers",
    )
    p.add_argument(
        "-w,--window-size",
        type=int,
        default=5000,
        help="Length of genome sequence windows",
    )

    mode = p.add_mutually_exclusive_group()
    mode.add_argument(
        "-C,--cre",
        action="store_true",
        help="Use the Conditional Relative Entropy measure"
    )
    mode.add_argument(
        "-I,--ivom",
        action="store_true",
        help="Use the Conditional Relative Entropy measure"
    )

    args = p.parse_args()
    return args


def kli_main():
    args = kli_args()
    return args
