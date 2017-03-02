# Copyright (c) 2017 Kevin Murray <kdmfoss@gmail.com>
# Copyright (c) 2017 Adam Taranto
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import sys
import logging


def cli_log(name="frisk-cli"):
    # prints log to stdout and also saves to specified log file
    logger = logging.getLogger(name)
    fmtclass = logging.Formatter
    if sys.stderr.isatty():
        try:
            import coloredlogs
            fmtclass = coloredlogs.ColoredFormatter
        except ImportError:
            pass
    formatter = fmtclass(fmt="%(asctime)s: %(message)s", datefmt="%H:%M:%S")
    stderrh = logging.StreamHandler()
    stderrh.setFormatter(formatter)
    logger.addHandler(stderrh)
    logger.setLevel(logging.DEBUG)
    return logger
