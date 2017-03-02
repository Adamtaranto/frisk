# Copyright (c) 2017 Kevin Murray <kdmfoss@gmail.com>
# Copyright (c) 2017 Adam Taranto
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.


def each_window(seq, size=5000, offset=0.5, coordinates=False):
    # FIXME: Make this deal with the leftovers at the end
    # FIXME: return coordinates if cordinates == True
    if offset < 1:
        offset *= size
    offset = int(offset)
    for start in range(0, len(seq), offset):
        if start + size + offset > len(seq):
            size = len(seq) - start
        yield seq[start:start+size]
