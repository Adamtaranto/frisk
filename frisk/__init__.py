# Copyright (c) 2017 Kevin Murray <kdmfoss@gmail.com>, 2017 Adam Taranto
# This Source Code Form is subject to the terms of the Mozilla Public License,
# v. 2.0. If a copy of the MPL was not distributed with this file, You can
# obtain one at http://mozilla.org/MPL/2.0/.
from __future__ import print_function, division, absolute_import

from ._version import get_versions
__version__ = get_versions()['version']
FRISK_VERSION = __version__
del get_versions
