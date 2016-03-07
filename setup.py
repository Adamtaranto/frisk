#!/usr/bin/env python
# Copyright (C) 2015 Adam Taranto <adam.p.taranto@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from setuptools import setup
import versioneer

desc = """
frisk: Detection of sequence composition anomalies using multiple order kmers.
"""

install_requires = [
    "numpy>=1.10",
    "scipy>=0.17",
    "cython>=0.23",
    "pandas>=0.17",
    "matplotlib>=1.5.0",
    "seaborn>=0.8",
    "hmmlearn>=0.2.0",
    "pysam>=0.9.0",
    "pybedtools>=0.7.6",
    "scikit-learn>=0.17.1",
]

test_requires = [
    'nose',
]

pypi_classifiers = [
    "Programming Language :: Python",
    "Programming Language :: Python :: 2",
    "Development Status :: 3 - Alpha",
    "Environment :: Console",
    "Intended Audience :: Developers",
    "Operating System :: OS Independent",
    "Topic :: Software Development :: Libraries :: Python Modules",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
    "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
]

setup(
    name="frisk",
    packages=['frisk', ],
    version=versioneer.get_version(),
    install_requires=install_requires,
    tests_require=test_requires,
    description=desc,
    author="Adam Taranto",
    author_email="adam.taranto@anu.edu.au",
    url="https://github.com/Adamtaranto/frisk",
    keywords=["kmer"],
    classifiers=pypi_classifiers,
    entry_points={
        'console_scripts': [
            'frisk=frisk:main',
        ],
    },
)
