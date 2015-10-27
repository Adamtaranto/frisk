#!/usr/bin/env python
from setuptools import setup
import versioneer

desc = """
frisk: Detection of sequence composition anomalies using multiple order kmers.
"""

install_requires = [
    "numpy",
    "scipy",
    "cython",
    "pandas",
    "matplotlib",
    "seaborn",
    "hmmlearn",
    "pybedtools",
    "scikit-learn",
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
