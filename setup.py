from setuptools import setup
import versioneer

versioneer.VCS = 'git'
versioneer.versionfile_source = 'frisk/_version.py'
versioneer.versionfile_build = 'frisk/_version.py'
versioneer.tag_prefix = 'v'  # tags are like v1.2.0
versioneer.parentdir_prefix = '.'
CMDCLASS = versioneer.get_cmdclass()

desc = """
frisk: Detection of sequence composition anomalies using multiple order kmers.
"""

install_requires = [
    "numpy",
    "pybedtools",
    "pandas",
    "matplotlib",
    "scipy",
    "seaborn",
    "versioneer",
]

test_requires = []

setup(
    name="frisk",
    py_modules=['frisk', ],
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    install_requires=install_requires,
    tests_require=test_requires,
    description=desc,
    author="Adam Taranto",
    author_email="adam.taranto@anu.edu.au",
    url="https://github.com/Adamtaranto/frisk",
    keywords=["kmer"],
    classifiers=[
        "Programming Language :: Python",
        "Programming Language :: Python :: 2",
        "Development Status :: 3 - Alpha",
        "Environment :: Console",
        "Intended Audience :: Developers",
        "Operating System :: OS Independent",
        "Topic :: Software Development :: Libraries :: Python Modules",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: GNU General Public License v3 or later " +
            "(GPLv3+)",
    ],
    test_suite="test",
)