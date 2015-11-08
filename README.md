Frisk
=====

[![Join the chat at https://gitter.im/Adamtaranto/frisk](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/Adamtaranto/frisk?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

Screen genomic scaffolds for regions of unusual k-mer composition.

Install
-------

Frisk is a python package. On almost any modern operating system, the following
should work:

    pip install cython numpy scipy  # Install depenencies
    pip install frisk

Please note that as our `setup.py` imports both numpy, scipy and cython, these
packages must be installed **before** you call either `pip install frisk` or
`python setup.py ...`.


Documentation
-------------

Is currently being written, and will appear at
[frisk.rtfd.org](http://frisk.readthedocs.org/) as it is completed.

License
-------

Frisk is licensed under the GNU GPL version 3, or (at your option) any later
version. See ./LICENSE for the full license text.
