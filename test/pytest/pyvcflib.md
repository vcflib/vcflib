% -*- coding: utf-8 -*-

# VCFLib Python FFI

We are building up a Python FFI for vcflib. Mostly for (our) testing purposes.

## Setting it up

First import the module. It may require setting the `PYTHONPATH` to the shared library `pyvcflib.cpython-39-x86_64-linux-gnu.so`.

    env PYTHONPATH=./build python3 -c 'import pyvcflib'

in a GNU Guix shell you may prepend `LD_LIBRARY_PATH` to find GLIBC etc.

    LD_LIBRARY_PATH=$LIBRARY_PATH

Now you should be able to use the `pyvcflib` module

```python
>>> from pyvcflib import *

>>> vcf = VariantCallFile()
>>> vcf.openFile("../samples/10158243.vcf")
True

```
