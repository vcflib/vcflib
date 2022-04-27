% -*- coding: utf-8 -*-

# VCFLib Python FFI

We are building up a Python FFI for vcflib using the brilliant Python pybind11 module. Mostly for (our) testing purposes, so we are not aiming for complete coverage. See below for adding new bindings.

## Setting it up and example

First import the module. It may require setting the `PYTHONPATH` to the shared library `pyvcflib.cpython-39-x86_64-linux-gnu.so`.

    env PYTHONPATH=./build python3 -c 'import pyvcflib'

in a GNU Guix shell you may prepend `LD_LIBRARY_PATH` to find GLIBC etc.

    LD_LIBRARY_PATH=$LIBRARY_PATH

Now you should be able to use the `pyvcflib` module. Let's try with a VCF [samples/10158243.vcf](../../samples/10158243.vcf) that has only one record:

```python
>>> from pyvcflib import *

>>> vcf = VariantCallFile()
>>> vcf.openFile("../samples/10158243.vcf")
True

# ...    list(rec.name,rec.pos,rec.ref,rec.alt)

>>> rec = Variant(vcf)
>>> while (vcf.getNextVariant(rec)):
...    [rec.name,rec.pos,rec.ref]
['grch38#chr4', 10158243, 'ACCCCCACCCCCACC']

>>> rec.alt[0]
'ACC'

>>> rec.alleles
['ACCCCCACCCCCACC', 'ACC', 'AC', 'ACCCCCACCCCCAC', 'ACCCCCACC', 'ACA']

```

This works fine! It may be the case you want additional bindings that we have not included yet. See vcflib's [pythonffi.cpp](../../src/pythonffi.cpp) for the existing bindings and [Variant.h](../../src/Variant.h) for the classes and accessors.
