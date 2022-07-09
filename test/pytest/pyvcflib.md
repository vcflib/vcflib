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

So the one input record shows it has a ref of 'ACCCCCACCCCCACC' and six alt alleles ['ACCCCCACCCCCACC', 'ACC', 'AC', 'ACCCCCACCCCCAC', 'ACCCCCACC', 'ACA'].

This works fine!

## Correcting for genotypes

When two VCF records get combined we need to combine the genotypes. So, say we have as input two variants at the same position the genotypes need to be updated:

```python
>>> var1 = "1/0 0/0 0/1 1/1".split()
>>> var2 = "0/0 1/1 1/0 0/0".split() # ['0/0', '1/1', '1/0', '0/0']

>> merge_genotypes(var1,var2)
>> "1/0 2/2 2/1 1/1".split()

```

## Another example

See [realign.py](../tests/realign.py) for examples of using the FFI in the form of a python unit test. We used that to develop the vcfwave module.

## Additional bindings

It may be the case you want additional bindings that we have not included yet. See vcflib's [pythonffi.cpp](../../src/pythonffi.cpp) for the existing bindings and [Variant.h](../../src/Variant.h) for the classes and accessors.
