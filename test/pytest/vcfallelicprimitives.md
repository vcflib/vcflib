% VCFALLELICPRIMITIVES(1) vcfallelicprimitives | Convert VCF to TSV
% Erik Garrison, Pjotr Prins and other vcflib contributors

# NAME

vcfallelicprimitives - Converts stdin or given VCF file and reduces alleles.

# SYNOPSIS

**vcfallelicprimitives**

# DESCRIPTION

**vcfallelicprimitives** converts stdin or given VCF file to tab-delimited format,

If multiple allelic primitives (gaps or mismatches) are specified in a
single VCF record, split the record into multiple lines, but drop all
INFO fields. Does not handle genotypes (yet). MNPs are split into
multiple SNPs unless the -m flag is provided.


## Options

-h, --help

: shows help message and exits.

# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# EXAMPLES


<!--

    >>> from pytest.rtest import run_stdout, head, cat

-->

```

>>> head("vcfallelicprimitives -h",1)
usage: vcfallelicprimitives [options] [file]

```

vcfallelicprimitives picks complex regions and simplifies nested alignments


## Source code

[vcfallelicprimitives.cpp](../../src/vcfallelicprimitives.cpp)

## Regression tests

```python
>>> run_stdout("vcfallelicprimitives -m -L 1000 ../samples/grch38#chr8_36353854-36453166.vcf", ext="vcf")
output in <a href="../data/regression/vcfallelicprimitives_2.vcf">vcfallelicprimitives_2.vcf</a>

>>> run_stdout("vcfallelicprimitives -m -L 1000 ../samples/grch38#chr4_10083863-10181258.vcf", ext="vcf")
output in <a href="../data/regression/vcfallelicprimitives_3.vcf">vcfallelicprimitives_3.vcf</a>

```

# LICENSE

Copyright 2022 (C) Erik Garrison, Pjotr Prins and vcflib contributors. MIT licensed.
