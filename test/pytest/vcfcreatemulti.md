% VCFCREATEMULTI(1) vcfcreatemulti (vcflib) | vcfcreatemulti (VCF transformation)
% Erik Garrison, Pjotr Prins and other vcflib contributors

# NAME

vcfcreatemulti - reduces complex alleles by pairwise alignment with BiWFA

# SYNOPSIS

**vcfcreatemulti**

# DESCRIPTION

**vcfcreatemulti** merges VCF records into one line by combining
  alleles.

Go through sorted VCF and if overlapping alleles are represented
across multiple records, merge them into a single record.  Currently
only for indels.

## Options

-h, --help

: shows help message and exits.

See more below.

# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# EXAMPLES


<!--

    >>> from rtest import run_stdout, head, cat, sh

-->

```

>>> head("vcfcreatemulti -h",25)
>
Usage: vcfcreatemulti [options] [file]
>
Go through sorted VCF and if overlapping alleles are represented
across multiple records, merge them into a single record.  Currently
only for indels.
>
Type: transformation
>

```



## Source code

[vcfcreatemulti.cpp](../../src/vcfcreatemulti.cpp)

## Regression tests

```python
# ./vcfcreatemulti ../samples/grch38#chr8_36353854-36453166.vcf > ../test/data/regression/vcfcreatemulti_2.vcf
>>> run_stdout("vcfcreatemulti ../samples/grch38#chr8_36353854-36453166.vcf", ext="vcf")
output in <a href="../data/regression/vcfcreatemulti_2.vcf">vcfcreatemulti_2.vcf</a>

>>> run_stdout("vcfcreatemulti ../samples/sample.vcf", ext="vcf")
output in <a href="../data/regression/vcfcreatemulti_3.vcf">vcfcreatemulti_3.vcf</a>

```

# LICENSE

Copyright 2022 (C) Erik Garrison, Pjotr Prins and vcflib contributors. MIT licensed.
