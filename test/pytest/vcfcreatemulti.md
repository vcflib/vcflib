% VCFCREATEMULTI(1) vcfcreatemulti (vcflib) | vcfcreatemulti (VCF transformation)
% Erik Garrison, Pjotr Prins and other vcflib contributors

# NAME

vcfcreatemulti - collates single ALT allele records into multi-allele records while tracking genotypes

# SYNOPSIS

**vcfcreatemulti**

# DESCRIPTION

**vcfcreatemulti** merges VCF records into one line by combining ALT alleles into a single VCF record. This tools is a great companion to [vcfwave](./vcfwave.md).

In 2022 **vcfcreatemulti** has been upgraded to track INFO records and genotypes (samples) so they are updated in the output.

See below EXAMPLES and [vcfwave](./vcfwave.md) for more information.

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
options:
     -n, --nextgen           next gen mode.
>
Type: transformation
>

```

vcfcreatemulti can combine overlapping alleles onto one record (VCF line), but it does not correct the INFO fields and sample (genotypes). For example:

```python

>>> sh("cat ../samples/10158243-after-vcfwave.vcf|grep -v ^\#")
grch38#chr4     10158244        >3655>3662_1    CCCCCACCCCCAC   C       60      .       AC=1;AF=0.011236;AN=89;AT=>3655>3656>3657>3660>3662;NS=45;LV=0;ORIGIN=grch38#chr4:10158243;LEN=12;INV=0;TYPE=del        GT      0|0     0|0     0|0     0|0     1|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0
grch38#chr4     10158244        >3655>3662_2    CCCCCACCCCCACC  C       60      .       AC=3;AF=0.033708;AN=89;AT=>3655>3656>3660>3662;NS=45;LV=0;ORIGIN=grch38#chr4:10158243;LEN=13;INV=0;TYPE=del     GT      0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     1|0     0|1     0|0     0|0     0|0     0|0     0|0     0|0     0|1     0|0     0
grch38#chr4     10158245        >3655>3662_3    CCCCACCCCCACC   C       60      .       AC=64;AF=0.719101;AN=89;AT=>3655>3656>3657>3658>3659>3660>3662;NS=45;LV=0;ORIGIN=grch38#chr4:10158243;LEN=12;INV=0;TYPE=del     GT      0|0     1|1     1|1     1|0     0|1     0|0     0|1     0|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     0|0     1|1     1|1     1|1     1|0     1|0     1|0     1|0     1|1     1|1     1|0     1|1     1|1     0|0     1|0     1|1     0|1     1|1     1|1     0|1     1|0     1|1     1|1     0|1     1|1     1|1     1|0     1|0     1|1     0
grch38#chr4     10158251        >3655>3662_4    CCCCACC C       60      .       AC=3;AF=0.033708;AN=89;AT=>3655>3656>3657>3658>3660>3662;NS=45;LV=0;ORIGIN=grch38#chr4:10158243;LEN=6;INV=0;TYPE=del    GT      0|0     0|0     0|0     0|0     0|0     0|1     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     1|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|1     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0
grch38#chr4     10158256        >3655>3662_5    CC      C       60      .       AC=2;AF=0.022472;AN=89;AT=>3655>3660>3662;NS=45;LV=0;ORIGIN=grch38#chr4:10158243;LEN=1;INV=0;TYPE=del   GT      0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|1     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     1|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0
grch38#chr4     10158257        >3655>3662_6    C       A       60      .       AC=1;AF=0.011236;AN=89;AT=>3655>3656>3657>3660>3662;NS=45;LV=0;ORIGIN=grch38#chr4:10158243;LEN=1;INV=0;TYPE=snp GT      0|0     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     0

```

gets converted into the following:

```python

>>> sh("../build/vcfcreatemulti -n ../samples/10158243-after-vcfwave.vcf|grep -v ^\#")
grch38#chr4     10158244        >3655>3662_1    CCCCCACCCCCACC  CC,C,CC,CCCCCACC,CCCCCACCCCCAC,CCCCCACCCCCACA   60      .       AC=1,3,64,3,2,1;AF=0.011236,0.033708,0.719101,0.033708,0.022472,0.011236;AN=89,89,89,89,89,89;AT=>3655>3656>3657>3660>3662,>3655>3656>3660>3662,>3655>3656>3657>3658>3659>3660>3662,>3655>3656>3657>3658>3660>3662,>3655>3660>3662,>3655>3656>3657>3660>3662;NS=45;LV=0;ORIGIN=grch38#chr4:10158243;LEN=12;INV=0,0,0,0,0,0;TYPE=del,del,del,del,del,snp;combined=10158244-10158257      GT      0|0     3|3     3|3     3|0     1|3     0|4     0|3     0|3     3|3     3|3     3|3     3|3     3|3     3|3     3|3     4|5     3|3     3|3     3|3     3|0     3|0     3|0     3|0     3|3     3|3     3|4     3|3     3|3     5|0     3|0     3|3     0|3     3|3     3|3     2|3     3|2     3|3     3|3     0|3     3|3     3|3     3|0     3|2     3|3     0

```

## Source code

[vcfcreatemulti.cpp](../../src/vcfcreatemulti.cpp)

## Regression tests

These tests mostly check for any major regressions between vcflib parser and outputter.
In the first example grch38#chr8 36377478,36394713,36409983 get combined

```python
# ./vcfcreatemulti -n ../samples/grch38#chr8_36353854-36453166.vcf > ../test/data/regression/vcfcreatemulti_2.vcf
>>> run_stdout("vcfcreatemulti ../samples/grch38#chr8_36353854-36453166.vcf", ext="vcf", uniq=2)
output in <a href="../data/regression/vcfcreatemulti_2.vcf">vcfcreatemulti_2.vcf</a>

>>> run_stdout("vcfcreatemulti -n ../samples/sample.vcf", ext="vcf", uniq=3)
output in <a href="../data/regression/vcfcreatemulti_3.vcf">vcfcreatemulti_3.vcf</a>

```

Check if the legacy version is still the same. Note it only retains the first genotype and has duplicate 'CC' alt alleles. INFO fields are not correct either.

```python
>>> sh("../build/vcfcreatemulti ../samples/10158243-after-vcfwave.vcf|grep -v ^\#")
grch38#chr4     10158244        >3655>3662_1    CCCCCACCCCCACC  CC,C,CC,CCCCCACC,CCCCCACCCCCAC,CCCCCACCCCCACA   60      .       AC=1;AF=0.011236;AN=89;AT=>3655>3656>3657>3660>3662;NS=45;LV=0;ORIGIN=grch38#chr4:10158243;LEN=12;INV=0;TYPE=del;combined=10158244-10158257     GT      0|0     0|0     0|0     0|0     1|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0

```

# LICENSE

Copyright 2022 (C) Erik Garrison, Pjotr Prins and vcflib contributors. MIT licensed.
