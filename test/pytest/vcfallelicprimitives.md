% VCFALLELICPRIMITIVES(1) vcfallelicprimitives (vcflib) | vcfallelicprimitives (VCF transformation)
% Erik Garrison, Pjotr Prins and other vcflib contributors

# NAME

vcfallelicprimitives - Converts stdin or given VCF file and reduces alleles.

# SYNOPSIS

**vcfallelicprimitives**

# DESCRIPTION

**vcfallelicprimitives** converts stdin or given VCF file to tab-delimited format,

Realign reference and alternate alleles with WFA, parsing out the primitive alleles
into multiple VCF records. New records have IDs that reference the source record ID.
Genotypes are handled. Deletion alleles will result in haploid (missing allele) genotypes.

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

    >>> from pytest.rtest import run_stdout, head, cat, sh

-->

```

>>> head("vcfallelicprimitives -h",28)
>
usage: ./vcfallelicprimitives [options] [file]
>
Note that this tool is considered legacy and will emit a warning!  Use
[vcfwave](./vcfwave.md) instead.
>
Realign reference and alternate alleles with WFA or SW, parsing out
the primitive alleles into multiple VCF records. New records have IDs
that reference the source record ID.  Genotypes are handled. Deletion
alleles will result in haploid (missing allele) genotypes.
>
options:
    -a, --algorithm TYPE    Choose algorithm (default) Wave front or (obsolete)
                            Smith-Waterman [WF|SW] algorithm
    -m, --use-mnps          Retain MNPs as separate events (default: false).
    -t, --tag-parsed FLAG   Annotate decomposed records with the source record
                            position (default: ORIGIN).
    -L, --max-length LEN    Do not manipulate records in which either the ALT or
                            REF is longer than LEN (default: unlimited).
    -k, --keep-info         Maintain site and allele-level annotations when
                            decomposing.  Note that in many cases,
                            such as multisample VCFs, these won't be
                            valid post decomposition.  For biallelic
                            loci in single-sample VCFs, they should be
                            used with caution.
    -d, --debug             debug mode.
>
Type: transformation

```

vcfallelicprimitives picks complex regions and simplifies nested alignments. For example:

```python

>>> sh("grep 10158243 ../samples/10158243.vcf")
grch38#chr4     10158243        >3655>3662      ACCCCCACCCCCACC ACC,AC,ACCCCCACCCCCAC,ACCCCCACC,ACA     60      .       AC=64,3,2,3,1;AF=0.719101,0.0337079,0.0224719,0.0337079,0.011236;AN=89;AT=>3655>3656>3657>3658>3659>3660>3662,>3655>3656>3660>3662,>3655>3660>3662,>3655>3656>3657>3658>3660>3662,>3655>3656>3657>3660>3662,>3655>3656>3661>3662;NS=45;LV=0     GT      0|0     1|1     1|1     1|0     5|1     0|4     0|1     0|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     4|3     1|1     1|1     1|1     1|0     1|0     1|0     1|0     1|1     1|1     1|4     1|1     1|1     3|0     1|0     1|1     0|1     1|1     1|1     2|1     1|2     1|1     1|1     0|1     1|1     1|1     1|0     1|2     1|1     0

```

After aligning it reduces into two records with variant alleles

```
10158243:ACCCCCA/A
10158243:ACCCCCACCCC/A
10158243:ACCCCCACCCCCA/A
10158243:ACCCCCACCCCCAC/A

10158255:AC/A
10158255:ACC/A
```

and adjusts the genotypes accordingly splitting into two records using the original (but arguably OBSOLETE) SW algorithm:

```python

>>> sh("../build/vcfallelicprimitives -a SW -m -L 1000 ../samples/10158243.vcf|grep -v ^\#")
grch38#chr4     10158243        >3655>3662_1    ACCCCCACCCCCAC  ACCCCCAC,ACAC,AC,A      60      .       AC=3,1,64,3;AF=0.0337079,0.011236,0.719101,0.0337079;LEN=6,10,12,13;ORIGIN=grch38#chr4:10158243,grch38#chr4:10158243,grch38#chr4:10158243,grch38#chr4:10158243;TYPE=del,del,del,del     GT      0|0     3|3     3|3     3|0     2|3     0|1     0|3     0|3     3|3     3|3     3|3     3|3     3|3     3|3     3|3     1|0     3|3     3|3     3|3     3|0     3|0     3|0     3|0     3|3     3|3     3|1     3|3     3|3     0|0     3|0     3|3     0|3     3|3     3|3     4|3     3|4     3|3     3|3     0|3     3|3     3|3     3|0     3|4     3|3     0
grch38#chr4     10158255        >3655>3662_2    ACC     AC,A    60      .       AC=2,1;AF=0.0224719,0.011236;LEN=1,2;ORIGIN=grch38#chr4:10158243,grch38#chr4:10158243;TYPE=del,del      GT      0|0     .|.     .|.     .|0     2|.     0|0     0|.     0|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     0|1     .|.     .|.     .|.     .|0     .|0     .|0     .|0     .|.     .|.     .|0     .|.     .|.     1|0     .|0     .|.     0|.     .|.     .|.     .|.     .|.     .|.     .|.     0|.     .|.     .|.     .|0     .|.     .|.     0

```

With the default wavefront algorithm we get a different result

```
1$10158244:CCCCCACCCCCAC/C
2$10158244:CCCCCACCCCCACC/C
3$10158245:CCCCACCCCCACC/C
4$10158251:CCCCACC/C
5$10158256:CC/C
```

```python

>>> sh("../build/vcfallelicprimitives -m -L 1000 ../samples/10158243.vcf|grep -v ^\#")
grch38#chr4     10158244        >3655>3662_1    CCCCCACCCCCACC  CC,C    60      .       AC=1,3;AF=0.011236,0.0337079;LEN=12,13;ORIGIN=grch38#chr4:10158243,grch38#chr4:10158243;TYPE=del,del    GT      0|0     0|0     0|0     0|0     1|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     2|0     0|2     0|0     0|0     0|0     0|0     0|0     0|0     0|2     0|0     0
grch38#chr4     10158245        >3655>3662_2    CCCCACCCCCACC   C       60      .       AC=64;AF=0.719101;LEN=12;ORIGIN=grch38#chr4:10158243;TYPE=del   GT      0|0     1|1     1|1     1|0     .|1     0|0     0|1     0|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     0|0     1|1     1|1     1|1     1|0     1|0     1|0     1|0     1|1     1|1     1|0     1|1     1|1     0|0     1|0     1|1     0|1     1|1     1|1     .|1     1|.     1|1     1|1     0|1     1|1     1|1     1|0     1|.     1|1     0
grch38#chr4     10158251        >3655>3662_3    CCCCACC C       60      .       AC=3;AF=0.0337079;LEN=6;ORIGIN=grch38#chr4:10158243;TYPE=del    GT      0|0     .|.     .|.     .|0     .|.     0|1     0|.     0|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     1|0     .|.     .|.     .|.     .|0     .|0     .|0     .|0     .|.     .|.     .|1     .|.     .|.     0|0     .|0     .|.     0|.     .|.     .|.     .|.     .|.     .|.     .|.     0|.     .|.     .|.     .|0     .|.     .|.     0
grch38#chr4     10158256        >3655>3662_4    CC      C       60      .       AC=2;AF=0.0224719;LEN=1;ORIGIN=grch38#chr4:10158243;TYPE=del    GT      0|0     .|.     .|.     .|0     .|.     0|.     0|.     0|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|1     .|.     .|.     .|.     .|0     .|0     .|0     .|0     .|.     .|.     .|.     .|.     .|.     1|0     .|0     .|.     0|.     .|.     .|.     .|.     .|.     .|.     .|.     0|.     .|.     .|.     .|0     .|.     .|.     0

```


## Source code

[vcfallelicprimitives.cpp](../../src/vcfallelicprimitives.cpp)

## Regression tests

Note the wave front version has no problem with longer sequences:

```python
>>> run_stdout("vcfallelicprimitives -a SW -m -L 1000 ../samples/grch38#chr8_36353854-36453166.vcf", ext="vcf")
output in <a href="../data/regression/vcfallelicprimitives_5.vcf">vcfallelicprimitives_5.vcf</a>

>>> run_stdout("vcfallelicprimitives -a SW -m -L 1000 ../samples/grch38#chr4_10083863-10181258.vcf", ext="vcf")
output in <a href="../data/regression/vcfallelicprimitives_6.vcf">vcfallelicprimitives_6.vcf</a>

>>> run_stdout("vcfallelicprimitives -L 10000 -m ../samples/grch38#chr8_36353854-36453166.vcf", ext="vcf")
output in <a href="../data/regression/vcfallelicprimitives_7.vcf">vcfallelicprimitives_7.vcf</a>

>>> run_stdout("vcfallelicprimitives -m ../samples/grch38#chr4_10083863-10181258.vcf", ext="vcf")
output in <a href="../data/regression/vcfallelicprimitives_8.vcf">vcfallelicprimitives_8.vcf</a>

```

Another diff example where the first is SW and the second WFA2 showing:

```python
>>> sh("diff data/regression/vcfallelicprimitives_6.vcf data/regression/vcfallelicprimitives_8.vcf|tail -6")
1661c1633,1635
< grch38#chr4   10180508        >4593>4597_1    CTT     CTTT,CT,C       60      .       AC=7,47,1;AF=0.0786517,0.52809,0.011236;LEN=1,1,2;ORIGIN=grch38#chr4:10180508,grch38#chr4:10180508,grch38#chr4:10180508;TYPE=ins,del,del        GT      2|0     0|2     2|2     2|0     2|2     0|0     0|2     0|2     2|2     2|2     2|0     2|0     0|2     2|0     2|2     2|0     2|2     2|2     2|2     2|0     0|1     0|0     2|1     2|2     0|2     2|2     2|0     0|2     0|3     2|1     0|2     0|0     2|0     1|2     2|2     0|1     2|2     0|0     0|0     1|0     0|1     2|0     0|0     2|2     2
---
> grch38#chr4   10180508        >4593>4597_1    CTT     C       60      .       AC=1;AF=0.011236;LEN=2;ORIGIN=grch38#chr4:10180508;TYPE=del     GT      0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|1     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0
> grch38#chr4   10180509        >4593>4597_2    TT      T       60      .       AC=47;AF=0.52809;LEN=1;ORIGIN=grch38#chr4:10180508;TYPE=del     GT      1|0     0|1     1|1     1|0     1|1     0|0     0|1     0|1     1|1     1|1     1|0     1|0     0|1     1|0     1|1     1|0     1|1     1|1     1|1     1|0     0|0     0|0     1|0     1|1     0|1     1|1     1|0     0|1     0|.     1|0     0|1     0|0     1|0     0|1     1|1     0|0     1|1     0|0     0|0     0|0     0|0     1|0     0|0     1|1     1
> grch38#chr4   10180510        >4593>4597_3    T       TT      60      .       AC=7;AF=0.0786517;LEN=1;ORIGIN=grch38#chr4:10180508;TYPE=ins    GT      .|0     0|.     .|.     .|0     .|.     0|0     0|.     0|.     .|.     .|.     .|0     .|0     0|.     .|0     .|.     .|0     .|.     .|.     .|.     .|0     0|1     0|0     .|1     .|.     0|.     .|.     .|0     0|.     0|.     .|1     0|.     0|0     .|0     1|.     .|.     0|1     .|.     0|0     0|0     1|0     0|1     .|0     0|0     .|.     .

```

shows how WFA2 is doing a better job at taking things apart.

Even so, this record is wrong. From grch38#chr4_10083863-10181258.vcf

```
grch38#chr4 10084924  >33>38  AC  CC,CT 60  . AC=25,1;AF=0.280899,0.011236;AN=89;AT=>33>36>37>38,>33>34>37>
38,>33>34>35>38;NS=45;LV=0 GT  1|1 1|0 0|0 0|1 0|0 1|0 1|1 1|0 0|0 0|0 0|1 0|0 0|0 0|1 0|0 1|0 0|0 0|0 0|0
0|1 1|2 0|1 0|0 0|0 0|0 0|1 0|1 0|0 0|1 0|1 1|0 1|0 0|0 0|0 1|0 0|1 0|0 1|1 0|0 0|0 0|0 0|1 0|0 0|0 0
```

SW vs WFA2 makes it in

```
diff vcfallelicprimitives_6.vcf vcfallelicprimitives_8.vcf
< grch38#chr4   10084924        >33>38_1        AC      CC,CT   60      .       AC=25,1;AF=0.280899,0.011236;LEN=2,2;ORIGIN=grch38#chr4:10084924,grch38#chr4:10084924;TYPE=mnp,mnp      GT      1|1     1|0     0|0
     0|1     0|0     1|0     1|1     1|0     0|0     0|0     0|1     0|0     0|0     0|1     0|0     1|0
     0|0     0|0     0|0     0|1     1|2     0|1     0|0     0|0     0|0     0|1     0|1     0|0     0|1
     0|1     1|0     1|0     0|0     0|0     1|0     0|1     0|0     1|1     0|0     0|0     0|0     0|1
     0|0     0|0     0
---
> grch38#chr4   10084924        >33>38_1        AC      CC,CC   60      .       AC=1,25;AF=0.011236,0.280899;LEN=1,2;ORIGIN=grch38#chr4:10084924,grch38#chr4:10084924;TYPE=snp,mnp      GT      2|2     2|0     0|0
     0|2     0|0     2|0     2|2     2|0     0|0     0|0     0|2     0|0     0|0     0|2     0|0     2|0
     0|0     0|0     0|0     0|2     2|1     0|2     0|0     0|0     0|0     0|2     0|2     0|0     0|2
     0|2     2|0     2|0     0|0     0|0     2|0     0|2     0|0     2|2     0|0     0|0     0|0     0|2
     0|0     0|0     0
```

Where does the `CC` come from? and biWFA does

```
grep 10084924 *
vcfallelicprimitives_6.vcf:grch38#chr4  10084924        >33>38_1        AC      CC,CT   60      .       AC=25,1;AF=0.280899,0.011236;LEN=2,2;ORIGIN=grch38#chr4:10084924,grch38#chr4:10084924;TYPE=mnp,mnp GT      1|11|0     0|0     0|1     0|0     1|0     1|1     1|0     0|0     0|0     0|1     0|0     0|0     0|1     0|01|0     0|0     0|0     0|0     0|1     1|2     0|1     0|0     0|0     0|0     0|1     0|1     0|0     0|10|1     1|0     1|0     0|0     0|0     1|0     0|1     0|0     1|1     0|0     0|0     0|0     0|1     0|00|0     0
vcfallelicprimitives_8.vcf:grch38#chr4  10084924        >33>38_1        AC      CC,CC   60      .       AC=1,25;AF=0.011236,0.280899;LEN=1,2;ORIGIN=grch38#chr4:10084924,grch38#chr4:10084924;TYPE=snp,mnp GT      2|22|0     0|0     0|2     0|0     2|0     2|2     2|0     0|0     0|0     0|2     0|0     0|0     0|2     0|02|0     0|0     0|0     0|0     0|2     2|1     0|2     0|0     0|0     0|0     0|2     0|2     0|0     0|20|2     2|0     2|0     0|0     0|0     2|0     0|2     0|0     2|2     0|0     0|0     0|0     0|2     0|00|0     0
vcfwave_5.vcf:grch38#chr4       10084924        >33>38_1        AC      CC,CC   60      .       AC=1,25;AF=0.011236,0.280899;INV=0,0;LEN=1,2;ORIGIN=grch38#chr4:10084924,grch38#chr4:10084924;TYPE=snp,mnp GT      2|22|0     0|0     0|2     0|0     2|0     2|2     2|0     0|0     0|0     0|2     0|0     0|0     0|2     0|02|0     0|0     0|0     0|0     0|2     2|1     0|2     0|0     0|0     0|0     0|2     0|2     0|0     0|20|2     2|0     2|0     0|0     0|0     2|0     0|2     0|0     2|2     0|0     0|0     0|0     0|2     0|00|0     0
```

where

* SW: vcfallelicprimitives_6.vcf
* WFA2: vcfallelicprimitives_8.vcf
* biWFA: vcfwave_5.vcf

# LICENSE

Copyright 2022 (C) Erik Garrison, Pjotr Prins and vcflib contributors. MIT licensed.
