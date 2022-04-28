% VCFWAVE(1) vcfwave (vcflib) | vcfwave (VCF transformation)
% Erik Garrison, Pjotr Prins and other vcflib contributors

# NAME

vcfwave - reduces complex alleles by pairwise alignment with BiWFA

# SYNOPSIS

**vcfwave**

# DESCRIPTION

**vcfwave** reduces complex alleles into simpler primitive representation using pairwise
alignment with BiWFA.

Realigns reference and alternate alleles with WFA, parsing out the primitive alleles
into multiple VCF records. New records have IDs that reference the source record ID.
Genotypes are handled. Deletion alleles will result in haploid (missing allele) genotypes
overlapping the deleted region.

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

>>> head("vcfwave -h",22)
usage: vcfwave [options] [file]
>
Realign reference and alternate alleles with WFA, parsing out the primitive alleles
into multiple VCF records. New records have IDs that reference the source record ID.
Genotypes are handled. Deletions generate haploid/missing genotypes at overlapping sites.
>
options:
    -p, --wf-params PARAMS  use the given BiWFA params (default: 0,19,39,3,81,1)
                            format=match,mismatch,gap1-open,gap1-ext,gap2-open,gap2-ext
    -m, --use-mnps          Retain MNPs as separate events (default: false).
    -t, --tag-parsed FLAG   Annotate decomposed records with the source record position
                            (default: ORIGIN).
    -L, --max-length LEN    Do not manipulate records in which either the ALT or
                            REF is longer than LEN (default: unlimited).
    -K, --inv-kmer K        Length of k-mer to use for inversion detection sketching (default: 17).
    -I, --inv-min LEN       Minimum allele length to consider for inverted alignment (default: 64).
    -k, --keep-info         Maintain site and allele-level annotations when decomposing.
                            Note that in many cases, such as multisample VCFs, these won't
                            be valid post-decomposition.  For biallelic loci in single-sample
                            VCFs, they should be usable with caution.
    -j, --threads N         use this many threads for variant decomposition
    -d, --debug             debug mode.
>
Type: transformation

```

vcfwave picks complex regions and simplifies nested alignments. For example:

```python

>>> sh("grep 10158243 ../samples/10158243.vcf")
grch38#chr4     10158243        >3655>3662      ACCCCCACCCCCACC ACC,AC,ACCCCCACCCCCAC,ACCCCCACC,ACA     60      .       AC=64,3,2,3,1;AF=0.719101,0.0337079,0.0224719,0.0337079,0.011236;AN=89;AT=>3655>3656>3657>3658>3659>3660>3662,>3655>3656>3660>3662,>3655>3660>3662,>3655>3656>3657>3658>3660>3662,>3655>3656>3657>3660>3662,>3655>3656>3661>3662;NS=45;LV=0     GT      0|0     1|1     1|1     1|0     5|1     0|4     0|1     0|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     4|3     1|1     1|1     1|1     1|0     1|0     1|0     1|0     1|1     1|1     1|4     1|1     1|1     3|0     1|0     1|1     0|1     1|1     1|1     2|1     1|2     1|1     1|1     0|1     1|1     1|1     1|0     1|2     1|1     0

```

This aligns and adjusts the genotypes accordingly splitting into multiple records, one for each unique allele found in the alignments:

```python

>> sh("../build/vcfwave -m -L 1000 ../samples/10158243.vcf|grep -v ^\#")
grch38#chr4     10158244        >3655>3662_1    CCCCCACCCCCACC  CC,C    60      .       AC=1,3;AF=0.011236,0.0337079;INV=0,0;LEN=12,13;ORIGIN=grch38#chr4:10158243,grch38#chr4:10158243;TYPE=del,del    GT      0|0   0|0      0|0     0|0     1|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0   0|0      0|0     0|0     0|0     0|0     0|0     2|0     0|2     0|0     0|0     0|0     0|0     0|0     0|0     0|2     0|0     0
grch38#chr4     10158245        >3655>3662_2    CCCCACCCCCACC   C       60      .       AC=64;AF=0.719101;INV=0;LEN=12;ORIGIN=grch38#chr4:10158243;TYPE=del     GT      0|0     1|1     1|1     1|0     .|1     0|0   0|1      0|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     0|0     1|1     1|1     1|1     1|0     1|0     1|0     1|0     1|1     1|1     1|0     1|1     1|1     0|0     1|0     1|1     0|1     1|1   1|1      .|1     1|.     1|1     1|1     0|1     1|1     1|1     1|0     1|.     1|1     0
grch38#chr4     10158251        >3655>3662_3    CCCCACC C       60      .       AC=3;AF=0.0337079;INV=0;LEN=6;ORIGIN=grch38#chr4:10158243;TYPE=del      GT      0|0     .|.     .|.     .|0     .|.     0|1     0|.   0|.      .|.     .|.     .|.     .|.     .|.     .|.     .|.     1|0     .|.     .|.     .|.     .|0     .|0     .|0     .|0     .|.     .|.     .|1     .|.     .|.     0|0     .|0     .|.     0|.     .|.     .|.   .|.      .|.     .|.     .|.     0|.     .|.     .|.     .|0     .|.     .|.     0
grch38#chr4     10158256        >3655>3662_4    CC      C       60      .       AC=2;AF=0.0224719;INV=0;LEN=1;ORIGIN=grch38#chr4:10158243;TYPE=del      GT      0|0     .|.     .|.     .|0     .|.     0|.     0|.   0|.      .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|1     .|.     .|.     .|.     .|0     .|0     .|0     .|0     .|.     .|.     .|.     .|.     .|.     1|0     .|0     .|.     0|.     .|.     .|.   .|.      .|.     .|.     .|.     0|.     .|.     .|.     .|0     .|.     .|.     0

```

These are currently not left aligned. Stay tuned...

We can also handle inversions.
This test case includes one that was introduced by building a variation graph with an inversion and then decomposing it into a VCF with `vg deconstruct` and finally "popping" the inversion variant with [`vcfbub`](https://github.com/pangenome/vcfbub).

```python

>> sh("../build/vcfwave ../test/data/regression/z.vcf|grep -v ^\#")
a       293     >1>9_1  A       T       60      .       AC=1;AF=1;INV=1;LEN=1;ORIGIN=a:281;TYPE=snp     GT      1
a       310     >1>9_2  T       C       60      .       AC=1;AF=1;INV=1;LEN=1;ORIGIN=a:281;TYPE=snp     GT      1
a       329     >1>9_3  T       A       60      .       AC=1;AF=1;INV=1;LEN=1;ORIGIN=a:281;TYPE=snp     GT      1

```

## Source code

[vcfwave.cpp](../../src/vcfwave.cpp)

## Regression tests

The bidirectional wavefront (BiWFA) version has no problem with longer sequences:

```python
>> run_stdout("vcfwave -m -L 1000 ../samples/grch38#chr8_36353854-36453166.vcf", ext="vcf")
output in <a href="../data/regression/vcfallelicprimitives_5.vcf">vcfallelicprimitives_5.vcf</a>

>> run_stdout("vcfwave -m -L 1000 ../samples/grch38#chr4_10083863-10181258.vcf", ext="vcf")
output in <a href="../data/regression/vcfallelicprimitives_6.vcf">vcfallelicprimitives_6.vcf</a>

```

# LICENSE

Copyright 2022 (C) Erik Garrison, Pjotr Prins and vcflib contributors. MIT licensed.
