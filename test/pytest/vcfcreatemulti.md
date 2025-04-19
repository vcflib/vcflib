% VCFCREATEMULTI(1) vcfcreatemulti (vcflib) | vcfcreatemulti (VCF transformation)
% Erik Garrison, Pjotr Prins and other vcflib contributors

# NAME

vcfcreatemulti - collates single ALT allele records into multi-allele records while tracking genotypes

# SYNOPSIS

**vcfcreatemulti**

# DESCRIPTION

**vcfcreatemulti** merges VCF records into one line by combining ALT alleles into a single VCF record. This tool is a great companion to [vcfwave](./vcfwave.md).

In 2022 **vcfcreatemulti** has been upgraded to track INFO records and genotypes (samples) so they are updated in the output.

Note that the tool is not perfect:
See below EXAMPLES, the caveat on 'too many variants' `MULTI=ALTPROBLEM`, and [vcfwave](./vcfwave.md) for more information.

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
Go through sorted VCF and when overlapping alleles are represented across multiple records, merge them into a single multi-ALT record. See the documentation for more information.
>
options:
>
    --quiet           no progress bar
    --legacy          legacy mode (old C++ implementation does not do genotypes)
>
Type: transformation
>

```

The original 'legacy' vcfcreatemulti can combine overlapping alleles onto one record (VCF line), but it does not correct the INFO fields and sample (genotypes). For example:

```python

>>> sh("cat ../samples/10158243-after-vcfwave.vcf|grep -v ^\#")
grch38#chr4     10158244        >3655>3662_1    CCCCCACCCCCAC   C       60      .       AC=1;AF=0.011236;AN=89;AT=>3655>3656>3657>3660>3662;NS=45;LV=0;ORIGIN=grch38#chr4:10158243;LEN=12;INV=0;TYPE=del        GT      0|0     0|0     0|0     0|0     1|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0
grch38#chr4     10158244        >3655>3662_2    CCCCCACCCCCACC  C       60      .       AC=3;AF=0.033708;AN=89;AT=>3655>3656>3660>3662;NS=45;LV=0;ORIGIN=grch38#chr4:10158243;LEN=13;INV=0;TYPE=del     GT      0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     1|0     0|1     0|0     0|0     0|0     0|0     0|0     0|0     0|1     0|0     0
grch38#chr4     10158245        >3655>3662_3    CCCCACCCCCACC   C       60      .       AC=64;AF=0.719101;AN=89;AT=>3655>3656>3657>3658>3659>3660>3662;NS=45;LV=0;ORIGIN=grch38#chr4:10158243;LEN=12;INV=0;TYPE=del     GT      0|0     1|1     1|1     1|0     0|1     0|0     0|1     0|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     0|0     1|1     1|1     1|1     1|0     1|0     1|0     1|0     1|1     1|1     1|0     1|1     1|1     0|0     1|0     1|1     0|1     1|1     1|1     0|1     1|0     1|1     1|1     0|1     1|1     1|1     1|0     1|0     1|1     0
grch38#chr4     10158251        >3655>3662_4    CCCCACC C       60      .       AC=3;AF=0.033708;AN=89;AT=>3655>3656>3657>3658>3660>3662;NS=45;LV=0;ORIGIN=grch38#chr4:10158243;LEN=6;INV=0;TYPE=del    GT      0|0     0|0     0|0     0|0     0|0     0|1     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     1|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|1     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0
grch38#chr4     10158256        >3655>3662_5    CC      C       60      .       AC=2;AF=0.022472;AN=89;AT=>3655>3660>3662;NS=45;LV=0;ORIGIN=grch38#chr4:10158243;LEN=1;INV=0;TYPE=del   GT      0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|1     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     1|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0
grch38#chr4     10158257        >3655>3662_6    C       A       60      .       AC=1;AF=0.011236;AN=89;AT=>3655>3656>3657>3660>3662;NS=45;LV=0;ORIGIN=grch38#chr4:10158243;LEN=1;INV=0;TYPE=snp GT      0|0     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     0

```

this now gets converted into the following:

```python

>>> sh("vcfcreatemulti ../samples/10158243-after-vcfwave.vcf|grep -v ^\#")
grch38#chr4     10158244        >3655>3662_1    CCCCCACCCCCACC  CC,C,CC,CCCCCACC,CCCCCACCCCCAC,CCCCCACCCCCACA   60      .       AC=1,3,64,3,2,1;AF=0.011236,0.033708,0.719101,0.033708,0.022472,0.011236;AN=89,89,89,89,89,89;AT=>3655>3656>3657>3660>3662,>3655>3656>3660>3662,>3655>3656>3657>3658>3659>3660>3662,>3655>3656>3657>3658>3660>3662,>3655>3660>3662,>3655>3656>3657>3660>3662;NS=45;LV=0;ORIGIN=grch38#chr4:10158243;LEN=12;INV=0,0,0,0,0,0;TYPE=del,del,del,del,del,snp;combined=10158244-10158257      GT      0|0     3|3     3|3     3|0     1|3     0|4     0|3     0|3     3|3     3|3     3|3     3|3     3|3     3|3     3|3     4|5     3|3     3|3     3|3     3|0     3|0     3|0     3|0     3|3     3|3     3|4     3|3     3|3     5|0     3|0     3|3     0|3     3|3     3|3     2|3     3|2     3|3     3|3     0|3     3|3     3|3     3|0     3|2     3|3     0

```

### Too many variants

*Or the MULTI=ALTPROBLEM.*

That looks proper. There is one caveat or blatant problem, however. If a variant sequence is long (a large 'bubble') and with the other alleles more (small) INDELs are scored than there are samples then the genotypes represent only the last match. Resulting in something ugly:

```
...,del,del,snp,ins,snp,snp,snp,snp,snp,snp,snp,snp,snp,snp,snp,snp,snp,snp,snp,snp,snp,snp,snp,snp,snp,snp,snp,snp,snp,snp,snp,snp,snp,snp,snp,snp,snp,snp,snp,snp,snp,snp,snp,snp,snp,snp,snp,snp,snp,snp,snp,s np,snp,snp,snp,snp,snp,snp,snp,snp,snp,snp,snp,snp,snp,snp,snp,snp,snp,snp,snp,snp,snp,snp,snp,snp,snp,snp,snp,s np,snp,snp,snp,snp,snp,snp,snp,snp,snp,snp,snp,snp,snp,snp,snp,snp,snp,snp;combined=36390210-36409660 GT

509|49 8 500|500 20|251  500|238 238|498 653|387 102|1 500|498 9|509 498|69  500|297 498|725 498|660 500|472 204|500 50 0|846 654|653 500|500 500|500 18|18 430|498 214|500 499|299 67|500  18|386  47|154  508|47  500|385 42|47 579|47 47|18 47|47 219|500 18|47 53|213  500|18  500|18  500|500 47|846  47|47 500|47  500|47  839|500 498|47  500
```

this is a simple artefact resulting from the fact that complex structures do not map (easily) on the simple VCF layout. Another problem is that multiple SNPs don't get incorporated in the ALTs - the algorithm uses the reference to build up the longer ALTs for every single SNP from the start.

In this example you see that the ALT for SNP2 does not contain SNP1 even though they may be in the same individual/sample:

```
                       sample
REF      ACTGACTGACTG
ALT-SNP1 ACTGC          1/0
ALT-SNP2 ACTGACTA       1/0
             ^
```

In words: the result is incorrect.

At this point, for analysis, there is little else to do but go to the original data (pangenome or VCF) and compare the results.
What `vcfcreatemulti` helps to do is point out that there is a complex region here with ample variation and the resulting layout is a problem (too many ALTs as in 'too many cooks'!).

To help vcflib show's a `WARNING: Too many ALT alleles to fit in sample(s)' and we add an INFO tag "MULTI=ALTPROBLEM". Searching for these will give an idea of this issue. E.g.

```
grep MULTI= ./test/tmp/vcfcreatemulti_2.vcf -c
```

Finds 3 marked records. One of them is derived from the combination:

```
grch38#chr8 36377478  >601>606  GTTTCTTGAAAAACCAAATGT GTTTCTTGAAAAACCAAAGGT,G 60  . AC=20,1;AF=0.224719,0.011236
;AN=89;AT=>601>602>603>605>606,>601>602>604>605>606,>601>606;NS=45;LV=0 GT  0|0 0|0 0|0 0|0 0|0 0|0 0|0 0|0 0|0
0|0 0|0 0|0 0|0 0|0 0|0 0|0 0|0 0|0 0|0 1|1 0|0 0|0 0|0 0|0 1|0 1|0 0|2 0|0 0|1 0|1 1|1 1|1 0|0 1|1 0|0 0|1 0|1
0|0 1|0 1|1 0|1 0|1 0|0 0|1 0
grch38#chr8 36377496  >602>605  T G 60  . AC=20;AF=0.227273;AN=88;AT=>602>603>605,>602>604>605;NS=45;LV=1;PS=>60
1>606 GT  0|0 0|0 0|0 0|0 0|0 0|0 0|0 0|0 0|0 0|0 0|0 0|0 0|0 0|0 0|0 0|0 0|0 0|0 0|0 1|1 0|0 0|0 0|0 0|0 1|0 1|
0 0|. 0|0 0|1 0|1 1|1 1|1 0|0 1|1 0|0 0|1 0|1 0|0 1|0 1|1 0|1 0|1 0|0 0|1 0
```

resulting in

```
grch38#chr8 36377478  >601>606  GTTTCTTGAAAAACCAAATGT GTTTCTTGAAAAACCAAAGGT,G,GTTTCTTGAAAAACCAAAGGT 60  . AC=20,
1,20;AF=0.224719,0.011236,0.227273;AN=89,89,88;AT=>601>602>603>605>606,>601>602>604>605>606,>601>602>603>605>606
,>601>606,>602>603>605,>602>604>605;NS=45;LV=0;MULTI=ALTPROBLEM;combined=36377478-36377496  GT  0|0 0|0 0|0 0|0
0|0 0|0 0|0 0|0 0|0 0|0 0|0 0|0 0|0 0|0 0|0 0|0 0|0 0|0 0|0 3|3 0|0 0|0 0|0 0|0 3|0 3|0 0|2 0|0 0|3 0|3 3|3 3|3
0|0 3|3 0|0 0|3 0|3 0|0 3|0 3|3 0|3 0|3 0|0 0|3 0
```

This is a combination of:

```
grch38#chr8 36377478  >601>606  GTTTCTTGAAAAACCAAATGT GTTTCTTGAAAAACCAAAGGT,G 60  . AC=20,1;AF=0.224719,0.011236
;AN=89;AT=>601>602>603>605>606,>601>602>604>605>606,>601>606;NS=45;LV=0 GT  0|0 0|0 0|0 0|0 0|0 0|0 0|0 0|0 0|0
0|0 0|0 0|0 0|0 0|0 0|0 0|0 0|0 0|0 0|0 1|1 0|0 0|0 0|0 0|0 1|0 1|0 0|2 0|0 0|1 0|1 1|1 1|1 0|0 1|1 0|0 0|1 0|1
0|0 1|0 1|1 0|1 0|1 0|0 0|1 0
grch38#chr8 36377496  >602>605  T G 60  . AC=20;AF=0.227273;AN=88;AT=>602>603>605,>602>604>605;NS=45;LV=1;PS=>60
1>606 GT  0|0 0|0 0|0 0|0 0|0 0|0 0|0 0|0 0|0 0|0 0|0 0|0 0|0 0|0 0|0 0|0 0|0 0|0 0|0 1|1 0|0 0|0 0|0 0|0 1|0 1|
0 0|. 0|0 0|1 0|1 1|1 1|1 0|0 1|1 0|0 0|1 0|1 0|0 1|0 1|1 0|1 0|1 0|0 0|1 0
```

Where the ALTs end up being a duplication and there is some overlap in the genotype calling.

One future solution might be to have vcfcreatemulti ignore SNPs, or only take the first one, but that somewhat would do away with pointing out complex arrangements. Another solution might be to edit the ALTs and merge ALT-SNP1 into ALT-SNP2 so we get `ACTGCCTA`.
Contributions and ideas are welcome!

Having a think about this: the safest approach is to backtrack on a conflict and leave it alone. So, when a variant comes up that conflicts with the combined record (so far) we should drop merging that variant and leave it alone. This will typically happen with a long ALT that overlaps many SNPs. We could come up with all types of solutions, but the point of this algorithm is to 'fix' the obvious cases. At this point we continue and show the MULTI=ALTPROBLEM info field. It is not satisfactory and it is slow too. We can have a stab at the backtrack in the future.

## Source code

[vcfcreatemulti.cpp](../../src/vcfcreatemulti.cpp)

## Regression tests

These tests mostly check for any major regressions between vcflib parser and outputter.
In the first example grch38#chr8 36377478,36394713,36409983 get combined

```python
# ./vcfcreatemulti ../samples/grch38#chr8_36353854-36453166.vcf > ../test/data/regression/vcfcreatemulti_2.vcf
>>> run_stdout("vcfcreatemulti ../samples/grch38#chr8_36353854-36453166-bcftools-normalised.vcf", ext="vcf", uniq=2)
output in <a href="../data/regression/vcfcreatemulti_2.vcf">vcfcreatemulti_2.vcf</a>

>>> run_stdout("vcfcreatemulti ../samples/sample.vcf", ext="vcf", uniq=3)
output in <a href="../data/regression/vcfcreatemulti_3.vcf">vcfcreatemulti_3.vcf</a>

```

Check if the legacy version is still the same. Note it only retains the first genotype and has duplicate 'CC' alt alleles. INFO fields are not correct either.

```python
>>> sh("vcfcreatemulti --legacy ../samples/10158243-after-vcfwave.vcf|grep -v ^\#")
grch38#chr4     10158244        >3655>3662_1    CCCCCACCCCCACC  CC,C,CC,CCCCCACC,CCCCCACCCCCAC,CCCCCACCCCCACA   60      .       AC=1;AF=0.011236;AN=89;AT=>3655>3656>3657>3660>3662;NS=45;LV=0;ORIGIN=grch38#chr4:10158243;LEN=12;INV=0;TYPE=del;combined=10158244-10158257     GT      0|0     0|0     0|0     0|0     1|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0

```

## Trouble shooting

If you get an error like

```
thread 502 panic: attempt to unwrap error: MultiAltNotSupported
```

It means the input file already contains multi-allele VCF records. To split these you can run a command such as `bcftools norm -m-` to normalise the VCF records and split out multiple ALT alleles into separate VCF records.
Finally use **vcfcreatemulti** to create multi-allele VCF records again.

### Warning: Too many ALT alleles to fit in sample(s)

See 'caveat' section [above](#Too-many-variants).

### Warning: This code only supports one ALT allele per record: bailing out --- try normalising the data with `bcftools norm -m-`

Your VCF already contains multi-allele entries - bring them back to one single ALT per record/line.

# LICENSE

Copyright 2022-2024 (C) Erik Garrison, Pjotr Prins and vcflib contributors. MIT licensed.
