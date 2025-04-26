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

Note that this tool is considered legacy and will emit a warning!  Use
[vcfwave](./vcfwave.md) instead.

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

>>> head("vcfallelicprimitives -h",29)
>
usage: ./vcfallelicprimitives [options] [file]
>
WARNING: this tool is considered legacy and is only retained for older
workflows.  It will emit a warning!  Even though it can use the WFA
you should use [vcfwave](./vcfwave.md) instead.
>
Realign reference and alternate alleles with SW or WF, parsing out
the primitive alleles into multiple VCF records. New records have IDs
that reference the source record ID.  Genotypes are handled. Deletion
alleles will result in haploid (missing allele) genotypes.
>
options:
    -a, --algorithm TYPE    Choose algorithm SW (Smith-Waterman) or WF wavefront
                            (default: WF)
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

>>> sh("vcfallelicprimitives -a SW -m -L 1000 ../samples/10158243.vcf|grep -v ^\#")
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

>>> sh("vcfallelicprimitives -m -L 1000 ../samples/10158243.vcf|grep -v ^\#")
grch38#chr4     10158244        >3655>3662_1    CCCCCACCCCCACC  CC,C    60      .       AC=1,3;AF=0.011236,0.0337079;LEN=12,13;ORIGIN=grch38#chr4:10158243,grch38#chr4:10158243;TYPE=del,del    GT      0|0     0|0     0|0     0|0     1|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     2|0     0|2     0|0     0|0     0|0     0|0     0|0     0|0     0|2     0|0     0
grch38#chr4     10158245        >3655>3662_2    CCCCACCCCCACC   C       60      .       AC=64;AF=0.719101;LEN=12;ORIGIN=grch38#chr4:10158243;TYPE=del   GT      0|0     1|1     1|1     1|0     .|1     0|0     0|1     0|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     0|0     1|1     1|1     1|1     1|0     1|0     1|0     1|0     1|1     1|1     1|0     1|1     1|1     0|0     1|0     1|1     0|1     1|1     1|1     .|1     1|.     1|1     1|1     0|1     1|1     1|1     1|0     1|.     1|1     0
grch38#chr4     10158251        >3655>3662_3    CCCCACC C       60      .       AC=3;AF=0.0337079;LEN=6;ORIGIN=grch38#chr4:10158243;TYPE=del    GT      0|0     .|.     .|.     .|0     .|.     0|1     0|.     0|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     1|0     .|.     .|.     .|.     .|0     .|0     .|0     .|0     .|.     .|.     .|1     .|.     .|.     0|0     .|0     .|.     0|.     .|.     .|.     .|.     .|.     .|.     .|.     0|.     .|.     .|.     .|0     .|.     .|.     0
grch38#chr4     10158256        >3655>3662_4    CC      C       60      .       AC=2;AF=0.0224719;LEN=1;ORIGIN=grch38#chr4:10158243;TYPE=del    GT      0|0     .|.     .|.     .|0     .|.     0|.     0|.     0|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|1     .|.     .|.     .|.     .|0     .|0     .|0     .|0     .|.     .|.     .|.     .|.     .|.     1|0     .|0     .|.     0|.     .|.     .|.     .|.     .|.     .|.     .|.     0|.     .|.     .|.     .|0     .|.     .|.     0
grch38#chr4     10158257        >3655>3662_5    C       A       60      .       AC=1;AF=0.011236;LEN=1;ORIGIN=grch38#chr4:10158243;TYPE=snp     GT      0|0     .|.     .|.     .|0     .|.     0|.     0|.     0|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|.     .|0     .|0     .|0     .|0     .|.     .|.     .|.     .|.     .|.     .|0     .|0     .|.     0|.     .|.     .|.     .|.     .|.     .|.     .|.     0|.     .|.     .|.     .|0     .|.     .|.     0

```


## Source code

[vcfallelicprimitives.cpp](../../src/vcfallelicprimitives.cpp)

## Regression tests

Note the wave front version has no problem with longer sequences:

```python
# ./vcfallelicprimitives -a SW -m -L 1000 ../samples/grch38#chr8_36353854-36453166.vcf > ../test/data/regression/vcfallelicprimitives_5.vcf
>>> run_stdout("vcfallelicprimitives -a SW -m -L 1000 ../samples/grch38#chr8_36353854-36453166.vcf", ext="vcf")
output in <a href="../data/regression/vcfallelicprimitives_5.vcf">vcfallelicprimitives_5.vcf</a>

# ./vcfallelicprimitives -a SW -m -L 1000 ../samples/grch38#chr4_10083863-10181258.vcf > ../test/data/regression/vcfallelicprimitives_6.vcf
/regression/vcfallelicprimitives_6.vcf
>>> run_stdout("vcfallelicprimitives -a SW -m -L 1000 ../samples/grch38#chr4_10083863-10181258.vcf", ext="vcf")
output in <a href="../data/regression/vcfallelicprimitives_6.vcf">vcfallelicprimitives_6.vcf</a>

# ./vcfallelicprimitives -L 10000 -m ../samples/grch38#chr8_36353854-36453166.vcf > ../test/data/regression/vcfallelicprimitives_7.vcf
>>> run_stdout("vcfallelicprimitives -L 10000 -m ../samples/grch38#chr8_36353854-36453166.vcf", ext="vcf")
output in <a href="../data/regression/vcfallelicprimitives_7.vcf">vcfallelicprimitives_7.vcf</a>

# ./vcfallelicprimitives -m ../samples/grch38#chr4_10083863-10181258.vcf > ../test/data/regression/vcfallelicprimitives_8.vcf
>>> run_stdout("vcfallelicprimitives -m ../samples/grch38#chr4_10083863-10181258.vcf", ext="vcf")
output in <a href="../data/regression/vcfallelicprimitives_8.vcf">vcfallelicprimitives_8.vcf</a>

```

Another diff example where the first is SW and the second WFA2 showing:

```python
>>> sh("diff data/regression/vcfallelicprimitives_6.vcf data/regression/vcfallelicprimitives_8.vcf|tail -6")
1670c1680,1682
< grch38#chr4   10180508        >4593>4597_1    CTT     CTTT,CT,C       60      .       AC=7,47,1;AF=0.0786517,0.52809,0.011236;LEN=1,1,2;ORIGIN=grch38#chr4:10180508,grch38#chr4:10180508,grch38#chr4:10180508;TYPE=ins,del,del        GT      2|0     0|2     2|2     2|0     2|2     0|0     0|2     0|2     2|2     2|2     2|0     2|0     0|2     2|0     2|2     2|0     2|2     2|2     2|2     2|0     0|1     0|0     2|1     2|2     0|2     2|2     2|0     0|2     0|3     2|1     0|2     0|0     2|0     1|2     2|2     0|1     2|2     0|0     0|0     1|0     0|1     2|0     0|0     2|2     2
---
> grch38#chr4   10180508        >4593>4597_1    CTT     C       60      .       AC=1;AF=0.011236;LEN=2;ORIGIN=grch38#chr4:10180508;TYPE=del     GT      0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|1     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0
> grch38#chr4   10180509        >4593>4597_2    TT      T       60      .       AC=47;AF=0.52809;LEN=1;ORIGIN=grch38#chr4:10180508;TYPE=del     GT      1|0     0|1     1|1     1|0     1|1     0|0     0|1     0|1     1|1     1|1     1|0     1|0     0|1     1|0     1|1     1|0     1|1     1|1     1|1     1|0     0|0     0|0     1|0     1|1     0|1     1|1     1|0     0|1     0|.     1|0     0|1     0|0     1|0     0|1     1|1     0|0     1|1     0|0     0|0     0|0     0|0     1|0     0|0     1|1     1
> grch38#chr4   10180510        >4593>4597_3    T       TT      60      .       AC=7;AF=0.0786517;LEN=1;ORIGIN=grch38#chr4:10180508;TYPE=ins    GT      .|0     0|.     .|.     .|0     .|.     0|0     0|.     0|.     .|.     .|.     .|0     .|0     0|.     .|0     .|.     .|0     .|.     .|.     .|.     .|0     0|1     0|0     .|1     .|.     0|.     .|.     .|0     0|.     0|.     .|1     0|.     0|0     .|0     1|.     .|.     0|1     .|.     0|0     0|0     1|0     0|1     .|0     0|0     .|.     .

```

shows how WFA2 is doing a better job at taking things apart.

Even so, this record is wrong. From grch38#chr4_10083863-10181258.vcf


Differences between WFA and biWFA:

```
wdiff vcfallelicprimitives_8.vcf vcfwave_5.vcf
grch38#chr4     10134337        >2103>2106_1    [-TTTTG AGGCA-] {+TTTTGGTGTACTGCCT      AGGCAGTACACCAAAA+}
grch38#chr4     10134492        >2125>2211_3    [-TTGATGGAGAATCCCAATTGATGGAGAATCCCAATTGATGGAGAATCCCAATTGATGGAGAATCCCAATTGATGGAGAATCCCAATTGATGGAGAATCCCAATTGATGGAGAATCCCAATTG    GTGATGTAGCATCCCAAGTGATGTAGCATCCCAAGTGATGTAGCATCCCAAGTGATGTAGCATCCCAAGTGATGTAGCATCCCAAGTGATGTAGCATCCCAAGTG-]     {+TTGATGGAGAATCCCAATTGATGGAGAATCCCAATTGATGGAGAATCCCAATTGATGGAGAATCCCAATTGATGGAGAATCCCAATTGATGGAGAATCCCAATTGATGGAGAATCCCAATTGATGGAGAATCCCAATTGATGGAG     GTGATGTAGCATCCCAAGTGATGTAGCATCCCAAGTGATGTAGCATCCCAAGTGATGTAGCATCCCAAGTGATGTAGCATCCCAAGTGATGTAGCATCCCAAGTGATGTAGCATCCCAAGTGATGTAG+}
grch38#chr4     10134498
grch38#chr4     [-10134501      >2125>2211_5    AATCCCAATTGATGGAGAATCCCAATTGATGGAGAATCCCAATTGATGGAGAATCCCAATTGATGGAGAATCCCAATTGATGGAGAATCCCAATTGATGGAGAATCCCAATTGATGGAGAATCCCAATT       CATCCCAATTGATGGAGCATCCCAATTGATGGAGCATCCCAATTGATGGAGCATCCCAATTGATGGAGCATCCCAATTGATGGAGCATCCCAATTGATGGAGCATCCCAATT-]      {+10134500      >2125>2211_6    GAATCCCAATTGATGGAG      G+}     60      .       [-AC=1;AF=0.011236;INV=0,0;LEN=17;ORIGIN=grch38#chr4:10134484;TYPE=complex-]    {+AC=1;AF=0.011236;INV=0;LEN=17;ORIGIN=grch38#chr4:10134484;TYPE=del+}
  GT      0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|.     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0
     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0
     0|0     0|1     0|0     0|0     0|0     0
grch38#chr4     {+10134518      >2125>2211_7    AATCCCAATTGATGGAGAATCCCAATTGATGGAGAATCCCAATTGATGGAGAATCCCAATTGATGGAGAATCCCAATTGATGGAGAATCCCAATTGATGGAGAATCCCAATT        CATCCCAATTGATGGAGCATCCCAATTGATGGAGCATCCCAATTGATGGAGCATCCCAATTGATGGAGCATCCCAATTGATGGAGCATCCCAATTGATGGAGCATCCCAATT        60      .       AC=1;AF=0.011236;INV=0;LEN=112;ORIGIN=grch38#chr4:10134484;TYPE=mnp     GT      0|0     0|0     0|0     0|0     0|0     0|0
     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|.     0|0     0|0     0|0     0|0     0|0
     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0
     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|.     0|0     0|0     0|0     0

```

Another problem case fixed with fb7365b07f832bcfcc6ddb693ddfc4a01a25cfbb

```
grch38#chr8_36353854-36453166.vcf:grch38#chr8   36382847        >721>726        GT      GC,AC

SW

vcfallelicprimitives_5.vcf:grch38#chr8  36382847        >721>726_1      GT      AC
vcfallelicprimitives_5.vcf:grch38#chr8  36382848        >721>726_2      T       C

WF

vcfwave_4.vcf:grch38#chr8       36382847        >721>726_1      G       A
```

Let's look at some longer sequences:

The original

```
grch38#chr4_10083863-10181258.vcf:grch38#chr4   10134514        >2136>2148      GGAGAATCCCAATTGATGG     GTAGCATCCCAAGTGATGT,GTAGAATCCCAATTGATGT,GGAGCATCCCAATTGATGG,GG     60      .       AC=11,7,1,3;AF=0.125,0.0795455,0.0113636,0.0340909;AN=88;AT=>2136>2138>2139>2141>2142>2144>2145>2147>2148,>2136>2137>2139>2140>2142>2143>2145>2146>2148,>2136>2137>2139>2141>2142>2144>2145>2146>2148,>2136>2138>2139>2140>2142>2144>2145>2147>2148,>2136>2138>2148;NS=45;LV=1;PS=>2125>2211     GT      0|1     1|0     0|0     0|1     0|0     1|0     1|01|0     0|0     0|0     0|0     0|0     0|0     0|.     0|0     2|2     0|0     4|0     0|0     0|1     0|10|1     0|2     0|0     4|0     0|2     0|0     0|0     2|0     0|0     0|0     0|0     0|0     0|0     2|04|1     0|0     0|0     0|0     0|0     0|3     0|0     0|2     0|0     1
# translates to SW
vcfallelicprimitives_6.vcf:grch38#chr4  10134514        >2136>2148_1    GGAGAATCCCAATTGATG      G       60.AC=3;AF=0.0340909;LEN=17;ORIGIN=grch38#chr4:10134514;TYPE=del   GT      0|0     0|0     0|0     0|0     0|00|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|.     0|0     0|0     0|0     1|0     0|00|0     0|0     0|0     0|0     0|0     1|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|00|0     0|0     1|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0
vcfallelicprimitives_6.vcf:grch38#chr4  10134515        >2136>2148_2    GAGAATCCCAATTGATGG      TAGAATCCCAATTGATGT,TAGCATCCCAAGTGATGT      60      .       AC=7,11;AF=0.0795455,0.125;LEN=18,18;ORIGIN=grch38#chr4:10134514,grch38#chr4:10134514;TYPE=mnp,mnp GT      0|2     2|0     0|0     0|2     0|0     2|0     2|0     2|00|0     0|0     0|0     0|0     0|0     0|.     0|0     1|1     0|0     .|0     0|0     0|2     0|2     0|20|1     0|0     .|0     0|1     0|0     0|0     1|0     0|0     0|0     0|0     0|0     0|0     1|0     .|20|0     0|0     0|0     0|0     0|0     0|0     0|1     0|0     2
vcfallelicprimitives_6.vcf:grch38#chr4  10134518        >2136>2148_3    AATCCCAATTGATGG CATCCCAATTGATGG 60.AC=1;AF=0.0113636;LEN=15;ORIGIN=grch38#chr4:10134514;TYPE=mnp   GT      0|0     0|0     0|0     0|0     0|00|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|.     0|0     0|0     0|0     .|0     0|00|0     0|0     0|0     0|0     0|0     .|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|00|0     0|0     .|0     0|0     0|0     0|0     0|0     0|1     0|0     0|0     0|0     0
# and WFA
vcfallelicprimitives_8.vcf:grch38#chr4  10134515        >2136>2148_1    GAGAATCCCAATTGATGG      TAGCATCCCAAGTGATGG,TAGAATCCCAATTGATGG,G    60      .       AC=11,7,3;AF=0.125,0.0795455,0.0340909;LEN=14,16,17;ORIGIN=grch38#chr4:10134514,grch38#chr4:10134514,grch38#chr4:10134514;TYPE=mnp,mnp,del GT      0|1     1|0     0|00|1     0|0     1|0     1|0     1|0     0|0     0|0     0|0     0|0     0|0     0|.     0|0     2|2     0|03|0     0|0     0|1     0|1     0|1     0|2     0|0     3|0     0|2     0|0     0|0     2|0     0|0     0|00|0     0|0     0|0     2|0     3|1     0|0     0|0     0|0     0|0     0|0     0|0     0|2     0|0     1
vcfallelicprimitives_8.vcf:grch38#chr4  10134518        >2136>2148_2    AATCCCAATTGATG  CATCCCAATTGATG  60.AC=1;AF=0.0113636;LEN=14;ORIGIN=grch38#chr4:10134514;TYPE=mnp   GT      0|0     0|0     0|0     0|0     0|00|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|.     0|0     0|0     0|0     .|0     0|00|0     0|0     0|0     0|0     0|0     .|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|00|0     0|0     .|0     0|0     0|0     0|0     0|0     0|1     0|0     0|0     0|0     0
```

To describe the bug we get for grch38#chr4 10134514 (or: original, sw:, wf:). It was fixed with commit fb7365b07f832bcfcc6ddb693ddfc4a01a25cfbb.

```
ori: GGAGAATCCCAATTGATGG->GG
sw:  GGAGAATCCCAATTGATG ->G
wf:   GAGAATCCCAATTGATGG->G
fix:  GAGAATCCCAATTGATGG->G

ori: GGAGAATCCCAATTGATGG->GTAGCATCCCAAGTGATGT
sw:   GAGAATCCCAATTGATGG-> TAGAATCCCAATTGATGT
wf:   GAGAATCCCAATTGATGG-> TAGAATCCCAATTGATGG <-
fix:  GAGAATCCCAATTGATGG-> TAGAATCCCAATTGATGT

ori: GGAGAATCCCAATTGATGG ->GTAGAATCCCAATTGATGT
sw:   GAGAATCCCAATTGATGG -> TAGCATCCCAAGTGATGT
wf:   GAGAATCCCAATTGATGG -> TAGCATCCCAAGTGATGG <-
fix:   GAGAATCCCAATTGATGG-> TAGCATCCCAAGTGATGT
       GAGAATCCCAATTGATGG-> TAGAATCCCAATTGATGG

ori: GGAGAATCCCAATTGATGG->GGAGCATCCCAATTGATGG
sw:    AATCCCAATTGATGG->      CATCCCAATTGATGG
wf:    AATCCCAATTGATG ->      CATCCCAATTGATG
fix:   AATCCCAATTGATGG->      CATCCCAATTGATGG
       A                      C

```

# vcfwave issue

Now where does the result TAGAATCCCAATTGATGG come from?

```
Original input record (see samples/10134514.vcf)

10134514 GGAGAATCCCAATTGATGG     GTAGCATCCCAAGTGATGT,GTAGAATCCCAATTGATGT,GGAGCATCCCAATTGATGG,GG

WF CIGARs:

10134514:1M1X2M1X7M1X5M1X:GGAGAATCCCAATTGATGG,GTAGCATCCCAAGTGATGT
10134514:1M1X16M1X:GGAGAATCCCAATTGATGG,GTAGAATCCCAATTGATGT
10134514:4M1X14M:GGAGAATCCCAATTGATGG,GGAGCATCCCAATTGATGG
10134514:2M17D:GGAGAATCCCAATTGATGG,GG

Decomposed alleles (return from parsedAlternates):

GGAGAATCCCAATTGATGG 10134514 GGAGAATCCCAATTGATGG -> GGAGAATCCCAATTGATGG
GGAGCATCCCAATTGATGG 10134514 GGAG -> GGAG
GGAGCATCCCAATTGATGG 10134518 A -> C
GGAGCATCCCAATTGATGG 10134519 ATCCCAATTGATGG -> ATCCCAATTGATGG
GTAGAATCCCAATTGATGT 10134514 G -> G
GTAGAATCCCAATTGATGT 10134515 G -> T
GTAGAATCCCAATTGATGT 10134516 AGAATCCCAATTGATG -> AGAATCCCAATTGATG
GTAGAATCCCAATTGATGT 10134532 G -> T
GTAGCATCCCAAGTGATGT 10134514 G -> G
GTAGCATCCCAAGTGATGT 10134515 G -> T
GTAGCATCCCAAGTGATGT 10134516 AG -> AG
GTAGCATCCCAAGTGATGT 10134518 A -> C
GTAGCATCCCAAGTGATGT 10134519 ATCCCAA -> ATCCCAA
GTAGCATCCCAAGTGATGT 10134526 T -> G
GTAGCATCCCAAGTGATGT 10134527 TGATG -> TGATG
GTAGCATCCCAAGTGATGT 10134532 G -> T

Final result (see test/data/regression/vcfwave_5.vcf):

10134515 GAGAATCCCAATTGATGG      TAGAATCCCAATTGATGG,G
10134518 A                       C
10134526 T                       G
10134532 G                       T
```

Now where does TAGAATCCCAATTGATGG come from?

Output produced by test/tests/realign.py



# LICENSE

Copyright 2011-2024 (C) Erik Garrison, Pjotr Prins and vcflib contributors. MIT licensed.
