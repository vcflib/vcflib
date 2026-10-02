% VCFCLASSIFY(1) vcfclassify (vcflib) | vcfclassify (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcfclassify**

# SYNOPSIS

**vcfclassify** <vcf file>

# DESCRIPTION

Creates a new VCF where each variant is tagged by allele class: snp, ts/tv, indel, mnp

# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# EXAMPLES

<!--

    >>> from rtest import run_stdout, head, cat, sh

-->

vcfclassify adds an INFO flag for every allele class in a record. TS and
TV apply to single base substitutions only. A transition is A<->G or
C<->T. Every other single base substitution is a transversion.
Insertions, deletions, MNPs and symbolic alleles are neither.

The test file holds two transitions, two transversions, an insertion, a
deletion, a symbolic allele, an MNP and two multi-allelic records:

```python

>>> sh("grep -v '^##' data/issue_171_vcfclassify.vcf")
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
chr1    10      .       A       G       100     .       .
chr1    20      .       C       T       100     .       .
chr1    30      .       A       C       100     .       .
chr1    40      .       T       G       100     .       .
chr1    50      .       A       AT      100     .       .
chr1    60      .       AT      A       100     .       .
chr1    70      .       A       <DEL>   100     .       .
chr1    80      .       AC      GT      100     .       .
chr1    90      .       A       G,AT    100     .       .
chr1    100     .       A       C,AT    100     .       .

```

Only the two transitions and the multi-allelic record at 90 get the TS flag:

```python

>>> sh("vcfclassify data/issue_171_vcfclassify.vcf|grep -v '^#'|grep TS")
chr1    10      .       A       G       100     .       SNP;TS
chr1    20      .       C       T       100     .       SNP;TS
chr1    90      .       A       G,AT    100     .       INS;SNP;TS

```

Only the two transversions and the multi-allelic record at 100 get the TV
flag. The insertion at 50, the deletion at 60, the symbolic allele at 70
and the MNP at 80 are not single base substitutions, so they get no TV
flag:

```python

>>> sh("vcfclassify data/issue_171_vcfclassify.vcf|grep -v '^#'|grep TV")
chr1    30      .       A       C       100     .       SNP;TV
chr1    40      .       T       G       100     .       SNP;TV
chr1    100     .       A       C,AT    100     .       INS;SNP;TV

```

# SEE ALSO

[vcflib](./vcflib.md)(1)

# OTHER

## Source code

[vcfclassify.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfclassify.cpp)

# LICENSE

Copyright 2011-2026 (C) Erik Garrison and vcflib contributors. MIT licensed.
Copyright 2020-2026 (C) Pjotr Prins.
