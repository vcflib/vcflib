# vcf2tsv

<!--

    >>> from pytest.rtest import run_stdout, head, cat

-->

```

>>> cat("vcf2tsv -h")
usage: vcf2tsv [-n null_string] [-g] [vcf file]
>
Converts stdin or given VCF file to tab-delimited format, using null string to replace empty values in the table.
Specifying -g will output one line per sample with genotype information.
When there is more than one alt allele there will be multiple rows, one for each allele and, the info will match the 'A' index

```

vcf2tsv converts a CSV to a tabulated test file, e.g.

```python

>>> head("vcf2tsv ../samples/sample.vcf")
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  AA      AC      AF      AN      DB      DP      H2      NS
19      111     .       A       C       9.6     .       .       .       .       .       .       .       .       .
19      112     .       A       G       10      .       .       .       .       .       .       .       .       .
20      14370   rs6054257       G       A       29      PASS    .       .       0.5     .       .       14      .       3

```

Use the `-g` switch to show genotypes

```python

>>> head("vcf2tsv -g ../samples/sample.vcf")
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  AA      AC      AF      AN      DB      DP      H2      NS      SAMPLE  DP      GQ      GT      HQ
19      111     .       A       C       9.6     .       .       .       .       .       .       .       .       .       NA00001 .       .       0|0     10,10
19      111     .       A       C       9.6     .       .       .       .       .       .       .       .       .       NA00002 .       .       0|0     10,10
19      111     .       A       C       9.6     .       .       .       .       .       .       .       .       .       NA00003 .       .       0/1     3,3

```

# Source code

[vcf2tsv.cpp](../../src/vcf2tsv.cpp)

# Regression tests

The following commands run full regression tests:

>>> run_stdout("vcf2tsv ../samples/sample.vcf", ext="tsv")
output in <a href="../data/regression/vcf2tsv_4.tsv">vcf2tsv_4.tsv</a>

>>> run_stdout("vcf2tsv -g ../samples/sample.vcf", ext="tsv")
output in <a href="../data/regression/vcf2tsv_5.tsv">vcf2tsv_5.tsv</a>
