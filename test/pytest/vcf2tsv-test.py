import rtest

datadir = "../samples"

vcf = datadir+"/sample.vcf"

def test01(n):
    """vcf2tsv converts a CSV to a tabulated test file, e.g.

       vcf2tsv samples/sample.vcf

    outputs

```
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  AA      AC      AF      AN      DB      DP      H2NS
19      111     .       A       C       9.6     .       .       .       .       .       .       .       . .
19      112     .       A       G       10      .       .       .       .       .       .       .       . .
20      14370   rs6054257       G       A       29      PASS    .       .       0.5     .       .       14.3
```

    >>> rtest.run_stdout(f"vcf2tsv {vcf}", ext="tsv")

    """
    pass


if __name__ == "__main__":
    import doctest
    rtest.setup()
    doctest.testmod()
