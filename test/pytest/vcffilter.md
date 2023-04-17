% VCFFILTER(1) vcffilter (vcflib) | vcffilter (VCF filter)
% Erik Garrison and vcflib contributors

# NAME

**vcffilter**

# SYNOPSIS

**vcffilter** [options] <vcf file>

# DESCRIPTION

VCF filter the specified vcf file using the set of filters.


# OPTIONS
<!--

    >>> from rtest import run_stdout, head, cat, sh

-->

Current command line options:

```

>>> head("vcffilter -h",39)
vcflib filter the specified vcf file using the set of filters
>
usage: vcffilter [options] <vcf file>
>
options:
    -f, --info-filter     specifies a filter to apply to the info fields of records,
                          removes alleles which do not pass the filter
    -g, --genotype-filter specifies a filter to apply to the genotype fields of records
    -k, --keep-info       used in conjunction with '-g', keeps variant info, but removes genotype
    -s, --filter-sites    filter entire records, not just alleles
    -t, --tag-pass        tag vcf records as positively filtered with this tag, print all records
    -F, --tag-fail        tag vcf records as negatively filtered with this tag, print all records
    -A, --append-filter   append the existing filter tag, don't just replace it
    -a, --allele-tag      apply -t on a per-allele basis.  adds or sets the corresponding INFO field tag
    -v, --invert          inverts the filter, e.g. grep -v
    -o, --or              use logical OR instead of AND to combine filters
    -r, --region          specify a region on which to target the filtering, requires a BGZF
                          compressed file which has been indexed with tabix.  any number of
                          regions may be specified.
>
Filter the specified vcf file using the set of filters.
Filters are specified in the form "<ID> <operator> <value>:
 -f "DP > 10"  # for info fields
 -g "GT = 1|1" # for genotype fields
 -f "CpG"  # for 'flag' fields
>
Operators can be any of: =, !, <, >, |, &
>
Any number of filters may be specified.  They are combined via logical AND
unless --or is specified on the command line.  Obtain logical negation through
the use of parentheses, e.g. "! ( DP = 10 )"
>
For convenience, you can specify "QUAL" to refer to the quality of the site, even
though it does not appear in the INFO fields.
>
type: filter
>

```



# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# EXAMPLES

Filter VCF records that have an allele count > 10 results in 1471 matches for

```python
>>> sh("../build/vcffilter -f 'AC > 10' ../samples/grch38#chr4_10083863-10181258.vcf|wc -l")
546

```

# SEE ALSO

[vcflib](./vcflib.md)(1)

Note that [bio-vcf](https://github.com/vcflib/bio-vcf) may easily give 5x better performance. E.g.

```
bio-vcf --filter 'r.info.ac>10'
```

instead of

```
vcfffilter 'AC > 10'
```

will filter all AC fields larger than 10.

# OTHER

## Source code

[vcffilter.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcffilter.cpp)

# LICENSE

Copyright 2011-2023 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
