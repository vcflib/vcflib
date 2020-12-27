% VCFFILTER(1) vcffilter (vcflib) | vcffilter (VCF filter)
% Erik Garrison and vcflib contributors

# NAME

**vcffilter**

# SYNOPSIS

**vcffilter** [options] <vcf file>

# DESCRIPTION

VCF filter the specified vcf file using the set of filters



# OPTIONS

```


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

Filter the specified vcf file using the set of filters.
Filters are specified in the form "<ID> <operator> <value>:
 -f "DP > 10"  # for info fields
 -g "GT = 1|1" # for genotype fields
 -f "CpG"  # for 'flag' fields

Operators can be any of: =, !, <, >, |, &

Any number of filters may be specified.  They are combined via logical AND
unless --or is specified on the command line.  Obtain logical negation through
the use of parentheses, e.g. "! ( DP = 10 )"

For convenience, you can specify "QUAL" to refer to the quality of the site, even
though it does not appear in the INFO fields.

type: filter

```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcffilter.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcffilter.cpp)

# LICENSE

Copyright 2011-2020 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
