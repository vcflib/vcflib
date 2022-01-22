% VCFSITESUMMARIZE(1) vcfsitesummarize (vcflib) | vcfsitesummarize (VCF statistics)
% Erik Garrison and vcflib contributors

# NAME

**vcfsitesummarize**

# SYNOPSIS

**vcfsitesummarize** <vcf file>

# DESCRIPTION

Summarize by site





# EXAMPLES

```

Example:

**vcfsitesummarize** samples/sample.vcf

CHROM   POS     ID      REF     QUAL    FILTER  AA      AC      AF      AN      DP      NS      DB      H2
19      111     .       A       9.6     .                                                       0       0
19      112     .       A       10      .                                                       0       0
20      14370   rs6054257       G       29      PASS                    0.5             14      3       1 1
20      17330   .       T       3       q10                     0.017           11      3       0       0
20      1110696 rs6040355       A       67      PASS    T                               10      2       1 0
20      1230237 .       T       47      PASS    T                               13      3       0       0
20      1234567 microsat1       G       50      PASS    G                       6       9       3       0 0
20      1235237 .       T       0       .                                                       0       0
X       10      rsTest  AC      10      PASS


Type: statistics

      

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

[vcfsitesummarize.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfsitesummarize.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
