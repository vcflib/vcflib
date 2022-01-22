% VCFGENOSUMMARIZE(1) vcfgenosummarize (vcflib) | vcfgenosummarize (VCF statistics)
% Erik Garrison and vcflib contributors

# NAME

**vcfgenosummarize**

# SYNOPSIS

**vcfgenosummarize** <[input file] >[output vcf]

# DESCRIPTION

Adds summary statistics to each record summarizing qualities reported in called genotypes. Uses: RO (reference observation count), QR (quality sum reference observations) AO (alternate observation count), QA (quality sum alternate observations)



# OPTIONS

```


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

[vcfgenosummarize.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfgenosummarize.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
