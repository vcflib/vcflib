% VCFKEEPSAMPLES(1) vcfkeepsamples (vcflib) | vcfkeepsamples (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcfkeepsamples**

# SYNOPSIS

**vcfkeepsamples** <vcf file> [SAMPLE1] [SAMPLE2] ...

# DESCRIPTION

outputs each record in the vcf file, removing samples not listed on the command line





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcfkeepsamples.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfkeepsamples.cpp)

# LICENSE

Copyright 2011-2024 (C) Erik Garrison and vcflib contributors. MIT licensed.
Copyright 2020-2024 (C) Pjotr Prins.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
