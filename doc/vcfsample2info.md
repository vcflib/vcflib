% VCFSAMPLE2INFO(1) vcfsample2info (vcflib) | vcfsample2info (VCF unknown)
% Erik Garrison and vcflib contributors

# NAME

**vcfsample2info**

# SYNOPSIS

./build/**vcfsample2info** [options] <vcf file>

# DESCRIPTION

options: -f, --field Add information about this field in samples to INFO column -i, --info Store the computed statistic in this info field -a, --average Take the mean of samples for field (default) -m, --median Use the median -n, --min Use the min -x, --max Use the max



# OPTIONS

```


Take annotations given in the per-sample fields and add the mean, median, min, or max
to the site-level INFO.

```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# OTHER

## Source code

[vcfsample2info.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfsample2info.cpp)

# LICENSE

Copyright 2011-2020 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
