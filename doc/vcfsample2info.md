% VCFSAMPLE2INFO(1) vcfsample2info (vcflib) | vcfsample2info (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcfsample2info**

# SYNOPSIS

**vcfsample2info** [options] <vcf file>

# DESCRIPTION

Take annotations given in the per-sample fields and add the mean, median, min, or max to the site-level INFO.



# OPTIONS

```

options:
    -f, --field         Add information about this field in samples to INFO column
    -i, --info          Store the computed statistic in this info field
    -a, --average       Take the mean of samples for field (default)
    -m, --median        Use the median
    -n, --min           Use the min
    -x, --max           Use the max


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

[vcfsample2info.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfsample2info.cpp)

# LICENSE

Copyright 2011-2024 (C) Erik Garrison and vcflib contributors. MIT licensed.
Copyright 2020-2024 (C) Pjotr Prins.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
