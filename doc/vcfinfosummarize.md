% VCFINFOSUMMARIZE(1) vcfinfosummarize (vcflib) | vcfinfosummarize (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcfinfosummarize**

# SYNOPSIS

**vcfinfosummarize** [options] <vcf file>

# DESCRIPTION

Take annotations given in the per-sample fields and add the mean, median, min, or max to the site-level INFO.



# OPTIONS

```

options:
    -f, --field         Summarize this field in the INFO column
    -i, --info          Store the computed statistic in this info field
    -a, --average       Take the mean for field (default)
    -m, --median        Use the median
    -n, --min           Use the min
    -x, --max           Use the max
    -h, --help          Print this message
    -v, --version       Print version


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

[vcfinfosummarize.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfinfosummarize.cpp)

# LICENSE

Copyright 2011-2021 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
