% VCFANNOTATE(1) vcfannotate (vcflib) | vcfannotate (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcfannotate**

# SYNOPSIS

**vcfannotate** [options] [<vcf file>]

# DESCRIPTION

Intersect the records in the VCF file with targets provided in a BED file. Intersections are done on the reference sequences in the VCF file. If no VCF filename is specified on the command line (last argument) the VCF read from stdin.



# OPTIONS

```


options:
    -b, --bed   use annotations provided by this BED file
    -k, --key   use this INFO field key for the annotations
    -d, --default  use this INFO field key for records without annotations

Type: transformation

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

[vcfannotate.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfannotate.cpp)

# LICENSE

Copyright 2011-2024 (C) Erik Garrison and vcflib contributors. MIT licensed.
Copyright 2020-2024 (C) Pjotr Prins.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
