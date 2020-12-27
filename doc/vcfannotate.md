% VCFANNOTATE(1) vcfannotate (vcflib) | vcfannotate (VCF unknown)
% Erik Garrison and vcflib contributors

# NAME

**vcfannotate**

# SYNOPSIS

**vcfannotate** [options] [<vcf file>]

# DESCRIPTION





# OPTIONS

```

options:
    -b, --bed   use annotations provided by this BED file
    -k, --key   use this INFO field key for the annotations
    -d, --default  use this INFO field key for records without annotations

Intersect the records in the VCF file with targets provided in a BED file.
Intersections are done on the reference sequences in the VCF file.
If no VCF filename is specified on the command line (last argument) the VCF
read from stdin.

```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# OTHER

## Source code

[vcfannotate.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfannotate.cpp)

# LICENSE

Copyright 2011-2020 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
