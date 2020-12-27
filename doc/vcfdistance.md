% VCFDISTANCE(1) vcfdistance (vcflib) | vcfdistance (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcfdistance**

# SYNOPSIS

**vcfdistance** [customtagname] < [vcf file]

# DESCRIPTION

adds a tag to each variant record which indicates the distance to the nearest variant. (defaults to BasesToClosestVariant if no custom tag name is given.



# OPTIONS

```


Type: transformation

```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# OTHER

## Source code

[vcfdistance.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfdistance.cpp)

# LICENSE

Copyright 2011-2020 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
