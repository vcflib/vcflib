% SEGMENTIHS(1) segmentIhs (vcflib) | segmentIhs (VCF statistics)
% Erik Garrison and vcflib contributors

# NAME

**segmentIhs**

# SYNOPSIS

segmentFst -s 2 -f iHS.normalized.output.txt

# DESCRIPTION

Creates genomic segments (bed file) for regions with high wcFst



# OPTIONS

```

Output : 8 columns :                 
     1. Seqid                        
     2. Start (zero based)           
     3. End   (zero based)           
     4. Average iHS                  
     5. Average high Fst (iHS > -s)  
     6. N iHS values in segment      
     7. N high iHS values in segment 
     8. Segment length               
required: -f            -- Output from normalizeIHS     
optional: -s            -- High absolute iHS cutoff [2] 

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

[segmentIhs.cpp](https://github.com/vcflib/vcflib/blob/master/src/segmentIhs.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
