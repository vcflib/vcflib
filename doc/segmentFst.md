% SEGMENTFST(1) segmentFst (vcflib) | segmentFst (VCF statistics)
% Erik Garrison and vcflib contributors

# NAME

**segmentFst**

# SYNOPSIS

**segmentFst** -s 0.7 -f wcFst.output.txt

# DESCRIPTION

**segmentFst** creates genomic segments (bed file) for regions with high wcFst



# OPTIONS

```


**segmentFst** provides a way to find continious regions with high Fst values.  It takes the output of wcFst and produces a BED file.  These high Fst region can be permutated with 'permuteGPATwindow'
Output : 8 columns :                 
     1. Seqid                        
     2. Start (zero based)           
     3. End   (zero based)           
     4. Average Fst                  
     5. Average high Fst (Fst > -s)  
     6. N Fst values in segment      
     7. N high fst values in segment 
     8. Segment length               
required: -f            -- Output from wcFst     
optional: -s            -- High Fst cutoff [0.8] 

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

[segmentFst.cpp](https://github.com/vcflib/vcflib/blob/master/src/segmentFst.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
