% ABBA-BABA(1) abba-baba (vcflib) | abba-baba (VCF genotype)
% Erik Garrison and vcflib contributors

# NAME

abba-baba

# SYNOPSIS

usage: abba-baba --tree 0,1,2,3 --file my.vcf --type PL

# DESCRIPTION

abba-baba calculates the tree pattern for four indviduals. This tool assumes reference is ancestral and ignores non abba-baba sites. The output is a boolian value: 1 = true , 0 = false for abba and baba. the tree argument should be specified from the most basal taxa to the most derived.

# OPTIONS

```


     Example:
     D   C  B   A 
     \ /  /    /  
      \  /    /   
       \    /    
        \  /     
         /        
        /         
 --tree A,B,C,D

Output : 4 columns :     
     1. seqid            
     2. position         
     3. abba             
     4. baba             
required: t,tree       -- a zero based comma separated list of target individuals corrisponding to VCF columns
required: f,file       -- a properly formatted VCF.                                                           
required: y,type       -- genotype likelihood format ; genotypes: GP,GL or PL;                                


type: genotype

```



# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# OTHER

## Source code

[abba-baba.cpp](https://github.com/vcflib/vcflib/blob/master/src/abba-baba.cpp)

# LICENSE

Copyright 2011-2020 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
