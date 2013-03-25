# vcflib
### a C++ library for parsing and manipulating VCF files.

#### author: Erik Garrison <erik.garrison@bc.edu>

#### license: MIT

---

## overview

The Variant Call Format (VCF) is a flat-file, tab-delimited textual format
intended to concisely describe reference-indexed variations between
individuals.  The current specification can be found on the 1000 Genomes wiki
(http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41)

This library provides a simple API to map each VCF record into a C++ class.
See included applications for example usage.

## installing

vcflib includes submodules, so to obtain vcflib you have to use:

    % git clone --recursive git://github.com/ekg/vcflib.git

or

    % git clone --recursive https://github.com/ekg/vcflib.git

To build, use Make:

    % cd vcflib
    % make

Executables are built into the root directory of the repository.

## executables

An overview of functions provided by executables in vcflib.

### vcf2tsv
    
    usage: vcf2tsv [-n null_string] [vcf file]
    Converts stdin or given VCF file to tab-delimited format, using null string to replace empty
    values in the table.

    
### vcfaddinfo
    
    usage: vcfaddinfo <vcf file> <vcf file>
    Adds info fields from the second file which are not present in the first vcf file.
    
    
### vcfafpath

Uses allele frequencies in the AF info column to estimate phylogeny at multiallelic sites.


### vcfallelicprimitives
    
    usage: vcfallelicprimitives [options] [file]
    
    options:
        -m, --use-mnps          Retain MNPs as separate events (default: false)
        -t, --tag-parsed FLAG   Tag records which are split apart of a complex allele
                                with this flag
    
    If multiple alleleic primitives (gaps or mismatches) are specified in
    a single VCF record, split the record into multiple lines, but drop all
    INFO fields.  "Pure" MNPs are split into multiple SNPs unless the -m
    flag is provided.  Genotypes are phased where complex alleles have been
    decomposed, provided genotypes in the input.

    
### vcfaltcount
    
Counts the number of alternate alleles in the record.
    
    
### vcfannotate
    
    usage: vcfannotate [options] [<vcf file>]
    
    options:
        -b, --bed   use annotations provided by this BED file
        -k, --key   use this INFO field key for the annotations
        -d, --default  use this INFO field key for records without annotations
    
    Intersect the records in the VCF file with targets provided in a BED file.
    Intersections are done on the reference sequences in the VCF file.
    If no VCF filename is specified on the command line (last argument) the VCF
    read from stdin.
    
    
### vcfannotategenotypes
    
    usage: vcfannotategenotypes <annotation-tag> <vcf file> <vcf file>
    annotates genotypes in the first file with genotypes in the second
    adding the genotype as another flag to each sample filed in the first file.
    annotation-tag is the name of the sample flag which is added to store the annotation.
    also adds a 'has_variant' flag for sites where the second file has a variant.
    
    
### vcfbreakmulti
    
    usage: vcfbreakmulti [options] [file]
    
    If multiple alleles are specified in a single record, break the record into
    multiple lines, preserving allele-specific INFO fields.
    
    
### vcfcheck
    
    usage: vcfcheck [options] <vcf file>
    
    options:
        -f, --fasta-reference  FASTA reference file to use to obtain primer sequences
    
    Verifies that the VCF REF field matches the reference as described.
    
    
    
### vcfcleancomplex
    
Removes reference-matching sequence from complex alleles and adjusts records to reflect positional change.
    
    
### vcfcommonsamples
    
    usage: vcfcommonsamples <vcf file> <vcf file>
    outputs each record in the first file, removing samples not present in the second
    
    
### vcfcountalleles
    
Counts the total number of alleles in the input.
    
    
### vcfdistance
    
Adds a value to each VCF record indicating the distance to the nearest variant in the file.
    
    
### vcfentropy
    
    usage: vcfentropy [options] <vcf file>
    
    options:
        -f, --fasta-reference  FASTA reference file to use to obtain primer sequences
        -w, --window-size      Size of the window over which to calculate entropy
    
    Anotates the output VCF file with, for each record, EntropyLeft, EntropyRight,
    EntropyCenter, which are the entropies of the sequence of the given window size to the
    left, right, and center  of the record.
    
    
    
### vcffilter
    
    usage: vcffilter [options] <vcf file>
    
    options:
        -f, --info-filter     specifies a filter to apply to the info fields of records,
                              removes alleles which do not pass the filter
        -g, --genotype-filter specifies a filter to apply to the genotype fields of records
        -s, --filter-sites    filter entire records, not just alleles
        -t, --tag             tag vcf records as filtered with this tag instead of suppressing them
        -R, --replace-filter  replace the existing filter tag, don't just append to it
        -a, --allele-tag      apply -t on a per-allele basis.  adds or sets the corresponding INFO field tag
        -v, --invert          inverts the filter, e.g. grep -v
        -o, --or              use logical OR instead of AND to combine filters
        -r, --region          specify a region on which to target the filtering, requires a BGZF
                              compressed file which has been indexed with tabix.  any number of
                              regions may be specified.
    
    Filter the specified vcf file using the set of filters.
    Filters are specified in the form "<ID> <operator> <value>:
     -f "DP > 10"  # for info fields
     -g "GT = 1|1" # for genotype fields
     -f "CpG"  # for 'flag' fields
    
    Operators can be any of: =, !, <, >, |, &
    
    Any number of filters may be specified.  They are combined via logical AND
    unless --or is specified on the command line.  Obtain logical negation through
    the use of parentheses, e.g. ! "( DP = 10 )"
    
    For convenience, you can specify "QUAL" to refer to the quality of the site, even
    though it does not appear in the INFO fields.
    
    
    
### vcffixup
    
Uses genotypes from the VCF file to correct AC, AF, NS, in the VCF records.  For example:

    % vcfkeepsamples file.vcf NA12878 | vcffixup - | vcffilter -f "AC > 0"

Would downsample file.vcf to only NA12878, removing sites for which the sample was not called as polymorphic.
    
    
### vcfflatten
    
    usage: vcfflatten [file]
    
    Removes multi-allelic sites by picking the most common alternate.  Requires
    allele frequency specification 'AF' and use of 'G' and 'A' to specify the
    fields which vary according to the Allele or Genotype. VCF file may be
    specified on the command line or piped as stdin.
    
    
### vcfgeno2haplo
    
    usage: vcfgeno2haplo [options] [<vcf file>]
    
    options:
        -w, --window-size N       compare records up to this many bp away (default 30)
        -r, --reference FILE      FASTA reference file, required with -i and -u
    
    Convert genotype-based phased alleles within --window-size into haplotype alleles.
    
    
    
### vcfgenotypecompare
    
    usage: vcfgenotypecompare <other-genotype-tag> <vcf file>
    adds statistics to the INFO field of the vcf file describing the
    amount of discrepancy between the genotypes (GT) in the vcf file and the
    genotypes reported in the <other-genotype-tag>.  use this after
    vcfannotategenotypes to get correspondence statistics for two vcfs.
    
    
### vcfgenotypes
    
Converts numerical representation of genotypes (standard in GT field) to the alleles provided in the call's ALT/REF fields.
    
    
### vcfglxgt
    
    usage: vcfglxgt [options] <vcf file>
    
    options:
        -n, --fix-null-genotypes   only apply to null and partly-null genotypes
    
    Set genotypes using the maximum genotype likelihood for each sample.
    
    
    
### vcfhetcount
    
Count the number of heterozygotes in the input VCF.
    
    
### vcfhethomratio
    
Provides the ratio between heterozygotes and homozygotes.

    
### vcfintersect
    
    usage: vcfintersect [options] [<vcf file>]
    
    options:
        -b, --bed FILE            use intervals provided by this BED file
        -v, --invert              invert the selection, printing only records which would
                                    not have been printed out
        -i, --intersect-vcf FILE  use this VCF for set intersection generation
        -u, --union-vcf FILE      use this VCF for set union generation
        -w, --window-size N       compare records up to this many bp away (default 30)
        -r, --reference FILE      FASTA reference file, required with -i and -u
        -l, --loci                output whole loci when one alternate allele matches
        -m, --ref-match           intersect on the basis of record REF string
        -t, --tag TAG             attach TAG to each record's info field if it would intersect
        -V, --tag-value VAL       use this value to indicate that the allele is passing
                                  '.' will be used otherwise.  default: 'PASS'
        -M, --merge-from FROM-TAG
        -T, --merge-to   TO-TAG   merge from FROM-TAG used in the -i file, setting TO-TAG
                                  in the current file.
    
    For bed-vcf intersection, alleles which fall into the targets are retained.
    
    For vcf-vcf intersection and union, unify on equivalent alleles within window-size bp
    as determined by haplotype comparison alleles.
    
    
### vcfkeepgeno
    
    usage: vcfkeepgeno <vcf file> [FIELD1] [FIELD2] ...
    outputs each record in the vcf file, removing FORMAT fields not listed
    on the command line from sample specifications in the output
    
    
### vcfkeepinfo
    
    usage: vcfkeepinfo <vcf file> [FIELD1] [FIELD2] ...
    outputs each record in the vcf file, removing INFO fields not listed on the command line
    
    
### vcfkeepsamples
    
    usage: vcfkeepsamples <vcf file> [SAMPLE1] [SAMPLE2] ...
    outputs each record in the vcf file, removing samples not listed on the command line
    
    
### vcflength
    
    terminate called after throwing an instance of 'std::out_of_range'
      what():  basic_string::substr
    
    
### vcfnumalt
    
Annotates the VCF stream on stdin with the number of alternate alleles at the site.
    
    
### vcfoverlay
    
    usage: vcfoverlay [options] [<vcf file> ...]
    
    options:
        -h, --help       this dialog
    
    Overlays records in the input vcf files in the order in which they appear.
    
    
### vcfparsealts
    
Demonstration of alternate allele parsing method.  This method uses pairwise alignment of REF and ALTs to determine component allelic primitives for each alternate allele.

Use `vcfallelicprimitives` to decompose records while preserving format.
    
    
### vcfprimers
    
    usage: vcfprimers [options] <vcf file>
    
    options:
        -f, --fasta-reference  FASTA reference file to use to obtain primer sequences
        -l, --primer-length    The length of the primer sequences on each side of the variant
    
    For each VCF record, extract the flanking sequences, and write them to stdout as FASTA
    records suitable for alignment.  This tool is intended for use in designing validation
    experiments.  Primers extracted which would flank all of the alleles at multi-allelic
    sites.  The name of the FASTA "reads" indicates the VCF record which they apply to.
    The form is >CHROM_POS_LEFT for the 3' primer and >CHROM_POS_RIGHT for the 5' primer,
    for example:
    
    >20_233255_LEFT
    CCATTGTATATATAGACCATAATTTCTTTATCCAATCATCTGTTGATGGA
    >20_233255_RIGHT
    ACTCAGTTGATTCCATACCTTTGCCATCATGAATCATGTTGTAATAAACA
    
    
    
### vcfrandomsample
    
    usage: vcfrandomsample [options] [<vcf file>]
    
    options:
        -r, --rate RATE      base sampling probability per locus
        -s, --scale-by KEY   scale sampling likelihood by this Float info field
        -p, --random-seed N  use this random seed
    
    Randomly sample sites from an input VCF file, which may be provided as stdin.
    Scale the sampling probability by the field specified in KEY.  This may be
    used to provide uniform sampling across allele frequencies, for instance.
    
    
### vcfremap
    
    usage: vcfremap [options] [<vcf file>]
    
    options:
        -w, --ref-window-size N      align using this many bases flanking each side of the reference allele
        -s, --alt-window-size N      align using this many flanking bases from the reference around each alternate allele
        -r, --reference FILE         FASTA reference file, required with -i and -u
        -m, --match-score N          match score for SW algorithm
        -x, --mismatch-score N       mismatch score for SW algorithm
        -o, --gap-open-penalty N     gap open penalty for SW algorithm
        -e, --gap-extend-penalty N   gap extension penalty for SW algorithm
        -z, --entropy-gap-open       use entropy scaling for the gap open penalty
        -R, --repeat-gap-extend N    penalize non-repeat-unit gaps in repeat sequence
        -a, --adjust-vcf TAG         supply a new cigar as TAG in the output VCF
    
    For each alternate allele, attempt to realign against the reference with lowered gap open penalty.
    If realignment is possible, adjust the cigar and reference/alternate alleles.
    
    
### vcfremoveaberrantgenotypes
    
Strips genotypes which are homozygous but have observations implying heterozygosity.  Requires RA (reference allele observation) and AA (alternate allele observation) for each genotype.
    
    
### vcfremovesamples
    
    usage: vcfremovesamples <vcf file> [SAMPLE1] [SAMPLE2] ...
    outputs each record in the vcf file, removing samples listed on the command line
    
    
### vcfroc
    
    usage: vcfroc [options] [<vcf file>]
    
    options:
        -t, --truth-vcf FILE      use this VCF as ground truth for ROC generation
        -w, --window-size N       compare records up to this many bp away (default 30)
        -r, --reference FILE      FASTA reference file
    
    Generates a pseudo-ROC curve using sensitivity and specificity estimated against
    a putative truth set.  Thresholding is provided by successive QUAL cutoffs.
    
    
### vcfsamplediff
    
    usage: vcfsamplediff <tag> <sample> <sample> [ <sample> ... ] <vcf file>
    tags each record where the listed sample genotypes differ with <tag>
    The first sample is assumed to be germline, the second somatic.
    Each record is tagged with <tag>={germline,somatic,loh} to specify the type of
    variant given the genotype difference between the two samples.
    
    
### vcfsamplenames
    
Prints the names of the samples in the VCF file.

    
### vcfsom
    
    usage: vcfsom [options] [vcf file]
    
    training: 
        vcfsom -s output.som -f "AF DP ABP" training.vcf
    
    application: 
        vcfsom -a output.som -f "AF DP ABP" test.vcf >results.vcf
    
    vcfsomtrains and/or applies a self-organizing map to the input VCF data
    on stdin, adding two columns for the x and y coordinates of the winning
    neuron in the network and an optional euclidean distance from a given
    node (--center).
    
    If a map is provided via --apply,  map will be applied to input without
    training.  Automated filtering to an estimated FP rate is 
    
    options:
    
        -h, --help             this dialog
    
    training:
    
        -f, --fields "FIELD ..."  INFO fields to provide to the SOM
        -a, --apply FILE       apply the saved map to input data to FILE
        -s, --save  FILE       train on input data and save the map to FILE
        -t, --print-training-results
                               print results of SOM on training input
                               (you can also just use --apply on the same input)
        -x, --width X          width in columns of the output array
        -y, --height Y         height in columns of the output array
        -i, --iterations N     number of training iterations or epochs
        -d, --debug            print timing information
    
    recalibration:
    
        -c, --center X,Y       annotate with euclidean distance from center
        -p, --paint-true VCF   use VCF file to annotate true variants (multiple)
        -f, --paint-false VCF  use VCF file to annotate false variants (multiple)
        -R, --paint-tag TAG    provide estimated FDR% in TAG in variant INFO
        -N, --false-negative   replace FDR% (false detection) with FNR% (false negative)
    
    
### vcfstats
    
    usage: vcfstats [options] <vcf file>
    
        -r, --region          specify a region on which to target the stats, requires a BGZF
                              compressed file which has been indexed with tabix.  any number of
                              regions may be specified.
        -a, --add-info        add the statistics intermediate information to the VCF file,
                              writing out VCF records instead of summary statistics
        -l, --no-length-frequency    don't out the indel and mnp length-frequency spectra
        -m, --match-score N          match score for SW algorithm
        -x, --mismatch-score N       mismatch score for SW algorithm
        -o, --gap-open-penalty N     gap open penalty for SW algorithm
        -e, --gap-extend-penalty N   gap extension penalty for SW algorithm
    
    Prints statistics about variants in the input VCF file.
    
    
### vcfstreamsort
    
Reads VCF on stdin and guarantees that the positional order is correct provided out-of-order
variants are no more than 100bp apart.
    
    
