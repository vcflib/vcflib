% vcflib(1) vcflib | vcflib (index)
% Erik Garrison and vcflib contributors

# NAME

**vcflib** index

# DESCRIPTION

vcflib contains tools and libraries for dealing with the Variant Call
Format (VCF) which is a flat-file, tab-delimited textual format
intended to describe reference-indexed variations between
individuals.

VCF provides a common interchange format for the description of
variation in individuals and populations of samples, and has become
the defacto standard reporting format for a wide array of genomic
variant detectors.

vcflib provides methods to manipulate and interpret sequence variation
as it can be described by VCF. It is both:

* an API for parsing and operating on records of genomic variation as it can be described by the VCF format,
* and a collection of command-line utilities for executing complex manipulations on VCF files.

The API itself provides a quick and extremely permissive method to
read and write VCF files. Extensions and applications of the library
provided in the included utilities (*.cpp) comprise the vast bulk of
the library's utility for most users.

<!--
  Created with ./scripts/bin2md.rb --index
-->


## filter

| filter command | description |
| :-------------- | :---------- |
 | [vcffilter](./vcffilter.md) | VCF filter the specified vcf file using the set of filters |
 | [vcfuniq](./vcfuniq.md) | List unique genotypes. Similar to GNU uniq, but aimed at VCF records. **vcfuniq** removes records which have the same position, ref, and alt as the previous record on a sorted VCF file. Note that it does not adjust/combine genotypes in the output, but simply takes the first record. See also vcfcreatemulti for combining records. |
 | [vcfuniqalleles](./vcfuniqalleles.md) | List unique alleles For each record, remove any duplicate alternate alleles that may have resulted from merging separate VCF files. |

## metrics

| metrics command | description |
| :-------------- | :---------- |
 | [vcfcheck](./vcfcheck.md) | Validate integrity and identity of the VCF by verifying that the VCF record's REF matches a given reference file. |
 | [vcfdistance](./vcfdistance.md) | Adds a tag to each variant record which indicates the distance to the nearest variant. (defaults to BasesToClosestVariant if no custom tag name is given. |
 | [vcfentropy](./vcfentropy.md) | Annotate VCF records with the Shannon entropy of flanking sequence. Anotates the output VCF file with, for each record, EntropyLeft, EntropyRight, EntropyCenter, which are the entropies of the sequence of the given window size to the left, right, and center of the record. Also adds EntropyRef and EntropyAlt for each alt. |
 | [vcfhetcount](./vcfhetcount.md) | Calculate the heterozygosity rate: count the number of alternate alleles in heterozygous genotypes in all records in the vcf file |
 | [vcfhethomratio](./vcfhethomratio.md) | Generates the het/hom ratio for each individual in the file |

## phenotype

| phenotype command | description |
| :-------------- | :---------- |
 | [permuteGPAT++](./permuteGPAT++.md) | **permuteGPAT++** is a method for adding empirical p-values to a GPAT++ score. |

## genotype

| genotype command | description |
| :-------------- | :---------- |
 | [abba-baba](./abba-baba.md) | **abba-baba** calculates the tree pattern for four indviduals. This tool assumes reference is ancestral and ignores non **abba-baba** sites. The output is a boolian value: 1 = true , 0 = false for abba and baba. the tree argument should be specified from the most basal taxa to the most derived. |
 | [hapLrt](./hapLrt.md) | HapLRT is a likelihood ratio test for haplotype lengths. The lengths are modeled with an exponential distribution. The sign denotes if the target has longer haplotypes (1) or the background (-1). |
 | [normalize-iHS](./normalize-iHS.md) | normalizes iHS or XP-EHH scores. |

## transformation

| transformation command | description |
| :-------------- | :---------- |
 | [dumpContigsFromHeader](./dumpContigsFromHeader.md) | Dump contigs from header |
 | [smoother](./smoother.md) | smoothes is a method for window smoothing many of the GPAT++ formats. |
 | [vcf2dag](./vcf2dag.md) | Modify VCF to be able to build a directed acyclic graph (DAG) |
 | [vcf2fasta](./vcf2fasta.md) | Generates sample_seq:N.fa for each sample, reference sequence, and chromosomal copy N in [0,1... ploidy]. Each sequence in the fasta file is named using the same pattern used for the file name, allowing them to be combined. |
 | [vcf2tsv](./vcf2tsv.md) | Converts VCF to per-allelle or per-genotype tab-delimited format, using null string to replace empty values in the table. Specifying -g will output one line per sample with genotype information. When there is more than one alt allele there will be multiple rows, one for each allele and, the info will match the 'A' index |
 | [vcfaddinfo](./vcfaddinfo.md) | Adds info fields from the second file which are not present in the first vcf file. |
 | [vcfafpath](./vcfafpath.md) | Display genotype paths |
 | [vcfallelicprimitives](./vcfallelicprimitives.md) | WARNING: this tool is considered legacy and is only retained for older workflows. It will emit a warning! Even though it can use the WFA you should use [vcfwave](./vcfwave.md) instead. |
 | [vcfannotate](./vcfannotate.md) | Intersect the records in the VCF file with targets provided in a BED file. Intersections are done on the reference sequences in the VCF file. If no VCF filename is specified on the command line (last argument) the VCF read from stdin. |
 | [vcfannotategenotypes](./vcfannotategenotypes.md) | Examine genotype correspondence. Annotate genotypes in the first file with genotypes in the second adding the genotype as another flag to each sample filed in the first file. annotation-tag is the name of the sample flag which is added to store the annotation. also adds a 'has_variant' flag for sites where the second file has a variant. |
 | [vcfbreakmulti](./vcfbreakmulti.md) | If multiple alleles are specified in a single record, break the record into multiple lines, preserving allele-specific INFO fields. |
 | [vcfcat](./vcfcat.md) | Concatenates VCF files |
 | [vcfclassify](./vcfclassify.md) | Creates a new VCF where each variant is tagged by allele class: snp, ts/tv, indel, mnp |
 | [vcfcleancomplex](./vcfcleancomplex.md) | Removes reference-matching sequence from complex alleles and adjusts records to reflect positional change. |
 | [vcfcombine](./vcfcombine.md) | Combine VCF files positionally, combining samples when sites and alleles are identical. Any number of VCF files may be combined. The INFO field and other columns are taken from one of the files which are combined when records in multiple files match. Alleles must have identical ordering to be combined into one record. If they do not, multiple records will be emitted. |
 | [vcfcommonsamples](./vcfcommonsamples.md) | Generates each record in the first file, removing samples not present in the second |
 | [vcfcreatemulti](./vcfcreatemulti.md) | Go through sorted VCF and when overlapping alleles are represented across multiple records, merge them into a single multi-ALT record. See the documentation for more information. |
 | [vcfecho](./vcfecho.md) | Echo VCF to stdout (simple demo) |
 | [vcfevenregions](./vcfevenregions.md) | Generates a list of regions, e.g. chr20:10..30 using the variant density information provided in the VCF file to ensure that the regions have even numbers of variants. This can be use to reduce the variance in runtime when dividing variant detection or genotyping by genomic coordinates. |
 | [vcffixup](./vcffixup.md) | Generates a VCF stream where AC and NS have been generated for each record using sample genotypes |
 | [vcfflatten](./vcfflatten.md) | Removes multi-allelic sites by picking the most common alternate. Requires allele frequency specification 'AF' and use of 'G' and 'A' to specify the fields which vary according to the Allele or Genotype. VCF file may be specified on the command line or piped as stdin. |
 | [vcfgeno2alleles](./vcfgeno2alleles.md) | modifies the genotypes field to provide the literal alleles rather than indexes |
 | [vcfgeno2haplo](./vcfgeno2haplo.md) | Convert genotype-based phased alleles within --window-size into haplotype alleles. Will break haplotype construction when encountering non-phased genotypes on input. |
 | [vcfgenosamplenames](./vcfgenosamplenames.md) | Get samplenames |
 | [vcfglbound](./vcfglbound.md) | Adjust GLs so that the maximum GL is 0 by dividing all GLs for each sample by the max. |
 | [vcfglxgt](./vcfglxgt.md) | Set genotypes using the maximum genotype likelihood for each sample. |
 | [vcfindex](./vcfindex.md) | Adds an index number to the INFO field (id=position) |
 | [vcfinfo2qual](./vcfinfo2qual.md) | Sets QUAL from info field tag keyed by [key]. The VCF file may be omitted and read from stdin. The average of the field is used if it contains multiple values. |
 | [vcfinfosummarize](./vcfinfosummarize.md) | Take annotations given in the per-sample fields and add the mean, median, min, or max to the site-level INFO. |
 | [vcfintersect](./vcfintersect.md) | VCF set analysis |
 | [vcfkeepgeno](./vcfkeepgeno.md) | Reduce file size by removing FORMAT fields not listed on the command line from sample specifications in the output |
 | [vcfkeepinfo](./vcfkeepinfo.md) | To decrease file size remove INFO fields not listed on the command line |
 | [vcfkeepsamples](./vcfkeepsamples.md) | outputs each record in the vcf file, removing samples not listed on the command line |
 | [vcfld](./vcfld.md) | Compute LD |
 | [vcfleftalign](./vcfleftalign.md) | Left-align indels and complex variants in the input using a pairwise ref/alt alignment followed by a heuristic, iterative left realignment process that shifts indel representations to their absolute leftmost (5') extent. |
 | [vcflength](./vcflength.md) | Add length info field |
 | [vcfnullgenofields](./vcfnullgenofields.md) | Makes the FORMAT for each variant line the same (uses all the FORMAT fields described in the header). Fills out per-sample fields to match FORMAT. Expands GT values of '.' with number of alleles based on ploidy (eg: './.' for dipolid). |
 | [vcfnumalt](./vcfnumalt.md) | outputs a VCF stream where NUMALT has been generated for each record using sample genotypes |
 | [vcfoverlay](./vcfoverlay.md) | Overlay records in the input vcf files with order as precedence. |
 | [vcfprimers](./vcfprimers.md) | For each VCF record, extract the flanking sequences, and write them to stdout as FASTA records suitable for alignment. |
 | [vcfqual2info](./vcfqual2info.md) | Puts QUAL into an info field tag keyed by [key]. |
 | [vcfremap](./vcfremap.md) | For each alternate allele, attempt to realign against the reference with lowered gap open penalty. If realignment is possible, adjust the cigar and reference/alternate alleles. Observe how different alignment parameters, including context and entropy-dependent ones, influence variant classification and interpretation. |
 | [vcfremoveaberrantgenotypes](./vcfremoveaberrantgenotypes.md) | strips samples which are homozygous but have observations implying heterozygosity. Remove samples for which the reported genotype (GT) and observation counts disagree (AO, RO). |
 | [vcfremovesamples](./vcfremovesamples.md) | outputs each record in the vcf file, removing samples listed on the command line |
 | [vcfsample2info](./vcfsample2info.md) | Take annotations given in the per-sample fields and add the mean, median, min, or max to the site-level INFO. |
 | [vcfsamplediff](./vcfsamplediff.md) | Establish putative somatic variants using reported differences between germline and somatic samples. Tags each record where the listed sample genotypes differ with <tag>. The first sample is assumed to be germline, the second somatic. Each record is tagged with <tag>={germline,somatic,loh} to specify the type of variant given the genotype difference between the two samples. |
 | [vcfsamplenames](./vcfsamplenames.md) | List sample names |
 | [vcfstreamsort](./vcfstreamsort.md) | Sorts the input (either stdin or file) using a streaming sort algorithm. Guarantees that the positional order is correct provided out-of-order variants are no more than 100 positions in the VCF file apart. |
 | [vcfwave](./vcfwave.md) | Realign reference and alternate alleles with WFA, parsing out the 'primitive' alleles into multiple VCF records. New records have IDs that reference the source record ID. Genotypes/samples are handled correctly. Deletions generate haploid/missing genotypes at overlapping sites. |

## statistics

| statistics command | description |
| :-------------- | :---------- |
 | [bFst](./bFst.md) | **bFst** is a Bayesian approach to Fst. Importantly **bFst** accounts for genotype uncertainty in the model using genotype likelihoods. For a more detailed description see: `A Bayesian approach to inferring population structure from dominant markers' by Holsinger et al. Molecular Ecology Vol 11, issue 7 2002. The likelihood function has been modified to use genotype likelihoods provided by variant callers. There are five free parameters estimated in the model: each subpopulation's allele frequency and Fis (fixation index, within each subpopulation), a free parameter for the total population's allele frequency, and Fst. |
 | [genotypeSummary](./genotypeSummary.md) | Generates a table of genotype counts. Summarizes genotype counts for bi-allelic SNVs and indel |
 | [iHS](./iHS.md) | **iHS** calculates the integrated haplotype score which measures the relative decay of extended haplotype homozygosity (EHH) for the reference and alternative alleles at a site (see: voight et al. 2006, Spiech & Hernandez 2014). |
 | [meltEHH](./meltEHH.md) |  |
 | [pFst](./pFst.md) | **pFst** is a probabilistic approach for detecting differences in allele frequencies between two populations. |
 | [pVst](./pVst.md) | **pVst** calculates vst, a measure of CNV stratification. |
 | [permuteSmooth](./permuteSmooth.md) | **permuteSmooth** is a method for adding empirical p-values smoothed wcFst scores. |
 | [plotHaps](./plotHaps.md) | **plotHaps** provides the formatted output that can be used with 'bin/plotHaplotypes.R'. |
 | [popStats](./popStats.md) | General population genetic statistics for each SNP |
 | [segmentFst](./segmentFst.md) | **segmentFst** creates genomic segments (bed file) for regions with high wcFst |
 | [segmentIhs](./segmentIhs.md) | Creates genomic segments (bed file) for regions with high wcFst |
 | [sequenceDiversity](./sequenceDiversity.md) | The **sequenceDiversity** program calculates two popular metrics of haplotype diversity: pi and extended haplotype homozygoisty (eHH). Pi is calculated using the Nei and Li 1979 formulation. eHH a convenient way to think about haplotype diversity. When eHH = 0 all haplotypes in the window are unique and when eHH = 1 all haplotypes in the window are identical. |
 | [vcfaltcount](./vcfaltcount.md) | count the number of alternate alleles in all records in the vcf file |
 | [vcfcountalleles](./vcfcountalleles.md) | Count alleles |
 | [vcfgenosummarize](./vcfgenosummarize.md) | Adds summary statistics to each record summarizing qualities reported in called genotypes. Uses: RO (reference observation count), QR (quality sum reference observations) AO (alternate observation count), QA (quality sum alternate observations) |
 | [vcfgenotypecompare](./vcfgenotypecompare.md) | adds statistics to the INFO field of the vcf file describing the amount of discrepancy between the genotypes (GT) in the vcf file and the genotypes reported in the <other-genotype-tag>. use this after vcfannotategenotypes to get correspondence statistics for two vcfs. |
 | [vcfgenotypes](./vcfgenotypes.md) | Report the genotypes for each sample, for each variant in the VCF. Convert the numerical represenation of genotypes provided by the GT field to a human-readable genotype format. |
 | [vcfparsealts](./vcfparsealts.md) | Alternate allele parsing method. This method uses pairwise alignment of REF and ALTs to determine component allelic primitives for each alternate allele. |
 | [vcfrandom](./vcfrandom.md) | Generate a random VCF file |
 | [vcfrandomsample](./vcfrandomsample.md) | Randomly sample sites from an input VCF file, which may be provided as stdin. Scale the sampling probability by the field specified in KEY. This may be used to provide uniform sampling across allele frequencies, for instance. |
 | [vcfroc](./vcfroc.md) | Generates a pseudo-ROC curve using sensitivity and specificity estimated against a putative truth set. Thresholding is provided by successive QUAL cutoffs. |
 | [vcfsitesummarize](./vcfsitesummarize.md) | Summarize by site |
 | [vcfstats](./vcfstats.md) | Prints statistics about variants in the input VCF file. |
 | [wcFst](./wcFst.md) | **wcFst** is Weir & Cockerham's Fst for two populations. Negative values are VALID, they are sites which can be treated as zero Fst. For more information see Evolution, Vol. 38 N. 6 Nov 1984. Specifically **wcFst** uses equations 1,2,3,4. |

# SOURCE CODE

See the source code repository at https://github.com/vcflib/vcflib

# CREDIT

Citations are the bread and butter of Science.  If you are using this
software in your research and want to support our future work, please
cite the following publication:

Please cite:

[A spectrum of free software tools for processing the VCF variant call format: vcflib, bio-vcf, cyvcf2, hts-nim and slivar](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009123).
Garrison E, Kronenberg ZN, Dawson ET, Pedersen BS, Prins P (2022), PLoS Comput Biol 18(5): e1009123. https://doi.org/10.1371/journal.pcbi.1009123


# LICENSE

Copyright 2011-2023 (C) Erik Garrison and vcflib contributors. MIT licensed.

