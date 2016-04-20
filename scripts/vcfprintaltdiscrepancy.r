#!/usr/bin/env Rscript --vanilla --slave

# get the input VCF tabular format, assert that sites must have AC > 0
vcf <- subset(read.table(pipe('cat /dev/stdin'), header=T), AC > 0)

tag <- commandArgs(TRUE)[1]

tag.genotypes_alternate_count <- paste(tag, '.genotypes.alternate_count', sep='')
tag.non_reference_discrepancy_count <- paste(tag, '.site.non_reference_discrepancy.count', sep='')
tag.non_reference_discrepancy_normalizer <- paste(tag, '.site.non_reference_discrepancy.normalizer', sep='')
tag.non_reference_sensitivity_count <- paste(tag, '.site.non_reference_sensitivity.count', sep='')
tag.non_reference_sensitivity_normalizer <- paste(tag, '.site.non_reference_sensitivity.normalizer', sep='')
tag.alternate_positive_discrepancy <- paste(tag, '.site.alternate_positive_discrepancy', sep='')
tag.alternate_negative_discrepancy <- paste(tag, '.site.alternate_negative_discrepancy', sep='')
tag.has_variant <- paste(tag, '.has_variant', sep='')

vcf.numberOfSites <- length(vcf[, tag.genotypes_alternate_count])
vcf.totalAltAlleles <- sum(vcf[, tag.genotypes_alternate_count])
vcf.positiveDiscrepancy <- sum(vcf[, tag.alternate_positive_discrepancy]) / sum(vcf[, tag.genotypes_alternate_count])
vcf.negativeDiscrepancy <- sum(vcf[, tag.alternate_negative_discrepancy]) / sum(vcf[, tag.genotypes_alternate_count])
vcf.sitesTruePositive <- sum(vcf[, tag.has_variant]) / nrow(vcf)

cat('number of sites', vcf.numberOfSites, '\n')
cat('total alternate alleles', vcf.totalAltAlleles, '\n')
cat('positive discrepancy', vcf.positiveDiscrepancy, '\n')
cat('negative discrepancy', vcf.negativeDiscrepancy, '\n')

x <- cbind(by(vcf, vcf$AC,
    function(x) {
        sum(x[, tag.alternate_positive_discrepancy]) / sum(x[, tag.genotypes_alternate_count])
    }))

byac <- data.frame(ac=as.numeric(rownames(x)), fdr=as.vector(x))

print(byac)


