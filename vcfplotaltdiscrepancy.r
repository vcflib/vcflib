#!/usr/bin/Rscript --vanilla --slave

# get the input VCF tabular format, assert that sites must have AC > 0
vcf <- subset(read.table(pipe('cat /dev/stdin'), header=T), AC > 0)

filename <- commandArgs(TRUE)[1]
tag <- commandArgs(TRUE)[2]

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

#cat('number of sites', vcf.numberOfSites, '\n')
#cat('total alternate alleles', vcf.totalAltAlleles, '\n')
#cat('positive discrepancy', vcf.positiveDiscrepancy, '\n')
#cat('negative discrepancy', vcf.negativeDiscrepancy, '\n')

x <- cbind(by(vcf, vcf$AC,
    function(x) {
        sum(x[, tag.alternate_positive_discrepancy]) / sum(x[, tag.genotypes_alternate_count])
    }))

byac <- data.frame(ac=as.numeric(rownames(x)), fdr=as.vector(x))

# false detection count
byac$fpc <- as.vector(cbind(by(vcf, vcf$AC,
    function(i) { sum(i[, tag.alternate_positive_discrepancy]) } )))

byac$alleles <- as.vector(cbind(by(vcf, vcf$AC,
    function(i) {
        sum(i[, tag.genotypes_alternate_count])
    })))

byac$sites <- as.vector(cbind(by(vcf$AC, vcf$AC, function(i) length(i))))

# count true positive sites
byac$site_tpc <- as.vector(cbind(by(vcf[, tag.has_variant], vcf$AC, function(i) sum(i))))
# fpc == false detection count
byac$site_fpc <- byac$sites - byac$site_tpc
# site detection fpr is 1 - true positive rate
byac$site_fpr <- 1 - ( byac$site_tpc / byac$sites )

#byac$site_fprlt <- as.vector(cbind(tapply(byac$ac, byac$ac, function(i) mean(subset(byac, ac <= i, select=site_fpr)))))
byac$site_fprlt <- as.vector(cbind(tapply(byac$ac, byac$ac, function(i) { 
    s <- subset(byac, ac <= i, select=c(site_fpc, sites))
    return(sum(s$site_fpc) / sum(s$sites))
})))

#byac$site_fprgt <- as.vector(cbind(tapply(byac$ac, byac$ac, function(i) mean(subset(byac, ac >= i, select=site_fpr)))))
byac$site_fprgt <- as.vector(cbind(tapply(byac$ac, byac$ac, function(i) { 
    s <- subset(byac, ac >= i, select=c(site_fpc, sites))
    return(sum(s$site_fpc) / sum(s$sites))
})))

byac$cfa <- as.vector(cbind(tapply(byac$ac, byac$ac, function(i) sum(subset(byac, ac <= i, select=alleles)) / sum(byac$alleles))))

byac$cfs <- as.vector(cbind(tapply(byac$ac, byac$ac, function(i) sum(subset(byac, ac <= i, select=sites)) / length(vcf$AC))))

# inappropriate collapse via averaging of fdr
#byac$alternate_pdlt <- as.vector(cbind(tapply(byac$ac, byac$ac, function(i) mean(subset(byac, ac <= i, select=fdr)))))
# use this one
byac$alternate_pdlt <- as.vector(cbind(tapply(byac$ac, byac$ac, function(i) { 
    s <- subset(byac, ac <= i, select=c(fpc, alleles))
    return(sum(s$fpc) / sum(s$alleles))
})))

#byac$alternate_pdgt <- as.vector(cbind(tapply(byac$ac, byac$ac, function(i) mean(subset(byac, ac >= i, select=fdr)))))
byac$alternate_pdgt <- as.vector(cbind(tapply(byac$ac, byac$ac, function(i) { 
    s <- subset(byac, ac >= i, select=c(fpc, alleles))
    return(sum(s$fpc) / sum(s$alleles))
})))

nan.to.zero <- function(n) {
    if (is.nan(n)) return(0) else return(n)
}

byac$nrs <- as.vector(cbind(by(vcf, vcf$AC, function(i) {
    return(nan.to.zero(sum(i[, tag.non_reference_sensitivity_count]) / sum(i[, tag.non_reference_sensitivity_normalizer])))
})))

byac$nrd <- as.vector(cbind(by(vcf, vcf$AC, function(i) { 
    return(nan.to.zero(sum(i[, tag.non_reference_discrepancy_count]) / sum(i[, tag.non_reference_discrepancy_normalizer])))
})))

byac$nrslt <- as.vector(cbind(tapply(byac$ac, byac$ac, function(i) {
    s <- subset(vcf, AC <= i, select=c(tag.non_reference_sensitivity_count, tag.non_reference_sensitivity_normalizer))
    return(nan.to.zero(sum(s[, tag.non_reference_sensitivity_count]) / sum(s[, tag.non_reference_sensitivity_normalizer])))
})))

byac$nrdlt <- as.vector(cbind(tapply(byac$ac, byac$ac, function(i) { 
    s <- subset(vcf, AC <= i, select=c(tag.non_reference_discrepancy_count, tag.non_reference_discrepancy_normalizer))
    return(nan.to.zero(sum(s[, tag.non_reference_discrepancy_count]) / sum(s[, tag.non_reference_discrepancy_normalizer])))
})))


pdf(paste(filename, '.', tag, '.FDR.vs.AC.pdf', sep=''))
par(cex=0.75)
par(mar=c(5,4,4,5) + 0.1)
plot(byac$alternate_pdgt, ylim=c(0,1.0),
    xlab='alternate allele count', xaxt='n',
    ylab='', yaxt='n', type='l', col='blue')
axis(2, at=seq(0,1,0.1), labels=seq(0,1,0.1))
axis(1, at=seq(0,max(byac$ac),10), labels=seq(0,max(byac$ac),10), cex=0.75)
grid(lty=5)
par(new=T)
title(paste(filename, 'putative false discovery rate versus', tag))
par(new=T)
countTicks <- seq(0,1,0.1) * vcf.numberOfSites
axis(4, at=seq(0,1,0.1), labels=round(countTicks))
mtext("number of sites", side=4, line=3, cex=0.75)
par(new=T)
plot(byac$alternate_pdlt, ylim=c(0,1.0), xlab='', xaxt='n', ylab='', yaxt='n', type='l', col='green')
par(new=T)
plot(byac$cfa, ylim=c(0,1.0), xlab='', xaxt='n', ylab='', yaxt='n', type='l', col='red')
par(new=T)
plot(byac$cfs, ylim=c(0,1.0), xlab='', xaxt='n', ylab='', yaxt='n', type='l', col='purple')
par(new=T)
plot(byac$site_fprgt, ylim=c(0,1.0), xlab='', xaxt='n', ylab='', yaxt='n', type='l', col='magenta')
par(new=T)
plot(byac$site_fprlt, ylim=c(0,1.0), xlab='', xaxt='n', ylab='', yaxt='n', type='l', col='cyan')
par(new=T)
plot(byac$nrslt, ylim=c(0,1.0), xlab='', xaxt='n', ylab='', yaxt='n', type='l', col='yellow')
par(new=T)
plot(byac$nrdlt, ylim=c(0,1.0), xlab='', xaxt='n', ylab='', yaxt='n', type='l', col='brown')
#par(new=T)
#plot(vcf$AC, vcf[, paste(tag, '.site.alternate_positive_discrepancy', sep='')] / vcf[, paste(tag, '.genotypes.alternate_count', sep='')],
#    pch=16, cex=0.5, col='grey',
#    ylim=c(0,1.0), xlab='', xaxt='n', ylab='', yaxt='n',)
par(new=T, cex=0.65)
legend('topleft', c('cumulative alt alleles', 'cumulative sites', 'alt alleles fdr >= alt count', 'alt alleles fdr <= alt count', 'sites fdr >= alt count', 'sites fdr <= alt count', 'non-ref sensitivity <= alt count', 'non-ref discrepancy <= alt count'),
    fill=c('red', 'purple', 'blue', 'green', 'magenta', 'cyan', 'yellow', 'brown'))
garbage <- dev.off()
