#!/usr/bin/env Rscript

# helper functions

nan.to.zero <- function(n) {
    if (is.nan(n)) return(0) else return(n)
}


# get the input VCF tabular format, assert that sites must have AC > 0
vcf <- subset(read.table(pipe('cat /dev/stdin'), header=T), AC > 0)

filename <- commandArgs(TRUE)[1]
tag <- commandArgs(TRUE)[2]

tag.genotypes_count <- paste(tag, '.genotypes.alternate_count', sep='')
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
vcf.sitesTruePositive <- mean(vcf[, tag.has_variant])

min_sites <- 5  # number of sites required for "simple plotting"

#library(ggplot2)
#vcf2 <- data.frame(QUAL=vcf$QUAL, AC=vcf$AC, has_variant=vcf[, tag.has_variant])
#qplot(AC, has_variant, group=AC, geom="boxplot", data=subset(vcf2, AC <= 20))
#ggsave(paste(filename, '.', tag, '.PD.vs.AC.boxplot.ac_lt_20.pdf', sep=''))


cat('number of sites', vcf.numberOfSites, '\n')
cat('total alternate alleles', vcf.totalAltAlleles, '\n')
cat('positive discrepancy', vcf.positiveDiscrepancy, '\n')
cat('negative discrepancy', vcf.negativeDiscrepancy, '\n')

#x <- cbind(tapply(vcf, as.list(seq(0,max(vcf$AC))),
#    function(x) {
#        sum(x[, tag.alternate_positive_discrepancy]) / sum(x[, tag.genotypes_alternate_count])
#    }))

byac <- data.frame(ac=as.vector(seq(1,max(vcf$AC)))) #, fdr=as.vector(x))


byac$fdr <- as.vector(cbind(by(byac$ac, byac$ac,
    function(x) {
        s <- subset(vcf, AC == x)
        return(nan.to.zero(sum(s[, tag.alternate_positive_discrepancy]) / sum(s[, tag.genotypes_alternate_count])))
    })))

# false detection count
byac$fpc <- as.vector(cbind(by(byac$ac, byac$ac,
    function(x) {
        s <- subset(vcf, AC == x)
        return(sum(s[, tag.alternate_positive_discrepancy]))
    })))

byac$alleles <- as.vector(cbind(by(byac$ac, byac$ac,
    function(x) {
        s <- subset(vcf, AC == x)
        return(sum(s[, tag.genotypes_alternate_count]))
    })))

byac$sites <- as.vector(cbind(by(byac$ac, byac$ac,
    function(x) {
        s <- subset(vcf, AC == x)
        return(length(s$AC))
    })))

# count true positive sites
byac$site_tpc <- as.vector(cbind(by(byac$ac, byac$ac,
    function(x) {
        s <- subset(vcf, AC == x)
        return(sum(s[, tag.has_variant]))
    })))

# fpc == false detection count
byac$site_fpc <- byac$sites - byac$site_tpc
# site detection fpr is 1 - true positive rate
byac$site_fpr <- 1 - ( byac$site_tpc / byac$sites )

summary(byac)

#print(byac$sites)
#print(byac$site_tpc)
#print(byac$site_fpc)
#print(byac$site_fpr)

#byac$site_fpr_gt0 <- as.vector(cbind(tapply(byac$ac, byac$ac, function(i) { 
#    s <- subset(byac, ac == i, select=c(site_fpr, sites))
#if (s$sites >= min_sites) {
#    return(s$site_fpr)
#} else {
#    return(NA)
#}
#})))

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

byac$alternate_pdr <- as.vector(cbind(tapply(byac$ac, byac$ac, function(i) { 
    s <- subset(byac, ac == i, select=c(fpc, alleles))
    return(sum(s$fpc) / sum(s$alleles))
})))

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

byac$nrs <- as.vector(cbind(by(byac$ac, byac$ac,
    function(x) {
        s <- subset(vcf, AC == x)
        return(nan.to.zero(sum(s[, tag.non_reference_sensitivity_count]) / sum(s[, tag.non_reference_sensitivity_normalizer])))
    })))

byac$nrd <- as.vector(cbind(by(byac$ac, byac$ac,
    function(x) {
        s <- subset(vcf, AC == x)
        return(nan.to.zero(sum(s[, tag.non_reference_discrepancy_count]) / sum(s[, tag.non_reference_discrepancy_normalizer])))
    })))

byac$nrslt <- as.vector(cbind(tapply(byac$ac, byac$ac, function(i) {
    s <- subset(vcf, AC <= i, select=c(tag.non_reference_sensitivity_count, tag.non_reference_sensitivity_normalizer))
    return(nan.to.zero(sum(s[, tag.non_reference_sensitivity_count]) / sum(s[, tag.non_reference_sensitivity_normalizer])))
})))

byac$nrdlt <- as.vector(cbind(tapply(byac$ac, byac$ac, function(i) { 
    s <- subset(vcf, AC <= i, select=c(tag.non_reference_discrepancy_count, tag.non_reference_discrepancy_normalizer))
    return(nan.to.zero(sum(s[, tag.non_reference_discrepancy_count]) / sum(s[, tag.non_reference_discrepancy_normalizer])))
})))

byac_gtsites <- subset(byac, sites >= min_sites)


if (FALSE) {
pdf(paste(filename, '.', tag, '.PD.vs.AC.smooth.pdf', sep=''))
par(cex=0.75)
par(mar=c(5,4,4,5) + 0.1)
plot(byac$cfa, ylim=c(0,1.0),
    xlab='alternate allele count (AC)', xaxt='n',
    ylab='', yaxt='n', type='l', col='red')
axis(2, at=seq(0,1,0.1), labels=seq(0,1,0.1))
axis(1, at=seq(0,max(byac$ac),10), labels=seq(0,max(byac$ac),10), cex=0.75)
grid(lty=5)
par(new=T)
title(paste(filename, 'positive discrepancy versus', tag, '(smoothed)'))
par(new=T)
countTicks <- seq(0,1,0.1) * vcf.numberOfSites
axis(4, at=seq(0,1,0.1), labels=round(countTicks))
mtext("number of sites", side=4, line=3, cex=0.75)
par(new=T)
#plot(byac$alternate_pdlt, ylim=c(0,1.0), xlab='', xaxt='n', ylab='', yaxt='n', type='l', col='green')
plot(byac$site_fpr, ylim=c(0,1.0), xlab='', xaxt='n', ylab='', yaxt='n')
par(new=T)
lines(byac$ac, predict(loess(byac$alternate_pdr ~ byac$ac, span=0.5)), col="green")
par(new=T)
plot(byac$cfa, ylim=c(0,1.0), xlab='', xaxt='n', ylab='', yaxt='n', type='l', col='red')
par(new=T)
plot(byac$cfs, ylim=c(0,1.0), xlab='', xaxt='n', ylab='', yaxt='n', type='l', col='purple')
par(new=T)
lines(byac$ac, predict(loess(byac$site_fpr ~ byac$ac, span=0.5)), col="blue")
par(new=T)
lines(byac$ac, predict(loess(byac$nrs ~ byac$ac, span=0.5)), col="magenta")
par(new=T)
lines(byac$ac, predict(loess(byac$nrd ~ byac$ac, span=0.5)), col="brown")
par(new=T, cex=0.65)
mtext(paste("alternate genotype PD: ", round(vcf.positiveDiscrepancy, digits=4), ", site PD: ", round(1 - vcf.sitesTruePositive, digits=4), sep=''))
par(new=T, cex=0.65)
legend('topleft', c('cumulative fraction of alt alleles', 'cumulative fraction of sites', 'alt genotypes PD', 'site PD', 'non-ref sensitivity', 'non-ref discrepancy', 'site PD at AC'),
    fill=c('red', 'purple', 'green', 'blue', 'magenta', 'brown', 'black'))
garbage <- dev.off()
}



pdf(paste(filename, '.', tag, '.PD.vs.AC.cumulative.pdf', sep=''))
par(cex=0.75)
par(mar=c(5,4,4,5) + 0.1)
plot(byac$cfa, ylim=c(0,1.0), 
    xlab='alternate allele count (AC)', xaxt='n',
    ylab='', yaxt='n', type='l', col='red')
axis(2, at=seq(0,1,0.1), labels=seq(0,1,0.1))
axis(1, at=seq(0,max(byac$ac),10), labels=seq(0,max(byac$ac),10), cex=0.75)
grid(lty=5)
par(new=T)
title(paste(filename, 'positive discrepancy versus', tag, '(cumulative)'))
par(new=T)
countTicks <- seq(0,1,0.1) * vcf.numberOfSites
axis(4, at=seq(0,1,0.1), labels=round(countTicks))
mtext("number of sites", side=4, line=3, cex=0.75)
par(new=T)
#plot(byac_gtsites$ac, byac_gtsites$site_fpr, ylim=c(0,1.0), xlab='', xaxt='n', ylab='', yaxt='n')
plot(byac$site_fpr, ylim=c(0,1.0),  xlab='', xaxt='n', ylab='', yaxt='n')
par(new=T)
plot(byac$alternate_pdlt, ylim=c(0,1.0),  xlab='', xaxt='n', ylab='', yaxt='n', type='l', col='green')
par(new=T)
plot(byac$cfa, ylim=c(0,1.0),  xlab='', xaxt='n', ylab='', yaxt='n', type='l', col='red')
par(new=T)
plot(byac$cfs, ylim=c(0,1.0),  xlab='', xaxt='n', ylab='', yaxt='n', type='l', col='purple')
par(new=T)
plot(byac$site_fprlt, ylim=c(0,1.0),  xlab='', xaxt='n', ylab='', yaxt='n', type='l', col='blue')
par(new=T)
plot(byac$nrslt, ylim=c(0,1.0),  xlab='', xaxt='n', ylab='', yaxt='n', type='l', col='magenta')
par(new=T)
plot(byac$nrdlt, ylim=c(0,1.0),  xlab='', xaxt='n', ylab='', yaxt='n', type='l', col='brown')
par(new=T, cex=0.65)
mtext(paste("alternate genotype PD: ", round(vcf.positiveDiscrepancy, digits=4), ", site PD: ", round(1 - vcf.sitesTruePositive, digits=4), sep=''))
par(new=T, cex=0.65)
legend('topleft', c('cumulative fraction of alt alleles', 'cumulative fraction of sites', 'alt genotypes PD', 'site PD', 'non-ref sensitivity', 'non-ref discrepancy', 'site PD at AC'),
    fill=c('red', 'purple', 'green', 'blue', 'magenta', 'brown', 'black'))
garbage <- dev.off()



pdf(paste(filename, '.', tag, '.PD.vs.AC.cumulative.simple.pdf', sep=''))
par(cex=0.75)
par(mar=c(5,4,4,5) + 0.1)
plot(byac$cfs, ylim=c(0,1.0), xlim=c(0,max(vcf$AC)),
    xlab='alternate allele count (AC)', xaxt='n',
    ylab='', yaxt='n', type='l', col='purple')
axis(2, at=seq(0,1,0.1), labels=seq(0,1,0.1))
axis(1, at=seq(0,max(byac$ac),10), labels=seq(0,max(byac$ac),10), cex=0.75)
grid(lty=5)
par(new=T)
title(paste(filename, 'positive discrepancy versus', tag, '(cumulative)'))
par(new=T)
countTicks <- seq(0,1,0.1) * vcf.numberOfSites
axis(4, at=seq(0,1,0.1), labels=round(countTicks))
mtext("number of sites", side=4, line=3, cex=0.75)
par(new=T)
plot(byac$site_fprlt, ylim=c(0,1.0), xlim=c(0,max(vcf$AC)), xlab='', xaxt='n', ylab='', yaxt='n', type='l', col='blue')
par(new=T)
plot(byac_gtsites$ac, byac_gtsites$site_fpr, ylim=c(0,1.0), xlim=c(0,max(vcf$AC)),   xlab='', xaxt='n', ylab='', yaxt='n')
par(new=T, cex=0.65)
mtext(paste("site PD: ", round(1 - vcf.sitesTruePositive, digits=4), sep=''))
par(new=T, cex=0.65)
legend('topleft', c('cumulative fraction of sites', 'cumulative site PD', paste('site PD at AC (>=', min_sites, 'sites)')),
    fill=c('purple', 'blue', 'black'))
garbage <- dev.off()




pdf(paste(filename, '.', tag, '.PD.vs.AC.instantaneous.ac_lt_20.pdf', sep=''))
#par(cex=0.75)
par(mar=c(5,4,4,5) + 0.1)
plot(byac$sites, ylim=c(0,max(byac$sites)), xlim=c(0,20),
    xlab='alternate allele count (AC)', xaxt='n',
    ylab='number of sites', type='l', pch=19, col='blue')
#axis(2, at=seq(0,1,0.1), labels=seq(0,1,0.1))
#countTicks <- round(seq(0,1,0.1) * max(byac$sites))
#axis(2, at=countTicks, labels=countTicks)
par(new=T)
axis(1, at=seq(0,max(byac$ac),1), labels=seq(0,max(byac$ac),1))
grid(lty=5)
par(new=T)
plot(byac$sites, ylim=c(0,max(byac$sites)), xlim=c(0,20), type='o', pch=19, col='blue', xlab='', xaxt='n', ylab='', yaxt='n')
par(new=T)
title(paste(filename, 'positive discrepancy versus', tag, '(instantaneous)'))
par(new=T)
mtext("number of sites", side=2, line=3) #, cex=0.75)
#par(new=T)
#plot(byac$site_fprlt, ylim=c(0,1.0), xlim=c(0,20),  xlab='', xaxt='n', ylab='', yaxt='n', type='l', col='blue')
par(new=T)
plot(byac_gtsites$ac, byac_gtsites$site_tpc, ylim=c(0,max(byac$sites)), xlim=c(0,20), xlab='', xaxt='n', ylab='', yaxt='n', col='red', pch=19, type='o')
par(new=T) #, cex=0.65)
mtext(paste("site PD: ", round(1 - vcf.sitesTruePositive, digits=4), sep=''))
par(new=T) #, cex=0.65)
#legend('topright', c('number of sites', 'site PD count'),
#    fill=c('blue', 'red'))
garbage <- dev.off()

# stratifying by QUAL

if (FALSE) {


x <- cbind(by(vcf, vcf$QUAL,
    function(x) {
        sum(x[, tag.alternate_positive_discrepancy]) / sum(x[, tag.genotypes_alternate_count])
    }))

byqual <- data.frame(qual=as.numeric(rownames(x)), fdr=as.vector(x))

# false detection count
byqual$fpc <- as.vector(cbind(by(vcf, vcf$QUAL,
    function(i) { sum(i[, tag.alternate_positive_discrepancy]) } )))

byqual$alleles <- as.vector(cbind(by(vcf, vcf$QUAL,
    function(i) {
        sum(i[, tag.genotypes_alternate_count])
    })))

byqual$sites <- as.vector(cbind(by(vcf$QUAL, vcf$QUAL, function(i) length(i))))

# count true positive sites
byqual$site_tpc <- as.vector(cbind(by(vcf[, tag.has_variant], vcf$QUAL, function(i) sum(i))))
# fpc == false detection count
byqual$site_fpc <- byqual$sites - byqual$site_tpc
# site detection fpr is 1 - true positive rate
byqual$site_fpr <- 1 - ( byqual$site_tpc / byqual$sites )

#byqual$site_fpr_gt0 <- as.vector(cbind(tapply(byqual$ac, byqual$ac, function(i) { 
#    s <- subset(byqual, ac == i, select=c(site_fpr, sites))
#if (s$sites >= min_sites) {
#    return(s$site_fpr)
#} else {
#    return(NA)
#}
#})))

#byqual$site_fprlt <- as.vector(cbind(tapply(byqual$qual, byqual$qual, function(i) mean(subset(byqual, qual <= i, select=site_fpr)))))
byqual$site_fprlt <- as.vector(cbind(tapply(byqual$qual, byqual$qual, function(i) { 
    s <- subset(byqual, qual <= i, select=c(site_fpc, sites))
    return(sum(s$site_fpc) / sum(s$sites))
})))

#byqual$site_fprgt <- as.vector(cbind(tapply(byqual$qual, byqual$qual, function(i) mean(subset(byqual, qual >= i, select=site_fpr)))))
byqual$site_fprgt <- as.vector(cbind(tapply(byqual$qual, byqual$qual, function(i) { 
    s <- subset(byqual, qual >= i, select=c(site_fpc, sites))
    return(sum(s$site_fpc) / sum(s$sites))
})))

byqual$cfa <- as.vector(cbind(tapply(byqual$qual, byqual$qual, function(i) sum(subset(byqual, qual <= i, select=alleles)) / sum(byqual$alleles))))

byqual$cfs <- as.vector(cbind(tapply(byqual$qual, byqual$qual, function(i) sum(subset(byqual, qual <= i, select=sites)) / length(vcf$QUAL))))

# inappropriate collapse via averaging of fdr
#byqual$alternate_pdlt <- as.vector(cbind(tapply(byqual$qual, byqual$qual, function(i) mean(subset(byqual, qual <= i, select=fdr)))))

byqual$alternate_pdr <- as.vector(cbind(tapply(byqual$qual, byqual$qual, function(i) { 
    s <- subset(byqual, qual == i, select=c(fpc, alleles))
    return(sum(s$fpc) / sum(s$alleles))
})))

# use this one
byqual$alternate_pdlt <- as.vector(cbind(tapply(byqual$qual, byqual$qual, function(i) { 
    s <- subset(byqual, qual <= i, select=c(fpc, alleles))
    return(sum(s$fpc) / sum(s$alleles))
})))

#byqual$alternate_pdgt <- as.vector(cbind(tapply(byqual$qual, byqual$qual, function(i) mean(subset(byqual, qual >= i, select=fdr)))))
byqual$alternate_pdgt <- as.vector(cbind(tapply(byqual$qual, byqual$qual, function(i) { 
    s <- subset(byqual, qual >= i, select=c(fpc, alleles))
    return(sum(s$fpc) / sum(s$alleles))
})))

nan.to.zero <- function(n) {
    if (is.nan(n)) return(0) else return(n)
}

byqual$nrs <- as.vector(cbind(by(vcf, vcf$QUAL, function(i) {
    return(nan.to.zero(sum(i[, tag.non_reference_sensitivity_count]) / sum(i[, tag.non_reference_sensitivity_normalizer])))
})))

byqual$nrd <- as.vector(cbind(by(vcf, vcf$QUAL, function(i) { 
    return(nan.to.zero(sum(i[, tag.non_reference_discrepancy_count]) / sum(i[, tag.non_reference_discrepancy_normalizer])))
})))

byqual$nrslt <- as.vector(cbind(tapply(byqual$qual, byqual$qual, function(i) {
    s <- subset(vcf, QUAL <= i, select=c(tag.non_reference_sensitivity_count, tag.non_reference_sensitivity_normalizer))
    return(nan.to.zero(sum(s[, tag.non_reference_sensitivity_count]) / sum(s[, tag.non_reference_sensitivity_normalizer])))
})))

byqual$nrdlt <- as.vector(cbind(tapply(byqual$qual, byqual$qual, function(i) { 
    s <- subset(vcf, QUAL <= i, select=c(tag.non_reference_discrepancy_count, tag.non_reference_discrepancy_normalizer))
    return(nan.to.zero(sum(s[, tag.non_reference_discrepancy_count]) / sum(s[, tag.non_reference_discrepancy_normalizer])))
})))

byqual_gt10 <- subset(byqual, sites >= min_sites)


if (FALSE) {
pdf(paste(filename, '.', tag, '.PD.vs.QUAL.smooth.pdf', sep=''))
par(cex=0.75)
par(mar=c(5,4,4,5) + 0.1)
plot(byqual$cfa, ylim=c(0,1.0),
    xlab='QUAL', xaxt='n',
    ylab='', yaxt='n', type='l', col='red')
axis(2, at=seq(0,1,0.1), labels=seq(0,1,0.1))
axis(1, at=seq(0,max(byqual$qual),10), labels=seq(0,max(byqual$qual),10), cex=0.75)
grid(lty=5)
par(new=T)
title(paste(filename, 'positive discrepancy versus', tag, '(smoothed)'))
par(new=T)
countTicks <- seq(0,1,0.1) * vcf.numberOfSites
axis(4, at=seq(0,1,0.1), labels=round(countTicks))
mtext("number of sites", side=4, line=3, cex=0.75)
par(new=T)
#plot(byqual$alternate_pdlt, ylim=c(0,1.0), xlab='', xaxt='n', ylab='', yaxt='n', type='l', col='green')
plot(byqual$site_fpr, ylim=c(0,1.0), xlab='', xaxt='n', ylab='', yaxt='n')
par(new=T)
lines(byqual$qual, predict(loess(byqual$alternate_pdr ~ byqual$qual, span=0.5)), col="green")
par(new=T)
plot(byqual$cfa, ylim=c(0,1.0), xlab='', xaxt='n', ylab='', yaxt='n', type='l', col='red')
par(new=T)
plot(byqual$cfs, ylim=c(0,1.0), xlab='', xaxt='n', ylab='', yaxt='n', type='l', col='purple')
par(new=T)
lines(byqual$qual, predict(loess(byqual$site_fpr ~ byqual$qual, span=0.5)), col="blue")
par(new=T)
lines(byqual$qual, predict(loess(byqual$nrs ~ byqual$qual, span=0.5)), col="magenta")
par(new=T)
lines(byqual$qual, predict(loess(byqual$nrd ~ byqual$qual, span=0.5)), col="brown")
par(new=T, cex=0.65)
mtext(paste("alternate genotype PD: ", round(vcf.positiveDiscrepancy, digits=4), ", site PD: ", round(1 - vcf.sitesTruePositive, digits=4), sep=''))
par(new=T, cex=0.65)
legend('topleft', c('cumulative fraction of alt alleles', 'cumulative fraction of sites', 'alt genotypes PD', 'site PD', 'non-ref sensitivity', 'non-ref discrepancy', 'site PD at QUAL'),
    fill=c('red', 'purple', 'green', 'blue', 'magenta', 'brown', 'black'))
garbage <- dev.off()
}



pdf(paste(filename, '.', tag, '.PD.vs.QUAL.cumulative.pdf', sep=''))
par(cex=0.75)
par(mar=c(5,4,4,5) + 0.1)
plot(byqual$cfa, ylim=c(0,1.0),
    xlab='QUAL', xaxt='n',
    ylab='', yaxt='n', type='l', col='red')
axis(2, at=seq(0,1,0.1), labels=seq(0,1,0.1))
axis(1, at=seq(0,max(byqual$qual),10), labels=seq(0,max(byqual$qual),10), cex=0.75)
grid(lty=5)
par(new=T)
title(paste(filename, 'positive discrepancy versus', tag, '(cumulative)'))
par(new=T)
countTicks <- seq(0,1,0.1) * vcf.numberOfSites
axis(4, at=seq(0,1,0.1), labels=round(countTicks))
mtext("number of sites", side=4, line=3, cex=0.75)
par(new=T)
plot(byqual_gt10$qual, byqual_gt10$site_fpr, ylim=c(0,1.0), xlab='', xaxt='n', ylab='', yaxt='n')
#plot(byqual$site_fpr, ylim=c(0,1.0), xlab='', xaxt='n', ylab='', yaxt='n')
par(new=T)
plot(byqual$alternate_pdlt, ylim=c(0,1.0), xlab='', xaxt='n', ylab='', yaxt='n', type='l', col='green')
par(new=T)
plot(byqual$cfa, ylim=c(0,1.0), xlab='', xaxt='n', ylab='', yaxt='n', type='l', col='red')
par(new=T)
plot(byqual$cfs, ylim=c(0,1.0), xlab='', xaxt='n', ylab='', yaxt='n', type='l', col='purple')
par(new=T)
plot(byqual$site_fprlt, ylim=c(0,1.0), xlab='', xaxt='n', ylab='', yaxt='n', type='l', col='blue')
par(new=T)
plot(byqual$nrslt, ylim=c(0,1.0), xlab='', xaxt='n', ylab='', yaxt='n', type='l', col='magenta')
par(new=T)
plot(byqual$nrdlt, ylim=c(0,1.0), xlab='', xaxt='n', ylab='', yaxt='n', type='l', col='brown')
par(new=T, cex=0.65)
mtext(paste("alternate genotype PD: ", round(vcf.positiveDiscrepancy, digits=4), ", site PD: ", round(1 - vcf.sitesTruePositive, digits=4), sep=''))
par(new=T, cex=0.65)
legend('topleft', c('cumulative fraction of alt alleles', 'cumulative fraction of sites', 'alt genotypes PD', 'site PD', 'non-ref sensitivity', 'non-ref discrepancy', 'site PD at QUAL (>= 10 sites)'),
    fill=c('red', 'purple', 'green', 'blue', 'magenta', 'brown', 'black'))
garbage <- dev.off()



pdf(paste(filename, '.', tag, '.PD.vs.QUAL.cumulative.simple.pdf', sep=''))
par(cex=0.75)
par(mar=c(5,4,4,5) + 0.1)
plot(byqual$cfs, ylim=c(0,1.0),
    xlab='QUAL', xaxt='n',
    ylab='', yaxt='n', type='l', col='purple')
axis(2, at=seq(0,1,0.1), labels=seq(0,1,0.1))
axis(1, at=seq(0,max(byqual$qual),10), labels=seq(0,max(byqual$qual),10), cex=0.75)
grid(lty=5)
par(new=T)
title(paste(filename, 'positive discrepancy versus', tag, '(cumulative)'))
par(new=T)
countTicks <- seq(0,1,0.1) * vcf.numberOfSites
axis(4, at=seq(0,1,0.1), labels=round(countTicks))
mtext("number of sites", side=4, line=3, cex=0.75)
par(new=T)
plot(byqual$site_fprlt, ylim=c(0,1.0), xlab='', xaxt='n', ylab='', yaxt='n', type='l', col='blue')
par(new=T)
plot(byqual_gt10$qual, byqual_gt10$site_fpr, ylim=c(0,1.0), xlab='', xaxt='n', ylab='', yaxt='n')
par(new=T, cex=0.65)
mtext(paste("site PD: ", round(1 - vcf.sitesTruePositive, digits=4), sep=''))
par(new=T, cex=0.65)
legend('topleft', c('cumulative fraction of sites', 'site PD', 'site PD at QUAL (>= 10 sites)'),
    fill=c('purple', 'blue', 'black'))
garbage <- dev.off()

}
