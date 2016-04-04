#!/usr/bin/env Rscript --vanilla --slave

# get the input VCF tabular format, assert that sites must have AC > 0
vcf <- subset(read.table(pipe('cat /dev/stdin'), header=T), AC > 0)

filename <- commandArgs(TRUE)[1]
tag <- commandArgs(TRUE)[2]

tag.has_variant <- paste(tag, '.has_variant', sep='')

vcf.numberOfSites <- length(vcf$AC)
vcf.sitesTruePositive <- mean(vcf[, tag.has_variant])

# false detection count
x <- cbind(by(vcf$AC, vcf$AC, function(i) length(i)))

byac <- data.frame(ac=as.numeric(rownames(x)), sites=as.vector(x))

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

byac$cfs <- as.vector(cbind(tapply(byac$ac, byac$ac, function(i) sum(subset(byac, ac <= i, select=sites)) / length(vcf$AC))))


pdf(paste(filename, '.', tag, '.site_FDR.vs.AC.smooth.pdf', sep=''))
par(cex=0.75)
par(mar=c(5,4,4,5) + 0.1)
plot(byac$cfs, ylim=c(0,1.0),
    xlab='alternate allele count (AC)', xaxt='n',
    ylab='false discovery rate (FDR)', yaxt='n', type='l', col='red')
axis(2, at=seq(0,1,0.1), labels=seq(0,1,0.1))
axis(1, at=seq(0,max(byac$ac),10), labels=seq(0,max(byac$ac),10), cex=0.75)
grid(lty=5)
par(new=T)
title(paste(filename, 'putative site false discovery rate versus', tag, '(smoothed)'))
par(new=T)
countTicks <- seq(0,1,0.1) * vcf.numberOfSites
axis(4, at=seq(0,1,0.1), labels=round(countTicks))
par(col='red')
mtext("number of sites", side=4, line=3, cex=0.75)
par(col='black')
par(new=T)
plot(byac$site_fpr, ylim=c(0,1.0), xlab='', xaxt='n', ylab='', yaxt='n')
par(new=T)
lines(byac$ac, predict(loess(byac$site_fpr ~ byac$ac, span=0.5)), col="blue")
par(new=T, cex=0.65)
mtext(paste("site FDR: ", round(1 - vcf.sitesTruePositive, digits=4), sep=''))
par(new=T, cex=0.65)
legend('topleft', c('cumulative sites', 'site FDR (loess smoothed)', 'FDR at AC'),
    fill=c('red', 'blue', 'black'))
garbage <- dev.off()



pdf(paste(filename, '.', tag, '.site_FDR.vs.AC.cumulative.pdf', sep=''))
par(cex=0.75)
par(mar=c(5,4,4,5) + 0.1)
plot(byac$cfs, ylim=c(0,1.0),
    xlab='alternate allele count (AC)', xaxt='n',
    ylab='false discovery rate (FDR)', yaxt='n', type='l', col='red')
axis(2, at=seq(0,1,0.1), labels=seq(0,1,0.1))
axis(1, at=seq(0,max(byac$ac),10), labels=seq(0,max(byac$ac),10), cex=0.75)
grid(lty=5)
par(new=T)
title(paste(filename, 'putative false discovery rate versus', tag, '(cumulative)'))
par(new=T)
countTicks <- seq(0,1,0.1) * vcf.numberOfSites
axis(4, at=seq(0,1,0.1), labels=round(countTicks))
par(col='red')
mtext("number of sites", side=4, line=3, cex=0.75)
par(col='black')
par(new=T)
plot(byac$site_fpr, ylim=c(0,1.0), xlab='', xaxt='n', ylab='', yaxt='n')
par(new=T)
plot(byac$site_fprgt, ylim=c(0,1.0), xlab='', xaxt='n', ylab='', yaxt='n', type='l', col='green')
par(new=T)
plot(byac$site_fprlt, ylim=c(0,1.0), xlab='', xaxt='n', ylab='', yaxt='n', type='l', col='blue')
par(new=T, cex=0.65)
mtext(paste("site FDR: ", round(1 - vcf.sitesTruePositive, digits=4), sep=''))
par(new=T, cex=0.65)
legend('topleft', c('cumulative sites', 'site FDR <= AC', 'site FDR >= AC', 'FDR at AC'),
    fill=c('red', 'blue', 'green', 'black'))
garbage <- dev.off()
