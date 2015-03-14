#!/usr/bin/env Rscript




require(plyr)
require(ggplot2)
require(pracma)
require(grid)

argv <- commandArgs(trailingOnly = TRUE)

prefix <- gsub("\\s","", argv[1])
print(prefix)
truthset <- argv[2]
print(truthset)
results <- argv[3]
print(results)
xmin <- as.numeric(argv[4])
xmax <- as.numeric(argv[5])
ymin <- as.numeric(argv[6])
ymax <- as.numeric(argv[7])

roc <- read.delim(results)

bests <- ddply(roc, .(set), function(x) { data.frame(best_snps=with(x, min(false_negative_snps + false_positive_snps)), best_snp_threshold=min(subset(x, (false_negative_snps + false_positive_snps) == with(x, min(false_negative_snps + false_positive_snps)))$threshold ), best_indels=with(x, min(false_negative_indels + false_positive_indels)), best_indel_threshold=min(subset(x, (false_negative_indels + false_positive_indels) == with(x, min(false_negative_indels + false_positive_indels)))$threshold )) })

write.table(bests, paste(prefix, ".bests.tsv", sep=""), row.names=FALSE, quote=FALSE, sep="\t")

#abs(trapz(c(1, roc$complexfpr), c(1, roc$complextpr)))

true_snps <- with(subset(roc, set==truthset), max(num_snps))
true_indels <- with(subset(roc, set==truthset), max(num_indels))

# get ROC AUC
auc <- ddply(roc, .(set),
      function(x) {
        data.frame(
                   snp_auc=ifelse(true_snps>0,
                     with(x,
                          abs(trapz(c(1,
                                      false_positive_snps/(false_positive_snps+ max(false_negative_snps + num_snps - false_positive_snps))),
                                    c(max(1- false_negative_snps/true_snps),
                                      1- false_negative_snps/true_snps)))),
                     0),
                   indel_auc=ifelse(true_indels>0,
                     with(x,
                          abs(trapz(c(1,
                                      false_positive_indels/(false_positive_indels+ max(false_negative_indels + num_indels - false_positive_indels))),
                                    c(max(1- false_negative_indels/true_indels),
                                      1- false_negative_indels/true_indels)))),
                     0)
                   )
      }
      )

write.table(auc, paste(prefix, ".auc.tsv", sep=""), row.names=FALSE, quote=FALSE, sep="\t")


rocsnps <- ddply(roc, .(set),
      function(x) {
        data.frame(
                   FPR=
                     with(x,
                          c(1,
                            false_positive_snps/(false_positive_snps+ max(false_negative_snps + num_snps - false_positive_snps)))),
                   TPR=
                      with(x,
                          c(max(1- false_negative_snps/true_snps),
                            1- false_negative_snps/true_snps)),
                   type=as.factor("snps")
                   )
          }
      )

rocindels <- ddply(roc, .(set),
      function(x) {
        data.frame(
                   FPR=
                     with(x,
                          c(1,
                            false_positive_indels/(false_positive_indels+ max(false_negative_indels + num_indels - false_positive_indels)))),
                   TPR=
                      with(x,
                          c(max(1- false_negative_indels/true_indels),
                            1- false_negative_indels/true_indels)),
                   type=as.factor("indels")
                   )
          }
      )


if (FALSE) {
if (true_snps>0) {
  ggplot(subset(roc, set != truthset),
         aes(false_positive_snps/(false_positive_snps+with(subset(roc, set==set), max(false_negative_snps + num_snps - false_positive_snps))),
             1- false_negative_snps/with(subset(roc, set==set), max(false_negative_snps + num_snps - false_positive_snps)),
             group=set,
             color=set)) + scale_x_continuous("false positive rate") + scale_y_continuous("true positive rate") + geom_path() + theme_bw()
            + coord_cartesian(xlim=c(xmin,xmax), ylim=c(ymin,ymax))
  ggsave(paste(prefix, ".snps.png", sep=""), height=6, width=9)
}

if (true_indels>0) {
  ggplot(subset(roc, set != truthset),
         aes(false_positive_indels/(false_positive_indels+with(subset(roc, set==set), max(false_negative_indels + num_indels - false_positive_indels))),
             1- false_negative_indels/with(subset(roc, set==set), max(false_negative_indels + num_indels - false_positive_indels)),
             group=set,
             color=set)) + scale_x_continuous("false positive rate") + scale_y_continuous("true positive rate") + geom_path() + theme_bw()
            + coord_cartesian(xlim=c(xmin,xmax), ylim=c(ymin,ymax))
  ggsave(paste(prefix, ".indels.png", sep=""), height=6, width=9)
}
}


# new versions
if (true_snps>0) {
  ggplot(subset(rocsnps, set != truthset),
         aes(FPR,
             TPR,
             group=set,
             color=set)) + scale_x_continuous("false positive rate") + scale_y_continuous("true positive rate") + geom_path() + theme_bw() + coord_cartesian(xlim=c(xmin,xmax), ylim=c(ymin,ymax))
  ggsave(paste(prefix, ".snps.png", sep=""), height=6, width=9)
}

if (true_indels>0) {
  ggplot(subset(rocindels, set != truthset),
         aes(FPR,
             TPR,
             group=set,
             color=set)) + scale_x_continuous("false positive rate") + scale_y_continuous("true positive rate") + geom_path() + theme_bw() + coord_cartesian(xlim=c(xmin,xmax), ylim=c(ymin,ymax))
  ggsave(paste(prefix, ".indels.png", sep=""), height=6, width=9)
}

if (true_indels>0 && true_snps>0) {

(
  ggplot(subset(rbind(rocsnps,rocindels), set != truthset),
         aes(FPR,
             TPR,
             group=set,
             color=set))
    + scale_x_continuous("false positive rate")
    + scale_y_continuous("true positive rate")
    + geom_path()
    + theme_bw()
    + coord_cartesian(xlim=c(xmin,xmax), ylim=c(ymin,ymax))
    + facet_grid(type ~ .)
    + theme(panel.margin = unit(1, "lines")) 
)
  ggsave(paste(prefix, ".both.png", sep=""), height=5, width=5)

}
