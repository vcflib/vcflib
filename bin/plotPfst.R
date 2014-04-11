#usage:  nohup R --vanilla < plotPfst --args pFst.txt

cmd_args <- commandArgs(trailingOnly = TRUE)

plotPfst<-function(x){
	require("ggplot2")
	dat<-read.table( x, header=FALSE )
	dat$V2<-1:length(dat$V2)
	dat<-dat[dat$V3 < 0.9,]
	theplot<-ggplot(dat, aes(x=V2, y=-log10(V3)))+geom_point()+theme_grey(15)+labs(x="SNP index", y="-log10(pFst)")
	pngName<-paste(c(x, format(Sys.time(), "%a%b%d_%H_%M_%S.png")), collapse="_")
	ggsave(filename=pngName, width=20, height=4, units="in", theplot)
}

plotPfst(cmd_args)