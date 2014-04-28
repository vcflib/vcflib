#usage:  nohup R --vanilla < plotPfst --args smoothedpFst.txt

cmd_args <- commandArgs(trailingOnly = TRUE)

plotPfst<-function(x){
	require("ggplot2")
	dat<-read.table( x, header=FALSE )
	dat$V2<-1:length(dat$V2)
	theplot<-ggplot(dat, aes(x=V2, y=-log10(V5), colour=V4))+geom_point()+theme_grey(15)+labs(x="SNP index", y="-log10(smoothed pFst)")+scale_colour_continuous(low="grey", high="red", name="variants in window")
	pngName<-paste(c(x, format(Sys.time(), "%a%b%d_%H_%M_%S.png")), collapse="_")
	ggsave(filename=pngName, width=20, height=4, units="in", theplot)
}

plotPfst(cmd_args)