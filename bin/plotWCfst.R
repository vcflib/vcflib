#usage:  nohup R --vanilla < plotPfst --args pFst.txt

cmd_args <- commandArgs(trailingOnly = TRUE)

plotPfst<-function(x){
	require("ggplot2")
	dat<-read.table( x, header=FALSE )
	dat$V2<-1:length(dat$V2)
	theplot<-ggplot(dat, aes(x=V2, y=V5))+geom_point()+theme_grey(15)+labs(x="SNP index", y="wcFst")+ylim(0,1.1)
	pngName<-paste(c(x, format(Sys.time(), "%a%b%d_%H_%M_%S.png")), collapse="_")
	ggsave(filename=pngName, width=20, height=4, units="in", theplot)
}

plotPfst(cmd_args)