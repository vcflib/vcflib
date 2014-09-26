#usage:  nohup R --vanilla < plotPfst --args pFst.txt

cmd_args <- commandArgs(trailingOnly = TRUE)

plotPfst<-function(x){
	require("ggplot2")
	dat<-read.table( x, header=FALSE )
	dat$V2 <-dat$V2 / 1e3

	theplot<-ggplot(dat, aes(x=V2, y=V9))+geom_point()+geom_segment(aes(x=V2, xend=V2, y=V10, yend=V11))+labs(x="KB position", y="Fst")+theme_grey(15)
	pngName<-paste(c(x, format(Sys.time(), "%a%b%d_%H_%M_%S.png")), collapse="_")
	ggsave(filename=pngName, width=20, height=4, units="in", theplot)
}

plotPfst(cmd_args)