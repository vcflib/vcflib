#usage:  nohup R --vanilla < plotPfst --args pFst.txt

cmd_args <- commandArgs(trailingOnly = TRUE)

plotPfst<-function(x){
	require("ggplot2")
	dat<-read.table( x, header=FALSE )
	dat<-dat[dat$V5 < 0.9,]
	theplot<-ggplot(dat, aes(x=V2/1e3, y=-log10(V5)*V6))+geom_point()+theme_grey(15)+labs(x="KB position", y="-log10(hapLRT * sign)")+geom_hline(aes(yintercept=0), linetype="dashed", colour="red")
	pngName<-paste(c(x, format(Sys.time(), "%a%b%d_%H_%M_%S.png")), collapse="_")
	ggsave(filename=pngName, width=20, height=4, units="in", theplot)
}

plotPfst(cmd_args)