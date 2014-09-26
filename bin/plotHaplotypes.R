#usage:  nohup R --vanilla < plotPfst --args plotHapOutput.txt

cmd_args <- commandArgs(trailingOnly = TRUE)


imageHap<-function(x){

	pngName<-paste(c(x, format(Sys.time(), "%a%b%d_%H_%M_%S.pdf")), collapse="_")

	dat<-read.table( x[1], header=FALSE )
	pos<-dat[,1]
	dat<-dat[,2:length(dat)]
	print(head(dat))
	hd<-dist(t(dat), method="binary")
	or<-hclust(hd)
	or$labels<-1:length(dat)
	print(or$labels)
	
	pdf(pngName, width=9, height=8)
	par(mfrow=c(2,1))
	image(1-as.matrix(dat[,or$order]), yaxt="n", xaxt="n", ylab="Haplotypes", xlab="SNP", cex.lab=1.1)
	plot(or, main="", cex=1.1, lwd=2, sub="", xlab="")
	dev.off()

}

imageHap(cmd_args)