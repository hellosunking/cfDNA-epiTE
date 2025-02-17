#
# Author: ahfyth
#
# R script for 
#
library(scales)
library(ggplot2)
library(patchwork)
argv = commandArgs(T);
if( length(argv) != 3 && length(argv) != 5 && length(argv) != 11 )
{
	print( 'usage: R --slave --args <out.pdf> <1.size> <1.Name> [<2.size> <2.Name>] < plot.R' );
	q();
}

outfileName=argv[1];
pdf( outfileName , width=10, height=4);

if( length(argv) == 3 ) {	## only one file
	b2r_colors = colorRampPalette(c("blue", "red"))(10);
	rand_color = b2r_colors[ sample.int(10,1) ];

	dat = read.table( argv[2] );
	plot( dat[,3]~dat[,1], type='l', lwd=2, xlim=c(0, 250), col=rand_color,
			xlab="Fragment size (bp)", ylab="Frequency (%)", cex.lab=1.2,
			main=paste("Size distribution for sample ", argv[3], sep=""), col.main=rand_color );
	plot( dat[,4]~dat[,1], type='l', lwd=2, xlim=c(0, 250), ylim=c(0, 1), col=rand_color,
			xlab="Fragment size (bp)", ylab="Cumulative frequency (%)", cex.lab=1.2,
			main=paste("Size distribution for sample ", argv[3], sep=""), col.main=rand_color );
	plot( dat[,3]~dat[,1], type='l', lwd=2, xlim=c(0, 600), col=rand_color,
			xlab="Fragment size (bp)", ylab="Frequency (%)", cex.lab=1.2,
			main=paste("Size distribution for sample ", argv[3], sep=""), col.main=rand_color );
	plot( dat[,4]~dat[,1], type='l', lwd=2, xlim=c(0, 600), ylim=c(0, 1), col=rand_color,
			xlab="Fragment size (bp)", ylab="Cumulative frequency (%)", cex.lab=1.2,
			main=paste("Size distribution for sample ", argv[3], sep=""), col.main=rand_color );
} else {	## five  files, which is the most common case
	sizeA = read.table( argv[2] );
	nameA = argv[3];
	sizeB = read.table( argv[4] );
	nameB = argv[5];
	sizeC = read.table( argv[6] );
	nameC = argv[7];
	sizeD = read.table( argv[8] );
	nameD = argv[9];
	sizeE = read.table( argv[10] );
	nameE = argv[11];

	ymax=max( sizeA[,3], sizeB[,3],sizeC[,3],sizeD[,3],sizeE[,3] );
	colour=rainbow(4,alpha = 0.4)

## 0~250 bp, for most cases
	df1 <- data.frame(len=sizeA[,1],A=sizeA[,3])
	df2 <- data.frame(len=sizeB[,1],B=sizeB[,3])
	df3 <- data.frame(len=sizeC[,1],C=sizeC[,3])
	df4 <- data.frame(len=sizeD[,1],D=sizeD[,3])
	df5 <- data.frame(len=sizeE[,1],E=sizeE[,3])
	df_list <- list(df1, df2, df3, df4, df5)  

	#data <- data.frame(cbind(sizeA[,1],sizeA[,3]),sizeB[,3])
	data<-Reduce(function(x, y) merge(x, y, all= TRUE ), df_list) 
	data[is.na(data)]=0

	p1<-ggplot(data,aes(x=data[,1]))+
	geom_line(aes(y = data[,2], color = nameA))+
	geom_line(aes(y = data[,3], color = nameB))+
	geom_line(aes(y = data[,4], color = nameC))+
	geom_line(aes(y = data[,5], color = nameD))+
	geom_line(aes(y = data[,6], color = nameE))+
	scale_color_manual(NULL,values=c('black',"#fa8080","#ffc076","#73abf5","#cb9bff"))+
	#scale_fill_discrete(breaks=c(nameA,nameB,nameC,nameD,nameE))+
	coord_cartesian(xlim=c(0, 500))+
	labs(title = 'Size distribution comparison', x = 'Fragment size (bp)', y = 'Frequency (%)')+
	theme(plot.title = element_text(hjust = 0.5)) +
	theme_bw(base_size=15)+theme(panel.grid = element_blank(),
								 axis.line = element_line(colour = "#000000",size=0.15),
								 axis.text = element_text(colour = "#000000"),
								 legend.position = c(0.9, 0.8),
								 legend.background = element_blank(),
								 legend.key = element_blank(),
								 plot.title = element_text(size=12,hjust=0.5));
	p2<-ggplot(data,aes(x=data[,1]))+
	geom_line(aes(y = data[,2], color = nameA))+
	geom_line(aes(y = data[,3], color = nameB))+
	geom_line(aes(y = data[,4], color = nameC))+
	geom_line(aes(y = data[,5], color = nameD))+
	geom_line(aes(y = data[,6], color = nameE))+
	scale_color_manual(NULL, values=c('black',"#fa8080","#ffc076","#73abf5","#cb9bff"))+
	coord_cartesian(xlim=c(0, 500))+
	labs(title = 'Size distribution comparison', x = 'Fragment size (bp)', y = 'Frequency (%) (logarithmic scale)')+
	theme(plot.title = element_text(hjust = 0.5)) +
		   #type='l', lwd=2, xlim=c(0, 250),col='red')
	#		xlab="Fragment size (bp)", ylab="Frequency (%)",
	#		main="Size distribution comparison" )
		    scale_y_log10(
			   breaks = scales::trans_breaks("log10", function(x) 10^x),
			      labels = scales::trans_format("log10", scales::math_format(10^.x))
			 ) + theme_bw(base_size=15)+theme(panel.grid = element_blank(),
			 axis.line = element_line(colour = "#000000",size=0.15),
			 axis.text = element_text(colour = "#000000"),
			 legend.position = c(0.9, 0.8),
			 legend.background = element_blank(),
			 legend.key = element_blank(),
			 plot.title = element_text(size=12,hjust=0.5)) ;
	print(p1)
	print(p2)
}

