#
# Author: ahfyth
#
# R script for 
#

argv = commandArgs(T);
if( length(argv) != 3 && length(argv) != 5 && length(argv) != 9 && length(argv) != 13 )
{
	print( 'usage: R --slave --args <out.pdf> <1.size> <1.Name> [<2.size> <2.Name>] < plot.R' );
	q();
}

outfileName=argv[1];
pdf( outfileName,w=6,h=4 );

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
} else if (length(argv) == 13) {	## six files, which is the most common case
	sizeA = read.table( argv[2],header=T);
	nameA = argv[3];
	sizeB = read.table( argv[4],header=T );
	nameB = argv[5];
	sizeC = read.table( argv[6],header=T );
	nameC = argv[7];
	sizeD = read.table( argv[8],header=T );
	nameD = argv[9];
	sizeE = read.table( argv[10],header=T );
	nameE = argv[11];
	sizeF = read.table( argv[12],header=T );
	nameF = argv[13];
	ymax=max( sizeA[,3],sizeB[,3],sizeC[,3],sizeD[,3],sizeE[,3],sizeF[,3]);

	colour=c('#f89588',"#9987ce","#efa666","#63b2ee","#7898e1","#0780cf")
	plot( sizeA[,3]~sizeA[,1], type='l', lwd=3, xlim=c(0, 500), ylim=c(0, ymax), col='#f89588',
			xlab="Fragment size (bp)", ylab="Frequency (%)",
			main="Size distribution comparison" );
	lines( sizeB[,3]~sizeB[,1], lwd=2, col=colour[2] );
	lines( sizeC[,3]~sizeC[,1], lwd=2, col=colour[3] );
	lines( sizeD[,3]~sizeD[,1], lwd=2, col=colour[4] );
	lines( sizeE[,3]~sizeE[,1], lwd=2, col=colour[5] );
	lines( sizeF[,3]~sizeF[,1], lwd=2, col=colour[6] );
	legend( 'topright', c(nameA, nameB, nameC, nameD, nameE,nameF), col=c('#f89588',"#9987ce","#efa666","#63b2ee","#7898e1","#0780cf"), lty=c(1,1), bty='n', cex=1.1 );

	plot( sizeA[,3]~sizeA[,1], type='l', lwd=3, xlim=c(0, 250), ylim=c(0, ymax), col='#f89588',
		 xlab="Fragment size (bp)", ylab="Frequency (%)",
		 main="Size distribution comparison" );
	lines( sizeB[,3]~sizeB[,1], lwd=2, col=colour[2] );
	lines( sizeC[,3]~sizeC[,1], lwd=2, col=colour[3] );
	lines( sizeD[,3]~sizeD[,1], lwd=2, col=colour[4] );
	lines( sizeE[,3]~sizeE[,1], lwd=2, col=colour[5] );
	lines( sizeF[,3]~sizeF[,1], lwd=2, col=colour[6] );
	legend( 'topleft', c(nameA, nameB, nameC,nameD,nameE,nameF), col=c('#f89588',"#9987ce","#efa666","#63b2ee","#7898e1","#0780cf"), lty=c(1,1), bty='n', cex=1.1 );
}

