# Programs and scripts used in Gong, Pan, and Lin et al. manuscript
Distributed under the [CC BY-NC-ND 4.0](https://creativecommons.org/licenses/by-nc-nd/4.0/ "CC BY-NC-ND")  license for **personal and academic usage only.**

The following software and packages are required:
- [ktrim v1.5.0](https://github.com/hellosunking/Ktrim/releases/tag/v1.5.0)
- [bowtie2 v2.3.5.1](https://github.com/BenLangmead/bowtie2/releases/tag/v2.3.5.1)
- [Msuite2 v2.1.0](https://github.com/hellosunking/Msuite2/releases/tag/v2.1.0)
- [bedtools v2.29.2](https://github.com/arq5x/bedtools2/releases/tag/v2.29.2)
- [R v4.2.0](https://cran.r-project.org/bin/windows/base/old/4.2.0), requires "ggplot2", "ggpubr", "ggprism", "reshape2", "forcats", "ggbeeswarm", "dplyr", "cowplot", "tidyverse", "verification", "gbm", "caret", "foreach", "doParallel", "binom", and "pROC" packages.

## 1. Prepare Transposon Element (TE) regions
```
## download TE regions annotated by repeatmasker from UCSC
wget -O hg38.rmsk.txt.gz  https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz
wget -O mm10.rmsk.txt.gz https://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/rmsk.txt.gz
wget -O canFam6.rmsk.txt.gz https://hgdownload.soe.ucsc.edu/goldenPath/canFam6/database/rmsk.txt.gz

## If you do not have genome data, download the files from UCSC
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz -O hg38.fa.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz -O mm10.fa.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/canFam6/bigZips/canFam6.fa.gz -O canFam6.fa.gz
## You also need to decompress the files
gzip -d hg38.fa.gz mm10.fa.gz canFam6.fa.gz

## extract regions for each family of TE
## we had included the selected TE list (with more than 10000 copies) for human, mouse, and dog.
## the compiled regions for hg38 are also included in this package
for genome in hg38 mm10 canFam6
do
    mkdir -p $genome.TE
    while read class family extra
    do
        zcat $genome.rmsk.txt.gz | perl -lane 'print "$F[5]\t$F[6]\t$F[7]" if $F[5]=~/^chr\d+$/ && $F[11] eq "'$class'" && $F[12] eq "'$family'"' | sort -k1,1 -k2,2n | gzip > $genome.TE/$class.$family.bed.gz
    done < 1.prepare.bed/$genome.TE.list

    ## generate TE.info
    cd $genome.TE
    python ../1.prepare.bed/prepare.TE.info.py TE.info *.bed.gz
    cd ../
done
```

## 2. Read alignment
For whole genome sequencing, cfChIP-seq, and ATAC-seq data:
```
## need to update sampleID and path to the FASTQ files
PRG=2.fragmentomics/0.alignment
sid=SampleID
FASTQ1=/path/to/$sid.R1.fq.gz
FASTQ2=/path/to/$sid.R2.fq.gz

## data preprocess. Need to update the sequencing kits in "-k" option
ktrim -1 $FASTQ1 -2 $FASTQ2 -t 8 -p 33 -o $sid.ktrim -k illumina > $sid.ktrim.log

## read alignment. Need to build index for bowtie2 first
hg38index=/path/to/hg38.bowtie2.index
bowtie2 --score-min L,0,-0.2 --ignore-quals --no-unal --no-head -p 32 --minins 0 --maxins 1000 --no-mixed --no-discordant -x $hg38index -1 $sid.ktrim.read1.fq -2 $sid.ktrim.read2.fq -S $sid.sam 2> $sid.bowtie2.log

## remove dupliate and convert to BED format, keep autosomal reads only
chrsize=2.fragmentomics/0.alignment/hg38.chr.size
$PRG/ksam_rmdup $chrsize PE $sid.sam $sid.rmdup | perl $PRG/pe_sam2bed.pl /dev/stdin /dev/stdout $sid.size 0 1 2>$sid.chr.count | sort -k1,1 -k2,2n -k3,3n | gzip > $sid.bed.gz
## The $sid.bed.gz file will be used in the following analyses.
```

For EM-seq data:
```
## Need to build hg38 index for Msuite2 first. The data is sequenced on illumina sequencers
sid=sampleID
FASTQ1=/path/to/$sid.R1.fq.gz
FASTQ2=/path/to/$sid.R2.fq.gz

PRG=2.fragmentomics/2.methylation/
genome=hg38

msuite2 -x $genome -1 $FASTQ1 -2 $FASTQ2 -o Msuite2.$sid -k illumina --cut-r1-tail 25 --cut-r2-head 25 --aligner hisat2 -p 16
cd Msuite2.$sid
make
make clean
cd ..

## convert alignment result to BED format
perl $PRG/Msuite2.bam2bed.pl Msuite2.$sid/Msuite2.final.bam $sid.bed.gz 30 0 25

## calculate DNA methylation levels for each TE copy, for TEANA-Dx model
while read TEclass TEfamily extra
do
    zcat $genome.TE/$TEclass.$TEfamily.bed.gz | $PRG/extract.meth.in.region $chrsize Msuite2.$sid/Msuite2.CpG.meth.call /dev/stdin > $sid.$TEfamily.meth &
done < 1.prepare.bed/$genome.TE.list
wait

echo -e "TEfamily\tDNAm" > $sid.meth.per.TE.txt
for TEm in `ls $sid.*.meth`
do
    len=`awk '{sum+=$5}END{print sum}' $TEm`
    per=`awk '{sum+=$6}END{print sum}' $TEm`
    m=`echo "scale=8; $len / ($per + $len) " | bc`
    echo -e "$TEm\t$m" >> $sid.meth.per.TE.txt
done
```

## 3. Extract reads in TEs and calculate fragmentomics per sample
```
## split reads in each TE and calculate fragmentomics, use aligned BED file as input
## need to update the path to genome fasta file
genomefasta=/path/to/hg38.fa
calcFragPRG=2.fragmentomics/1.fragmentomics
genome=hg38
thread=16

while read TEclass TEfamily extra
do
    T=$TEclass.$TEfamily
    bedtools intersect -f 0.5 -a $sid.bed.gz -b $genome.TE/$T.bed.gz -wao -sorted | perl $calcFragPRG/deal.overlap.pl - >$sid.$T.bed &
    echo -e "$sid.$T\t$sid.$T.bed" >> $sid.bed.list
done < 1.prepare.bed/$genome.TE.list
wait
echo -e "$sid.Overall\t$sid.bed.gz" > $sid.bed.list
## analyze size and motif
perl $calcFragPRG/analyze.motif.size.multithread.pl $genomefasta $sid.bed.list $thread
perl $calcFragPRG/extract.size.pl      $sid 150  >$sid.size.pool
perl $calcFragPRG/extract.motif.pl     $sid CCCA >$sid.motif.pool
perl $calcFragPRG/extract.diversity.pl $sid 4    >$sid.diversity.pool
perl $calcFragPRG/extract.coverage.pl $genome.TE/TE.info $chrsize $sid >$sid.coverage.pool

## calculate E-index, for human data only
## Need to update the path to the end model, which is extremely large (~12GB).
## You can build it following the instructions in https://github.com/hellosunking/molecular-cfDNA-fragmentomics,
## or approach the corresponding author for the file.
EindexModel=/path/to/Eindex.model
hg38blacklist=$calcFragPRG/ENCODE.blacklist.hg38.bed
$calcFragPRG/calc.E-index.multi-thread $chrsize $EindexModel $hg38blacklist $sid.bed.list $thread y > $sid.E-index
```

## 4. Data visualizations and ROC analysis per cohort
```
## we had provided the compiled data for these plots
PRG=2.fragmentomics/3.visualization.ROC/
## boxplot for control samples (Fig. 1a)
## the fragmentomics results for 24 control samples are provided in "Processed.files" directory
cohort=humanNC24
Rscript $PRG/control.bar.plot.R $PRG/Processed.files/$cohort.size.pool $cohort.control.size
Rscript $PRG/control.bar.plot.R $PRG/Processed.files/$cohort.motif.pool $cohort.control.motif
Rscript $PRG/control.bar.plot.R $PRG/Processed.files/$cohort.diversity.pool $cohort.control.diversity
Rscript $PRG/control.coverage.bar.plot.R $PRG/Processed.files/$cohort.coverage.pool $cohort.control.coverage

## size distribution (Fig. 1b)
exampleSid=L01_2   ## this sample is choosen as an example
R --slave --args $cohort.size.pdf $PRG/Processed.files/$exampleSid.Overall.size $exampleSid.all \
$PRG/Processed.files/$exampleSid.LINE.L1.size $exampleSid.LINE.L1 \
$PRG/Processed.files/$exampleSid.SINE.Alu.size $exampleSid.SINE.Alu \
$PRG/Processed.files/$exampleSid.SINE.MIR.size $exampleSid.SINE.MIR \
$PRG/Processed.files/$exampleSid.LTR.ERVL.size $exampleSid.LTR.ERVL < $PRG/plot.size.log.R

## boxplot for cancer vs controls and ROC analyses (Fig. 4,5)
for cohortID in HCC Liang Cristiano Bie
do
    Rscript $PRG/box.plot.R $PRG/Processed.files/$cohortID.size.pool $cohortID.size
    Rscript $PRG/box.plot.R $PRG/Processed.files/$cohortID.motif.pool $cohortID.motif
	Rscript $PRG/box.plot.R $PRG/Processed.files/$cohortID.diversity.pool $cohortID.diversity
	Rscript $PRG/box.plot.R $PRG/Processed.files/$cohortID.coverage.pool $cohortID.coverage
	Rscript $PRG/box.plot.R $PRG/Processed.files/$cohortID.E-index.pool $cohortID.E-index

    Rscript $PRG/do.ROC.R $PRG/Processed.files/$cohortID.size.pool $cohortID.size
    Rscript $PRG/do.ROC.R $PRG/Processed.files/$cohortID.motif.pool $cohortID.motif
    Rscript $PRG/do.ROC.R $PRG/Processed.files/$cohortID.diversity.pool $cohortID.diversity
    Rscript $PRG/do.ROC.R $PRG/Processed.files/$cohortID.coverage.pool $cohortID.coverage
    Rscript $PRG/do.ROC.R $PRG/Processed.files/$cohortID.E-index.pool $cohortID.E-index
done
```

## 5. CfDNA fragmentomics associated with epigenetic features per sample
```
## DNA methylation
## we selected 20 control samples with highest sequencing depth from the Bie et al. dataset
## and merged the data to call cfDNA methylome.
## The list of these 20 samples are provided in "Bie.selected.20.controls.list"
cfDNAmethylome=/path/to/pooled.control.call
PRG=3.fragmentomics.epigenetics/1.DNAm/

## calculate DNA methylation in the given TE copies
while read TEclass TEfamily extra
do
    zcat $genome.TE/$TEclass.$TEfamily.bed.gz | $PRG/extract.meth.in.region $chrsize $cfDNAmethylome /dev/stdin > cfDNApool.$TEfamily.meth

    ## split TE copies into 3 categories based on the average methylation density of covered CpG sites
    perl -lane 'print if $F[3]>=1 && $F[4]/($F[5]+$F[4]) <= 0.2' cfDNApool.$TEfamily.meth > cfDNApool.$TEfamily.DNAm.low.bed
    perl -lane 'print if $F[3]>=1 && $F[4]/($F[5]+$F[4]) > 0.2 && $F[4]/($F[5]+$F[4]) <=0.9' cfDNApool.$TEfamily.meth > cfDNApool.$TEfamily.DNAm.medium.bed
    perl -lane 'print if $F[3]>=1 && $F[4]/($F[5]+$F[4]) > 0.9' cfDNApool.$TEfamily.meth > cfDNApool.$TEfamily.DNAm.high.bed
    python $PRG/split.TE.group.py $TEfamily.DNAm.group.size cfDNApool.$TEfamily.DNAm.*.bed 
done < 1.prepare.bed/$genome.TE.list

## for each TE family, split reads into the 3 groups defined by DNA methylation level
## here we use SINE/Alu as an example
TEclass=SINE
TEfamily=Alu

calcFragPRG=2.fragmentomics/1.fragmentomics
genomefasta=/path/to/hg38.fa
genome=hg38
thread=16

bedtools intersect -f 0.5 -a $sid.bed.gz -b cfDNApool.$TEfamily.DNAm.high.bed -wao -sorted | perl $calcFragPRG/deal.overlap.pl - >$sid.$TEclass.$TEfamily.DNAm.high.bed &
bedtools intersect -f 0.5 -a $sid.bed.gz -b cfDNApool.$TEfamily.DNAm.medium.bed -wao -sorted | perl $calcFragPRG/deal.overlap.pl - >$sid.$TEclass.$TEfamily.DNAm.medium.bed &
bedtools intersect -f 0.5 -a $sid.bed.gz -b cfDNApool.$TEfamily.DNAm.low.bed -wao -sorted | perl $calcFragPRG/deal.overlap.pl - >$sid.$TEclass.$TEfamily.DNAm.low.bed &
wait

echo -e "$sid.cfDNApool.$TEfamily.DNAm.high\t$sid.$TEclass.$TEfamily.DNAm.high.bed" > $sid.$TEfamily.DNAm.bed.list
echo -e "$sid.cfDNApool.$TEfamily.DNAm.medium\t$sid.$TEclass.$TEfamily.DNAm.medium.bed" >> $sid.$TEfamily.DNAm.bed.list
echo -e "$sid.cfDNApool.$TEfamily.DNAm.low\t$sid.$TEclass.$TEfamily.DNAm.low.bed" >> $sid.$TEfamily.DNAm.bed.list
echo -e "$sid.Overall\t$sid.bed.gz" >> $sid.$TEfamily.DNAm.bed.list
perl $calcFragPRG/analyze.motif.size.multithread.pl $genomefasta $sid.$TEfamily.DNAm.bed.list 3
perl $calcFragPRG/extract.size.pl      $sid.$TEfamily.DNAm 150  > $sid.$TEfamily.DNAm.size.pool
perl $calcFragPRG/extract.motif.pl     $sid.$TEfamily.DNAm CCCA > $sid.$TEfamily.DNAm.motif.pool
perl $calcFragPRG/extract.diversity.pl $sid.$TEfamily.DNAm 4    > $sid.$TEfamily.DNAm.diversity.pool
perl $calcFragPRG/extract.coverage.epigenetics.pl $TEfamily.DNAm.group.size $chrsize $sid.$TEfamily.DNAm $sid > $sid.$TEfamily.DNAm.coverage.pool

mkdir -p $sid.$TEclass.$TEfamily.DNAm
mv $sid.$TEfamily.DNAm.bed.list $sid.$TEclass.$TEfamily.*bed $sid.$TEclass.$TEfamily.*motif $sid.$TEclass.$TEfamily.*size $sid.$TEfamily.*pool $sid.$TEclass.$TEfamily.DNAm

## boxplots (Fig. 2a) after calculating the features for all control samples.
## we have provided the compiled matrices in "3.fragmentomics.epigenetics/1.DNAm" within this package
PRG=3.fragmentomics.epigenetics/1.DNAm
TEtoplot=SINE.Alu
Rscript $PRG/boxplot.R $PRG/Processed.files/NC24.Alu.DNAm.size.pool NC24.Alu.DNAm.size $TEtoplot
Rscript $PRG/boxplot.R $PRG/Processed.files/NC24.Alu.DNAm.motif.pool NC24.Alu.DNAm.motif $TEtoplot
Rscript $PRG/boxplot.R $PRG/Processed.files/NC24.Alu.DNAm.diversity.pool NC24.Alu.DNAm.diversity $TEtoplot
Rscript $PRG/boxplot.R $PRG/Processed.files/NC24.Alu.DNAm.coverage.pool NC24.Alu.DNAm.coverage $TEtoplot

## histone modifications (H3K27ac, H3K4me1, H3K4me3, H3K36me3, H3K27me3, and H3K9me3)
## and ATAC-seq (megakaryocytes, T-cells, and neutrophils)
## here we use SINE/Alu and H3K27ac (need to update path to cfChIP-seq data) as an example
targetEpi=H3K27ac
epigenomeBEDfile=/path/to/H3K27ac.bed
PRG=3.fragmentomics.epigenetics/2.histone.modification.ATAC/

## split TE into 2 groups based on the signal in epigenome
bedtools intersect -a $genome.TE/$TEclass.$TEfamily.bed.gz -F 0.5 -b $epigenomeBEDfile -wa -sorted -c > $TEfamily.$targetEpi.count
perl $PRG/split.TE.by.signal.pl $TEfamily.$targetEpi.count $TEfamily.$targetEpi
## the above program will write 3 files: "$TEfamily.$targetEpi.low.bed", "$TEfamily.$targetEpi.high.bed", "$TEfamily.$targetEpi.split.log"
python $PRG/split.TE.group.py $TEfamily.$targetEpi.group.size $TEfamily.$targetEpi.*.bed 

## for each TE family, split reads into the 2 groups defined by the epigenomic marker
genome=hg38
thread=16
calcFragPRG=2.fragmentomics/1.fragmentomics

bedtools intersect -f 0.5 -a $sid.bed.gz -b $TEfamily.$targetEpi.high.bed -wao -sorted | perl $calcFragPRG/deal.overlap.pl - >$sid.$TEfamily.$targetEpi.high.bed &
bedtools intersect -f 0.5 -a $sid.bed.gz -b $TEfamily.$targetEpi.low.bed -wao -sorted | perl $calcFragPRG/deal.overlap.pl - >$sid.$TEfamily.$targetEpi.low.bed &

echo -e "$sid.$TEfamily.$targetEpi.high\t$sid.$TEfamily.$targetEpi.high.bed" > $sid.$TEfamily.$targetEpi.bed.list
echo -e "$sid.$TEfamily.$targetEpi.low\t$sid.$TEfamily.$targetEpi.low.bed" >> $sid.$TEfamily.$targetEpi.bed.list
echo -e "$sid.Overall\t$sid.bed.gz" >> $sid.$TEfamily.$targetEpi.bed.list
perl $PRG/analyze.motif.size.multithread.pl $genomefasta $sid.$TEfamily.$targetEpi.bed.list 2
perl $PRG/extract.size.pl      $sid.$TEfamily.$targetEpi 150  >$sid.$TEfamily.$targetEpi.size.pool 
perl $PRG/extract.motif.pl     $sid.$TEfamily.$targetEpi CCCA >$sid.$TEfamily.$targetEpi.motif.pool 
perl $PRG/extract.diversity.pl $sid.$TEfamily.$targetEpi 4    >$sid.$TEfamily.$targetEpi.diversity.pool 
perl $PRG/extract.coverage.epigenetics.pl $TEfamily.$targetEpi.group.size $chrsize $sid.$TEfamily.$targetEpi $sid >$sid.$TEfamily.$targetEpi.coverage.pool 

mkdir -p $sid.$TEclass.$TEfamily.$targetEpi
mv $sid.$TEfamily.$targetEpi.bed.list  $sid.$TEfamily.*bed $sid.$TEfamily.*motif $sid.$TEfamily.*size $sid.$TEfamily.*pool $sid.$TEclass.$TEfamily.$targetEpi

## boxplots (Fig. 2a), we had provided the compiled matrix within this package
PRG=3.fragmentomics.epigenetics/2.histone.modification.ATAC/
Rscript $PRG/boxplot.R $PRG/Processed.files/NC24.H3K27ac.Alu.size.pool NC24.H3K27ac.Alu.size
Rscript $PRG/boxplot.R $PRG/Processed.files/NC24.H3K27ac.Alu.motif.pool NC24.H3K27ac.Alu.motif
Rscript $PRG/boxplot.R $PRG/Processed.files/NC24.H3K27ac.Alu.diversity.pool NC24.H3K27ac.Alu.diversity
Rscript $PRG/boxplot.R $PRG/Processed.files/NC24.H3K27ac.Alu.coverage.pool NC24.H3K27ac.Alu.coverage
```

## 6. Chromatin state analysis
Chromatin states were defined on cfChIP-seq data of 6 histone modifications (H3K27ac, H3K4me1, H3K4me3, H3K36me3, H3K27me3, and H3K9me3) segmented using chromHMM software.
```
## define chromatin states
## put (or link) the cfChIP-seq and cfDNA WGS data (as input) to "3.fragmentomics.epigenetics/3.chromatin.state/bed.file",
## and make sure the file names are consistent with that in "markers.table" file
PRG=3.fragmentomics.epigenetics/3.chromatin.state
java -mx32G -jar $PRG/ChromHMM.jar BinarizeBed -b 200 -f 2 -t out.signal $PRG/hg38.txt $PRG/bed.file/ $PRG/bed.file/markers.table $PRG/BinarizedTable
java -mx32G -jar $PRG/ChromHMM.jar LearnModel -r 1000 -b 200 -noautoopen -p 16 $PRG/BinarizedTable $PRG/State15 15 hg38
## Note that the chromatin states were manually annotated based on their histone profiles
## Our segmentation results and annotation are both provided within this package
ChromatinState15=3.fragmentomics.epigenetics/3.chromatin.state/cfDNA_15_segments.bed.gz
ChromatinStateAnno=3.fragmentomics.epigenetics/3.chromatin.state/cfDNA_15_segments.anno

## split each chromatin state
for CS in E{1..15}
do
    zgrep -w $CS $ChromatinState15 | sort -k1,1 -k2,2n >$PRG/State.$CS.bed
done
cd $PRG
python split.TE.group.py cfDNA_15_segments.size State.*.bed
cd ../../

## For each sample, split the reads into each chromatin state and calculate fragmentomics
PRG=3.fragmentomics.epigenetics/3.chromatin.state/
genome=hg38
genomefasta=/path/to/hg38.fa
thread=16

for CS in E{1..15}
do
    bedtools intersect -f 0.5 -a $sid.bed.gz -b $PRG/State.$CS.bed -wao -sorted | perl $calcFragPRG/deal.overlap.pl - > $sid.$CS.bed &
    echo -e "$sid.State.$CS\t$sid.$CS.bed" >> $sid.Chromatin.State.bed.list
done
wait
echo -e "$sid.Overall\t$sid.bed.gz" >> $sid.Chromatin.State.bed.list

## Extract fragmentomics per Chromatin State
calcFragPRG=2.fragmentomics/1.fragmentomics
perl $calcFragPRG/analyze.motif.size.multithread.pl $genomefasta $sid.Chromatin.State.bed.list $thread
perl $calcFragPRG/extract.size.pl      $sid.Chromatin.State 150  >$sid.Chromatin.State.size.pool 
perl $calcFragPRG/extract.motif.pl     $sid.Chromatin.State CCCA >$sid.Chromatin.State.motif.pool 
perl $calcFragPRG/extract.diversity.pl $sid.Chromatin.State 4    >$sid.Chromatin.State.diversity.pool 
perl $calcFragPRG/extract.coverage.epigenetics.pl $PRG/cfDNA_15_segments.size $chrsize $sid.Chromatin.State $sid >$sid.Chromatin.State.coverage.pool

## boxplots (Fig. 3b) on 24 controls, the compiled matrices are provided
Rscript $PRG/boxplot.R $PRG/Processed.files/NC24.chromatin.state.size.pool NC24.chromatin.state.size
Rscript $PRG/boxplot.R $PRG/Processed.files/NC24.chromatin.state.motif.pool NC24.chromatin.state.motif
Rscript $PRG/boxplot.R $PRG/Processed.files/NC24.chromatin.state.diversity.pool NC24.chromatin.state.diversity
Rscript $PRG/boxplot.R $PRG/Processed.files/NC24.chromatin.state.coverage.pool NC24.chromatin.state.coverage

## size distribution for different chromatin state (Fig. 3c)
exampleSid=L01_2
R --slave --args $exampleSid.size.pdf $PRG/Processed.files/$exampleSid.E5.size TssA \
$PRG/$exampleSid/$exampleSid.E2.size Enh1 \
$PRG/exampleSid/$exampleSid.E9.size Tx \
$PRG/exampleSid/$exampleSid.E13.size ReprPC \
$PRG/exampleSid/$exampleSid.E15.size Het \
$PRG/exampleSid/$exampleSid.E11.size Quies1 < $PRG/plot.mix.size.R
```

## 7. TEANA-Dx and TEANA-Top models
For whole genome sequencing data, we used the following cfDNA fragmentomic features in each TE family to build TEANA models: fraction of short fragments, CCCA end motif usage, end motif diversity, RSD, and E-index values. For Bie et al. dataset, DNA methylation densities were also included. The compiled feature matrices were provided within this package under "4.TEANA.models/" directory.
```
## TEANA-Dx diagnostic models
PRG=4.TEANA.models/

## Cristiano et al. dataset (Fig. 6a,b,c)
mkdir -p $PRG/Cristinao.Dx
Rscript $PRG/TEANA-Dx.R $PRG/Cristinao.matrix $PRG/Cristinao.Dx > Cristiano.Dx.log
python $PRG/average.py  $PRG/Cristinao.Dx/test_pred[0-9]* Cristiano.predict
Rscript $PRG/roc.R Cristiano.predict $PRG/Cristinao.Dx Cristinao.final > Cristiano.mean.roc.log
## statistics and ROCs
Rscript $PRG/TEANA-Dx.stage.R  Cristiano.predict $PRG/Cristinao.metadata
Rscript $PRG/TEANA-Dx.tissue.R Cristiano.predict $PRG/Cristinao.metadata

## Bie et al. dataset (Fig. 6d,e,f)
mkdir -p $PRG/Bie.Dx
Rscript $PRG/TEANA-Dx.R $PRG/Bie.matrix $PRG/Bie.Dx > Bie.Dx.log
python $PRG/average.py  $PRG/Bie.Dx/test_pred[0-9]* Bie.predict
Rscript $PRG/roc.R Bie.predict $PRG/Bie.Dx Bie.final > Bie.mean.roc.log
## statistics and ROCs
Rscript $PRG/TEANA-Dx.stage.R  Bie.predict $PRG/Bie.metadata
Rscript $PRG/TEANA-Dx.tissue.R Bie.predict $PRG/Bie.metadata

## TEANA-Top model
mkdir -p $PRG/Cristinao.Top
Rscript $PRG/TEANA-Top.R $PRG/Cristinao.matrix $PRG/Cristinao.Top > Cristiano.Top.log
python $PRG/get.multi.matrix.py Cristinao.multi.matrix Cristinao.cancer.matrix Bile Breast Colorectal Gastric Lung Ovarian Pancreatic
Rscript $PRG/stat.TEANA.top.R Cristinao.multi.matrix  > TEANA-Top.predict.stat
```
