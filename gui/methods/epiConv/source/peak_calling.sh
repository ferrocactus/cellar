#!/bin/bash

##$1: head
##$2: extsize
##$3: fraction of reads retained
##$4: whether clean temporary file

HEAD=$1
EXTSIZE=$2
FRAC=$3
CLEAN=$4

zcat ${HEAD}_frag.bed.gz |\
	gawk 'BEGIN{ORS="";OFS="\t"}{print $1,$2,$2+1,$4,".","+\n"$1,$3-1,$3,$4,".","-\n"}' >${HEAD}_inst.bed
macs2 pileup -i ${HEAD}_inst.bed -o ${HEAD}_pileup.txt -f BED -B --extsize $EXTSIZE

LIB=$(zcat ${HEAD}_frag.bed.gz|wc -l|gawk '{print $1}')
sort -k4nr,4 ${HEAD}_pileup.txt |\
	gawk 'BEGIN{sum=0}{sum+=($3-$2)*$4/'"$[$EXTSIZE*2]"';print $0"\t"sum/'"$[$LIB*2]"'}' >${HEAD}_density.bedgraph

gawk 'BEGIN{OFS="\t"}{if($5<'"$FRAC"') {$2=$2-'"$EXTSIZE"';$3=$3+'"$EXTSIZE"';if($2<0) $2=0;print $1,$2,$3}}' ${HEAD}_density.bedgraph |\
	sort -k1,1 -k2n,2 |\
	bedtools merge -d 0 -i stdin |\
	gawk '{print $0"\t"NR}' >${HEAD}_peak.bed

rm ${HEAD}_inst.bed \
   ${HEAD}_pileup.txt \
   ${HEAD}_density.bedgraph

