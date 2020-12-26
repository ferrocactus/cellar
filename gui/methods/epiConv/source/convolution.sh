#!/bin/bash

##$1: head
##$2: tail
##$3: standard deviation

DIR=$(dirname $0)
HEAD=$1
TAIL=$2
SD=$3

zcat ${HEAD}_frag.bed.gz |\
	perl $DIR/bed_filter.pl ${HEAD}_ident.tsv - |\
	gawk 'BEGIN{ORS="";OFS="\t"}{print $1,$2,$2+1,$4"\n"$1,$3-1,$3,$4"\n"}' |\
	bedtools intersect -wo -a stdin -b ${HEAD}_peak.${TAIL} |\
	gawk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$8}' |\
	sort -k1,1 -k2n,2 \
	> ${HEAD}_inst.${TAIL}

NCELL=$(wc -l ${HEAD}_ident.tsv|gawk '{print $1}')
seq 1 $NCELL |\
	perl $DIR/lib_count.pl ${HEAD}_inst.${TAIL} - |\
	gawk '{print $2}' > ${HEAD}_lib.${TAIL}

bedtools window -l 0 -r $(echo "3*$SD"|bc) \
	-a ${HEAD}_inst.${TAIL} \
	-b ${HEAD}_inst.${TAIL} |\
	gawk '{OFS="\t"}{if($4!=$9&&$5==$10){
			if($4<$9) print $4,$9,$10,$7-$2;
			else print $9,$4,$10,$7-$2}}' |\
	$DIR/matrix_sum $NCELL $SD ${HEAD}_cmat.${TAIL}

$DIR/matrix_sampl ${HEAD}_cmat.${TAIL} ${HEAD}_ident.tsv $NCELL ${HEAD}_sampl.mtx \
	>${HEAD}_sampled.${TAIL}

rm ${HEAD}_inst.${TAIL}
