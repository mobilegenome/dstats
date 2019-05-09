#!/bin/bash

F=$1
GENIC="/home/flammers/teddy/data/2015-02-05/umaritimus/polar_bear.filtered_1Mb.gene.sorted.fixed.genic_regions.sorted.gff"
INTERGENIC="/home/flammers/teddy/data/2015-02-05/umaritimus/polar_bear.filtered_1Mb.gene.sorted.fixed.intergenic.sorted.gff"
GENOME="/home/flammers/teddy/data/2015-02-05/umaritimus/BGI.scaf.filtered_1Mb.genome"

grep "h1-h2-h3" $F | grep -e "ABBA" -e "BABA"| bedtools sample -n 10000 -i ->  ${F/.bed/.h123.bed} &&
grep "h1-h2-h3" $F | grep  -e "BBAA" -e "SNV" | bedtools sample -n 10000 -i - > ${F/.bed/.h123.snv.bed} &&

#bedtools complement  -i ${F/.bed/.h123.bed} -g $GENOME | bedtools intersect -v -a - -b /home/flammers/teddy/data/2015-02-05/umaritimus/BGI.scaf.gaps.sorted.bed >  ${F/.bed/.other.bed} &&
#bedtools random -l 1 -g $GENOME |  bedtools intersect -v -a - -b /home/flammers/teddy/data/2015-02-05/umaritimus/BGI.scaf.gaps.sorted.bed  |  bedtools intersect -v -a - -b  ${F/.bed/.h123.bed}  >  ${F/.bed/.other.bed} &&
	
ab_genic=$(bedtools intersect -a ${F/.bed/.h123.bed} -b $GENIC  | wc -l) ;
ab_intergenic=$(bedtools intersect -a ${F/.bed/.h123.bed} -b $INTERGENIC  | wc -l) ; 
other_genic=$(bedtools intersect -a ${F/.bed/.h123.snv.bed} -b $GENIC  | wc -l) ; 
other_intergenic=$(bedtools intersect -a ${F/.bed/.h123.snv.bed} -b $INTERGENIC  | wc -l) ; 

echo "ABBABABA	genic	$ab_genic"
echo "ABBABABA	intergenic	$ab_intergenic"
echo "BBAA	genic	$other_genic"
echo "BBAA	intergenic	$other_intergenic"
echo "x <- matrix(c($ab_genic, $ab_intergenic, $other_genic, $other_intergenic),2,2)"

