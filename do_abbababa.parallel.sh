#!/usr/bin/python

outgroup=$5
echo $1 $2 $3 $4 $outgroup
outgroup_name=$(echo $outgroup | cut -f 1 -d ".")
python dstat.py -i $2 \
		-i $3 \
		-i $4 \
		-a $outgroup \
		-o $1.$outgroup_name.abbababa 
		-e $1.$outgroup_name.abbabababbaa 
		-b $1.$outgroup_name.bed; 

#Rscript jackKnife_fritjof.R file=$1.$outgroup_name.abbabababbaa indNames=$1.flist outfile=$1.$outgroup_name.bbaa.jackknife ;
#Rscript ~/programs/angsd/R/jackKnife.R file=$1.$outgroup_name.abbababa indNames=$1.flist outfile=$1.$outgroup_name.jackknife ;
