# dstats

An implementation of the D-statistics test for gene flow signal that works on consensus sequences. 

In addition to a file containing the number of ABBA and BABA sites, the script returns a second file that includes the number of BBAA sites (`-e`). To calculate the summary statistics from this file use the `jackknife_ext.R`.

The jackknife scripts are forked from the ANGSD repository, where also the original implementation of the D-statistics is implemented (https://github.com/ANGSD/angsd/).

Finally, the script returns a list of variable sites in the genomes as well as the coordinates of the ABBA, BABA and BBAA sites in BED format (`-b`).

**Note** Fasta identifiers in consensus sequences need to be identicical across samples/taxa.


Simple syntax for few sequences. 

```

python dstat.py -i A.fasta \                                                                                                                                        │
                -i B.fasta \                                                                                                                                        │
                -i C.fasta\                                                                                                                                        │Simple syntax for few sequences.                                                                                                                              
                -a OUT.fasta \                                                                                                                                 │
                -o out.abbababa                                                                                                                  │```
                -e out.abbabababbaa                                                                                                              │                                                                                                                                                              
                -b out.bed; 
```

## Parallization

In order to parallize the analyses for >4 taxa, I provide the `do_abbababa_parallel.sh`.

For each outgroup taxon, prepare a text file containing lines that contain the testname, samples and outgroup. 

E.g.

testfile
```
test1 taxonA.fasta taxonB.fasta taxonC.fasta outgroup.fasta
test2 taxonB.fasta taxonC.fasta taxonD.fasta outgroup.fasta
test3 taxonA.fasta taxonC.fasta taxonD.fasta outgroup.fasta
test4 taxonA.fasta taxonB.fasta taxonD.fasta outgroup.fasta
```

Then start a loop:

```

cat testfile  \
| while read line;
do
	nohup bash do_abbababa.parallel.sh $line&; 
done 
```

