all: *.balMys.abbababa *.Egl00.abbababa *.Bac00.abbababa *.Bmu.abbababa  *.Bbo02.abbababa 
#	*.Egl00.jackknife.txt *.Bac00.jackknife.txt *.Bmu.jackknife.txt *.balMys.jackknife.txt *.Bbo02.jackknife.txt

.precious: *.jackknife.txt

*.balMys.abbababa *.balMys.jackknife.txt: outgroup_bowhead
	while read line; do nohup bash do_abbababa.parallel.sh $$line balMys.scaf.filtered_100kb.fa & done < $<

*.Egl00.abbababa *.Egl00.jackknife.txt: outgroup_natr
	while read line; do nohup bash do_abbababa.parallel.sh $$line Egl00.sorted.mkdup.consensus.fasta & done < $<

*.Bac00.abbababa *.Bac00.jackknife.txt: outgroup_minke
	while read line; do nohup bash do_abbababa.parallel.sh $$line Bac00.sorted.mkdup.consensus.fasta & done < $<

*.Bbo02.abbababa *.Bbo02.jackknife.txt: outgroup_sei
	while read line; do nohup bash do_abbababa.parallel.sh $$line Bbo02.sorted.mkdup.consensus.fasta & done < $<

*.Bmu.abbababa *.Bmu.jackknife.txt: outgroup_blue
	while read line; do nohup bash do_abbababa.parallel.sh $$line Bmu.merged.consensus.fasta & done < $<


