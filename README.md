# Haplotype-diversity
This repository contains code for the great haplotype diveristy in potato.

1. Genome survey for haplotype-resolved potato genome.

##for summary HiFi and Hic reads
seqkit stats $hifi.fastq.fa
seqkit stats $hic.fastq.R1.fa  $hic.fastq.R2.fa 
##genome size estimation (you need software of jellyfish and genomescope)
jellyfish count -C -m 21 -s 1000000000 -o $sample.reads.jf -t 10 <(zcat $sample.ccs.fastq.gz) 
jellyfish histo -t 10 $sample/$sample.reads.jf > $sample.reads.histo
genomescope.R  ../$sample/$sample.reads.histo 21 15000 ./

