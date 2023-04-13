# Haplotype-diversity
This repository contains code for the great haplotype diveristy in potato.

1. Genome survey .

##for summary HiFi and Hic reads (you need software ofseqkit )
```
seqkit stats $hifi.fastq.fa
seqkit stats $hic.fastq.R1.fa  $hic.fastq.R2.fa 
```
##genome size estimation (you need software of jellyfish and genomescope)

```jellyfish count -C -m 21 -s 1000000000 -o $sample.reads.jf -t 10 <(zcat $sample.ccs.fastq.gz) 
jellyfish histo -t 10 $sample/$sample.reads.jf > $sample.reads.histo
genomescope.R  ../$sample/$sample.reads.histo 21 15000 ./
```

2. Genome assembly.
```
##for initial assembly
hifiasm -o $assembly -t 20 $input.hifi.fastq.fa --h1 $hic_1 --h2 $hic_2
cat assembly.hifiasm.H1.fa assembly.hifiasm.H2.fa > assembly.Hapall.fa

##for merged_nodups.txt from 3ddna
    bwa index $assembly.fa
    python2 juicer-1.6/misc/generate_site_positions.py MboI assembly.Hapall  assembly.Hapall.fa
    awk 'BEGIN{OFS=\"\t\"}{print \$1,\$NF}' assembly.Hapall_MboI.txt > assembly.Hapall.chrom.size
    mkdir reference
    mkdir fastq 
    cd fastq
    ln -s $hic_dir/$name/${name}_R1.fq.gz ${name}_R1.fastq.gz
    ln -s $hic_dir/$name/${name}_R2.fq.gz ${name}_R2.fastq.gz
    cd ..
    scripts=/home/chenglin/softwares/juicer-1.6/CPU/common
    sh juicer-1.6/CPU/juicer.sh -s MboI -g assembly.Hapall -D $scripts -z assembly.Hapall.fa -p assembly.Hapall.chrom.sizes -y assembly.Hapall_MboI.txt

##Anchoring of initial contigs to the DM reference genome for easy manual correction of contigs.
ref=DM.reference .fa
ragtag.py scaffold $ref $assembly.hifiasm.H1.fa -o ragtag_H1
ragtag.py scaffold $ref $assembly.hifiasm.H2.fa -o ragtag_H2
##converting ragtag apg file to assembly format (agp2assembly.py from 3d-dna)
agp2assembly.py ragtag_H1/ragtag.scaffold.agp  ragtag_H1/$sample.H1.assembly
##Merging twao assembly into one for hic contact map visulization
python 01_2assemblyto1.py assembly.H1.assembly assembly.H2.assembly assembly.Hapall.assembly
/Users/cl/Desktop/genomics/code/01_sv/python/01_2assemblyto1.sort.py 
```



