# Haplotype-diversity
This repository contains workflow and related code for the analysis of great haplotype diveristy in potato.

## 1. Genome survey.

- For summary HiFi and Hic reads (https://github.com/shenwei356/seqkit)
```
seqkit stats $hifi.fastq.fa
seqkit stats $hic.fastq.R1.fa  $hic.fastq.R2.fa 
```
- Genome size,heterozygousity,repeat elements estimation (https://github.com/gmarcais/Jellyfish and https://github.com/schatzlab/genomescope)

```
jellyfish count -C -m 21 -s 1000000000 -o $sample.reads.jf -t 10 <(zcat $sample.ccs.fastq.gz) 
jellyfish histo -t 10 $sample/$sample.reads.jf > $sample.reads.histo
genomescope.R  ../$sample/$sample.reads.histo 21 15000 ./
```

## 2. Genome assembly.

- For initial assembly (https://github.com/chhylp123/hifiasm)
```
hifiasm -o $assembly -t 20 $input.hifi.fastq.fa --h1 $hic_1 --h2 $hic_2
cat assembly.hifiasm.H1.fa assembly.hifiasm.H2.fa > assembly.Hapall.fa
```
- for merged_nodups.txt from 3ddna (https://github.com/aidenlab/3d-dna)
```
    bwa index $assembly.fa
    python2 juicer-1.6/misc/generate_site_positions.py MboI assembly.Hapall  assembly.Hapall.fa
    awk 'BEGIN{OFS="\t"}{print $1,$NF}' assembly.Hapall_MboI.txt > assembly.Hapall.chrom.size
    mkdir reference
    mkdir fastq 
    cd fastq
    ln -s $hic_dir/$name/${name}_R1.fq.gz ${name}_R1.fastq.gz
    ln -s $hic_dir/$name/${name}_R2.fq.gz ${name}_R2.fastq.gz
    cd ..
    scripts=/home/chenglin/softwares/juicer-1.6/CPU/common
    sh juicer-1.6/CPU/juicer.sh -s MboI -g assembly.Hapall -D $scripts -z assembly.Hapall.fa -p assembly.Hapall.chrom.sizes -y assembly.Hapall_MboI.txt
```
- Anchoring of initial contigs to the DM reference genome for easy manual correction of contigs.(https://github.com/malonge/RagTag)
```
ref=DM.reference .fa
ragtag.py scaffold $ref $assembly.hifiasm.H1.fa -o ragtag_H1
ragtag.py scaffold $ref $assembly.hifiasm.H2.fa -o ragtag_H2
```
- Converting ragtag apg file to assembly format (agp2assembly.py from 3d-dna)
```
agp2assembly.py ragtag_H1/ragtag.scaffold.agp  ragtag_H1/$sample.H1.assembly
agp2assembly.py ragtag_H2/ragtag.scaffold.agp  ragtag_H2/$sample.H2.assembly
```
- Merging twao assembly into one for hic contact map visulization(the script of 01_2assemblyto1.sort.py was writen based on 12 chromosme ,you may 
need to change based on your species)
```
python 01_2assemblyto1.sort.py assembly.H1.assembly assembly.H2.assembly assembly.Hapall.assembly
```

- Hic contact map for haplotype resolved genome (you may nned to change -q to 1 for high heterozygous genome assembly)
```
3d-dna/visualize/run-assembly-visualizer.sh -q 0 assembly.Hapall.assembly merged_nodups.txt
```

- Then you  need to correct your hic contact map mannually by using juicerbox.
- Here are some methods that you may interseting in the correction of your assembly.


https://github.com/baozg/phased-assembly-check

https://www.youtube.com/watch?v=Nj7RhQZHM18&t=378s   ## Or bilibili  BV1B4411r77A


- Finally, you will get a haplotype-resolved genome.
```
3d-dna/run-asm-pipeline-post-review.sh -q 0 -r assembly.Hapall.review.assembly assembly.Hapall.fa assembly.Hapall.merged_nodups.txt
```


## 2.1 Assembly assessment.

- N50,gap and size summary ([software assembly-stats](https://github.com/sanger-pathogens/assembly-stats))
```
assembly-stats  $dir/$sample.$hap.chr.fa >01_ind/$sample.$hap.N.stats
```
- BUSCO  (software busco)
```
busco -m genome -i $asm -o busco_${sample}  --offline -l embryophyta_odb10 -c $th
```

- Switch and hamming errors (https://github.com/tangerzhang/calc_switchErr)
- You could also used  WhatsHap compare function for this purpose.



## 3. Genome annotation.
- For genome annotation, we used pipline from Zhigui Bao's article below:
``
Bao Z, Li C, Li G, Wang P, Peng Z, Cheng L, Li H, Zhang Z, Li Y, Huang W, Ye M, Dong D, Cheng Z, VanderZaag P, Jacobsen E, Bachem CWB, Dong S, Zhang C, Huang S, Zhou Q. Genome architecture and tetrasomic inheritance of autotetraploid potato. Mol Plant. 2022 Jul 4;15(7):1211-1226. doi: 10.1016/j.molp.2022.06.009. Epub 2022 Jun 22. PMID: 35733345.
``

## 4. Variation calling (minimap2,samtools and syri)
- SV simulation based on diploid potato genome


- Alignment for BAM file
```
ref=DM.reference.fa
asm=assembly.hap1.fa
sample=C001
minimap2 -ax asm5 -t $th --eqx $ref $asm | samtools sort -O BAM - > DM.$sample.$hap.bam
samtools index 01_bam/$sample.dm.$sample.$hap1.bam
```
- Variation calling from syri
```
$py38 $syri -c DM.$sample.$hap.bam -r $ref -q $asm -F B --prefix DM.$sample.$hap.
$py38 $plotsr DM.$sample.$hap.syri.out $ref $asm -H 8 -W 5
```
- Variant  filtering based on INS DEL INV DUP at column 11.
```

```





