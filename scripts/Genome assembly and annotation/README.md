# The general workflow for Genome assembly and annotation
<img width="600" alt="image" src="https://github.com/Chenglin20170390/Haplotype-diversity/assets/33062118/56ce3e4c-4ddc-4e00-894c-5a354d7489dd">

## 1. Genome survey.

- For summary HiFi and Hic reads (https://github.com/shenwei356/seqkit)
```
for sample in $(cat list);do
seqkit stats $sample.hifi.fastq.fa
seqkit stats $sample.hic.fastq.R1.fa  $sample.hic.fastq.R2.fa
done
```
- Genome size,heterozygousity,repeat elements estimation (https://github.com/gmarcais/Jellyfish and https://github.com/schatzlab/genomescope)

```
for sample in $(cat list);do
jellyfish count -C -m 21 -s 1000000000 -o $sample.reads.jf -t 10 <(zcat $sample.hifi.fastq.gz) 
jellyfish histo -t 10 $sample.reads.jf > $sample.reads.histo
genomescope.R  ../$sample.reads.histo 21 15000 ./
done
```

## 2. Genome assembly.

- For initial assembly (https://github.com/chhylp123/hifiasm)
```
for sample in $(cat list);do
hifiasm -o $sample.assembly -t 20 $sample.hifi.fastq.fa --h1 $sample.hic.fastq.R1.gz --h2 $sample.hic.fastq.R1.gz
cat $sample.assembly.hifiasm.H1.fa $sample.assembly.hifiasm.H2.fa > $sample.assembly.Hapall.fa
done
```
- For Hic contact information by alignment hic to phased assemblies before scaffording contig into phased pseudo-chromosome (https://github.com/aidenlab/3d-dna)
```
    for sample in $(cat list);do
    bwa index $sample.assembly.Hapall.fa
    python2 juicer-1.6/misc/generate_site_positions.py MboI $sample.assembly.Hapall  $sample.assembly.Hapall.fa
    awk 'BEGIN{OFS="\t"}{print $1,$NF}' $sample.assembly.Hapall_MboI.txt > $sample.assembly.Hapall.chrom.size
    mkdir reference
    mkdir fastq 
    cd fastq   ##put hic files in this directory
    ln -s $sample.hic.fastq.R1.gz ${sample}_R1.fastq.gz
    ln -s $sample.hic.fastq.R2.gz ${sample}_R2.fastq.gz
    cd ..
    scripts=/path to /juicer-1.6/CPU/common
    sh juicer-1.6/CPU/juicer.sh -s MboI -g $sample.assembly.Hapall -D $scripts -z $sample.assembly.Hapall.fa -p $sample.assembly.Hapall.chrom.sizes -y $sample.assembly.Hapall_MboI.txt
    done
```
- Anchoring of initial contigs to the DM reference genome for easy manual correction of contigs.(https://github.com/malonge/RagTag)
```
for sample in $(cat list);do
ref=DM.reference .fa  ##reference genome or relative species genome
ragtag.py scaffold $ref $sample.assembly.hifiasm.H1.fa -o ragtag_H1
ragtag.py scaffold $ref $asample.ssembly.hifiasm.H2.fa -o ragtag_H2
done
```
- Converting ragtag apg file to assembly format (agp2assembly.py from 3d-dna)
```
agp2assembly.py ragtag_H1/ragtag.scaffold.agp  ragtag_H1/$sample.H1.assembly
agp2assembly.py ragtag_H2/ragtag.scaffold.agp  ragtag_H2/$sample.H2.assembly
```
- Merging two assembly into one for hic contact map visulization(the script of 01_2assemblyto1.sort.py was writen based on 12 chromosme ,you may 
need to change based on your species)
```
python 01_2assemblyto1.sort.py assembly.H1.assembly assembly.H2.assembly assembly.Hapall.assembly
```

- Hic contact map for haplotype resolved genome (you may need to change -q to 1 for high heterozygous genome assembly)
```
3d-dna/visualize/run-assembly-visualizer.sh -q 0 assembly.Hapall.assembly merged_nodups.txt
```

- Then you  need to correct your hic contact map mannually by using juicerbox.
- Here are some methods that you may interseting in the correction of your assembly.


https://github.com/baozg/phased-assembly-check

https://www.youtube.com/watch?v=Nj7RhQZHM18&t=378s   ## Or bilibili  BV1B4411r77A


- Finally, you will get a haplotype-resolved genome by below command.
```
3d-dna/run-asm-pipeline-post-review.sh -q 0 -r assembly.Hapall.review.assembly assembly.Hapall.fa assembly.Hapall.merged_nodups.txt
```


## 2.1 Assembly assessment.

- Summary of N50,Gap, genome size et al  ([software assembly-stats](https://github.com/sanger-pathogens/assembly-stats))
```
for sample in $(cat list);do
assembly-stats  $sample.Hapall.fa > $sample.stats
done
```
- BUSCO  (software busco)
```
for sample in $(cat list);do
busco -m genome -i $sample.Hapall,fa -o busco_${sample}  --offline -l embryophyta_odb10 -c 10
done
```
-QV(https://github.com/lh3/yak)
```
for sample in $(cat list);do
yak count -b37 -t32 -o $sample.yak <(zcat ${sample}.R*.fq.gz) <(zcat ${sample}.R*.fq.gz)
yak qv -t32 -p -K3.2g -l100k $sample.yak $asm > $sample.qv.txt
done

```
- Switch and Hamming errors (https://github.com/tangerzhang/calc_switchErr)
- Generally, switch and Hamming errors were caculated by two phased vcf from reads and phased assemblies, respectively. https://github.com/tangerzhang/calc_switchErr pipline.
Then, we compared two types phased errors from whatshap compare fuction
```
whatshap compare $sample.reads.vcf $sample.assembly.vcf > $sample.phased_error.txt
```
- Reliable analysis of Flagger (https://github.com/mobinasri/flagger)，Here we not do optional step in Flagger piplines.


## 3. Genome annotation.
- For genome annotation, we used [the pipline](https://github.com/baozg/assembly-annotation-pipeline/tree/main) with a minor modification from Zhigui Bao's article below:

``
Bao Z, Li C, Li G, Wang P, Peng Z, Cheng L, Li H, Zhang Z, Li Y, Huang W, Ye M, Dong D, Cheng Z, VanderZaag P, Jacobsen E, Bachem CWB, Dong S, Zhang C, Huang S, Zhou Q. Genome architecture and tetrasomic inheritance of autotetraploid potato. Mol Plant. 2022 Jul 4;15(7):1211-1226. doi: 10.1016/j.molp.2022.06.009. Epub 2022 Jun 22. PMID: 35733345.
``

  
