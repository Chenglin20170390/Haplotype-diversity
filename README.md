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


- Alignment（Genomes and Reads） for BAM file
```
ref=DM.reference.fa
asm=$sample.hap1.assembly.fa
ccs=$sample.hap1.ccs.fa
sample=C001
minimap2 -ax asm5 -t $th --eqx $ref $asm | samtools sort -O BAM - > DM.$sample.$hap.bam
samtools index DM.$sample.$hap.bam
minimap2 -ax asm5 -t $th --eqx $ref $ccs | samtools sort -O BAM - > DM.$sample.$hap.ccs.bam
samtools index DM.$sample.$hap.ccs.bam
```
- Variation calling from syri based on BAM from assembly
```
py38=/home/softwares/miniconda3/envs/py38/bin/python
syri=/home/softwares/syri1.5/bin/syri
plotsr=/home/softwares/syri/bin/plotsr
$py38 $syri -c DM.$sample.$hap.bam -r $ref -q $asm -F B --prefix DM.$sample.$hap.
$py38 $plotsr DM.$sample.$hap.syri.out $ref $asm -H 8 -W 5
```
- Variant  filtering based on INS DEL INV DUP at column 11.
```
mkdir 01_bed
for sample in $(cat  list30);do
    for hap in H1 H2;do
    awk '$11=="DEL" && $3-$2>49'  $sample.$hap.syri.out  > 01_bed/$sample.$hap.50.bed
    awk '$11=="INS" && $8-$7>49'  $sample.$hap.syri.out  >> 01_bed/$sample.$hap.50.bed
    awk '$11=="INV" && $3-$2>49'  $sample.$hap.syri.out   >> 01_bed/$sample.$hap.50.bed
    awk '$11=="DUP" && $3-$2>49'  $sample.$hap.syri.out   >> 01_bed/$sample.$hap.50.bed
    done
done

```
- Prepare SV position file for reads coverage detection from DM.$sample.$hap.ccs.bam
```
mkdir -p 05_ccs_cov/01_cov
for sample in $(cat list30);do
    for bin in 0 ;do
    size=$bin
        for hap in H1 H2;do
        ##for all sv
        awk -v bin=$size  'OFS="\t"{print $1,$2,$9}' 01_bed/$sample.$hap.50.bed > 05_ccs_cov/01_cov/$sample.$hap.sv.up.$bin.bed
        awk -v bin=$size  'OFS="\t"{print $1,$3,$9}' 01_bed/$sample.$hap.50.bed  > 05_ccs_cov/01_cov/$sample.$hap.sv.dn.$bin.bed
        done
    done
done
```
- read coverage detection from read.BAM (reads map quality >=50 are used for coverage caculation)
```
for sample in $(cat  /home/chenglin/work/01_sv/01_genome/02_hap/list30);do 
sample='$sample'
bin=0
export PATH=/public/software/env01/bin/:$PATH
for hap in H1 H2;do
samtools depth -b 01_cov/$sample.$hap.sv.up.$bin.bed -Q 50 -a  01_bam/$sample.dm.$sample.bam -@ 2 > 04_bk_cov/$sample.$hap.up.123.Q50.cov
samtools depth -b 01_cov/$sample.$hap.sv.dn.$bin.bed  -Q 50 -a 01_bam/$sample.dm.$sample.bam -@ 2 > 04_bk_cov/$sample.$hap.dn.123.Q50.cov
done
```

- merge up.bed and dn.bed into one file
```
for sample in $(cat  list30);do 
    for hap in H1 H2;do
    cat 04_bk_cov/$sample.$hap.up.123.Q50.cov 04_bk_cov/$sample.$hap.dn.123.Q50.cov  |sort |uniq > 04_bk_cov/$sample.$hap.all.123.Q50.cov
    done
done
```

- filter pass SV (reads coverage >=3)
```
mkdir 05_ftsv_bed
for sample in $(cat  list30);do 
    for hap in H1 H2;do
    python 40.sv.cov.py 04_bk_cov/$sample.$hap.all.123.Q50.cov  ../01_bed/$sample.$hap.bed 05_ftsv_bed/$sample.$hap.ft.bed 
    done
done
```
- Prepare merge_list for SV merge
```
for sample in $(cat list30);do 
export PATH=/home//softwares/SURVIVOR/Debug:$PATH
    for type in INS DEL INV DUP;do
        for hap in H1 H2;do
        ##bed to vcf
        SURVIVOR  bedtovcf 01_bed/$sample.$hap.50.$type.bed $type 02_bed2vcf/$sample.$hap.$type.survivor.vcf
        echo $sample $hap
        sed  -i 's/Sample/'${sample}_$hap'/g' 02_bed2vcf/$sample.$hap.$type.survivor.vcf
        done
    sed -i 's/.\/./1|0/g' 02_bed2vcf/$sample.H1.$type.survivor.vcf
    sed -i 's/.\/./0|1/g' 02_bed2vcf/$sample.H2.$type.survivor.vcf
    ##merge SV based on 80% overlap   
    echo  02_bed2vcf/$sample.H1.$type.survivor.vcf >> merge_list_$type
    echo  02_bed2vcf/$sample.H2.$type.survivor.vcf >>  merge_list_$type
    done
done
```
- merge different SV types
```
for type in INS DEL INV DUP;do
SURVIVOR merge merge_list_$type $key  1 1 1 -1 50 01_merge/61.Hapall.$key.$type.vcf
sed -i 's/.\/./0|0/g'  01_merge/61.Hapall.$key.$type.vcf
done
```

- merge ALL SV
```
for type in INS DEL INV DUP;do
vcf-sort  61.Hapall.$type.vcf > 61.Hapall.$type.sort.vcf
rm 61.Hapall.$type.sort.vcf.gz
bgzip  61.Hapall.$type.sort.vcf
tabix 61.Hapall.$type.sort.vcf.gz 
bcftools view -S samples.txt 61.Hapall.$type.sort.vcf.gz  > 61.Hapall.$type.sort.order.vcf
done
cat 61.Hapall.INS.sort.order.vcf > 61.Hapall.mergeall.vcf 
grep -v "#"   61.Hapall.DEL.sort.order.vcf   >>61.Hapall.mergeall.vcf 
grep -v "#"   61.Hapall.INV.sort.order.vcf   >>61.Hapall.mergeall.vcf 
grep -v "#"   61.Hapall.DUP.sort.order.vcf   >>61.Hapall.mergeall.vcf 
vcf-sort 61.Hapall.mergeall.vcf > 61.Hapall.mergeall.sort.vcf
```







