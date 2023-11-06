# Building a phased potato pangenome

## Rename by Pan-SN
Follow the PanSN nameing scheme, we rename all the sequence to `A157#1#chr01` such as this.

## Estimation of divergence

```bash
for i in `seq 1 12`;do mash triangle panPotato.chr${i}.fa.gz > chr${i}.mdist;done
for i in `seq 1 12`;do sed 1,1d chr${i}.mdist | tr '\t' '\n' | grep -v "chr"|sort -k 1gr | head -n 5 >> all.top;done

sort -k1gr all.top|head -n 5
0.0588483
0.0586096
0.0586096
0.0583726
0.0583726

```

## Chromosome community

```bash
singularity run -B ${PWD}:/data pggb-0.5.4.sif /bin/bash -c "wfmash -m -p 90 -s 20000 -n 61 -t 128 /data/panPotato.chunked.1Mb.fasta.gz > /data/panPotato.p90s20kn61.paf"

python ~/software/pggb/scripts/paf2net.py -p panPotato.p90s20kn61.paf

# color and clustering in Gephi

```

## Building pangenome graph

### PGGB
```bash
#f=panPotato.chr1.fa.gz
f=$1
p=90
s=10000
n=61
k=47
G=700,900,1100
ref="DM"
out=$(basename "$f" .fa.gz)_p$p.s$s.n$n.k$k.G$(echo $G | tr ',' '-').$ref
t=128
O=0.001
POA=asm20

# PGGB
singularity run -B ${PWD}/data:/data pggb-0.5.4.sif /bin/bash -c "pggb -i /data/$f -p $p -s $s -n $n -k $k -P $POA -O $O -G $G -t $t -v -o /data/$out"


# vis

ch=$1
fmt="%C\n%Us user %Ss system %P cpu %es total %MKb max memory"
timer=$(which time)

prefix_smoothed_output="panPotato.chr${ch}.fa.gz.47c7dfd.e34d4cd.98ebc75.smooth"
log_file="${prefifx_smmothed_output}.log"

#color by nucl pos in the path
singularity run -B ${PWD}:/data pggb-0.5.4.sif /bin/bash -c "odgi viz -i /data/${prefix_smoothed_output}.final.og -t 24 -o /data/${prefix_smoothed_output}.final.og.viz_pos_multiqc.png -x 1500 -y 500 -a 10 -u -d -I \"Consensus_\" 2> >(tee -a /data/\"$log_file\")"

# color by mean depth
singularity run -B ${PWD}:/data pggb-0.5.4.sif /bin/bash -c "odgi viz -i /data/${prefix_smoothed_output}.final.og -t 24 -o /data/${prefix_smoothed_output}.final.og.viz_depth_multiqc.png -x 1500 -y 500 -a 10 -m -I \"Consensus_\" 2> >(tee -a /data/\"$log_file\")"

# color by mean inversion rate
singularity run -B ${PWD}:/data pggb-0.5.4.sif /bin/bash -c "odgi viz -i /data/${prefix_smoothed_output}.final.og -t 24 -o /data/${prefix_smoothed_output}.final.og.viz_inv_multiqc.png -x 1500 -y 500 -a 10 -z -I \"Consensus_\" 2> >(tee -a /data/\"$log_file\")"

# uncalled
#singularity run -B ${PWD}:/data pggb-0.5.4.sif /bin/bash -c "odgi viz -i /data/${prefix_smoothed_output}.final.og -o /data/${prefix_smoothed_output}.final.og.viz_uncalled_multiqc.png -x 1500 -y 500 -a 10 -N -I \"Consensus_\" 2> >(tee -a /data/\"$log_file\")"

# pav
#singularity run -B ${PWD}:/data pggb-0.5.4.sif /bin/bash -c "odgi pav -t 48 -i /data/${prefix_smoothed_output}.final.og -b /data/DM.chr${ch}.100kb.bed > /data/${prefix_smoothed_output}.w100kbp.pavs.tsv":w

# stats
singularity run -B ${PWD}:/data pggb-0.5.4.sif /bin/bash -c "odgi stats -t 48 -S -i /data/${prefix_smoothed_output}.final.og > /data/${prefix_smoothed_output}.odgi.stats.tsv"G```

#### Variants


```bash
# deconstruct
vg deconstruct -P DM -H '#' -e -a -t 24 ${prefix_smoothed_output}.final.gfa|bgzip -@ 48 > ${prefix_smoothed_output}.vcf.gz

# split SNPs and SVs
parallel -j 12 'bcftools view --threads 4 -V snps panPotato.chr{}.fa.gz.47c7dfd.e34d4cd.98ebc75.smooth.vcf.gz -O z -o SVs/panPotato.chr{}.vcf.gz' ::: `seq 1 12`
parallel -j 12 'bcftools view --threads 4 -v snps panPotato.chr{}.fa.gz.47c7dfd.e34d4cd.98ebc75.smooth.vcf.gz -O z -o snps/panPotato.chr{}.vcf.gz' ::: `seq 1 12`

# SVs

## InDels (<50bp)
parallel -j 12 'vcfbub -i SVs/panPotato.chr{}.vcf.gz -A 1 -a 50 -l 0 |bgzip -@ 6 -c > panPotato.chr{}.indels.vcf.gz' ::: `seq 1 12`
parallel -j 12 'bcftools stats -s "-" panPotato.chr{}.indels.vcf.gz > panPotato.indels.chr{}.stats' ::: `seq 1 12`


## SVs (50bp < SVs < 1Mb)
parallel -j 12 'vcfbub -i SVs/panPotato.chr{}.vcf.gz -A 50 -a 1000000 -l 0 |bgzip -@ 6 -c > panPotato.chr{}.vcfbub_50bp_1Mb.vcf.gz' ::: `seq 1 12`
parallel -j 12 'bcftools stats -s "-" panPotato.chr{}.vcfbub_50bp_1Mb.vcf.gz > panPotato.vcfbub_50bp_1Mb.chr{}.stats' ::: `seq 1 12`

## 1Mb < SVs < 100Mb (large inversion)
parallel -j 12 'vcfbub -i SVs/panPotato.chr{}.vcf.gz -A 1000000 -a 100000000 -l 0 |bgzip -@ 6 -c > panPotato.chr{}.vcfbub_lt1Mb.vcf.gz' ::: `seq 1 12`
parallel -j 12 'bcftools stats -s "-" panPotato.chr{}.vcfbub_lt1Mb.vcf.gz > panPotato.lt1Mb.chr{}.stats' ::: `seq 1 12`

```

#### Non-ref paths

```bash

# non-ref paths
singularity exec -B ${PWD}:/data odgi-non-ref.sif /bin/bash -c "odgi paths -i /data/${prefix_smoothed_output}.final.og -t 64 --non-reference-ranges /data/DM.chr${ch}.list" > /data/panPotato.nonref.chr${ch}.bed

# split by sample and hap
perl split.pl panPotato.nonref.chr${ch}.bed

# >50bps nodes
for sam in `ls *_H*.non-ref.bed|sed "s/.non-ref.bed//g"`;do awk '$3-$2>50' ${sam}.non-ref.bed|sort -k1V -k2n > sort/${sam}.50bp.bed


# overlap

sample=$1
bedtools intersect -a sort/${sample}_H1.50bp.bed -b filter/${sample}_H1.TE.gff3 -wa -wb|cut -f4-20 > sort/${sample}_H1.TE.tmp.gff3
cat sort/${sample}_H1.TE.tmp.gff3|./gff2out.pl > sort/${sample}_H1.TE.tmp.out
perl ~/miniconda3/envs/EDTA/share/EDTA/util/buildSummary.pl sort/${sample}_H1.TE.tmp.out > sort/${sample}_H1.TE.tbl

rm sort/${sample}_H1.TE.tmp.gff3 sort/${sample}_H1.TE.tmp.out

bedtools intersect -a sort/${sample}_H2.50bp.bed -b filter/${sample}_H2.TE.gff3 -wa -wb|cut -f4-20 > sort/${sample}_H2.TE.tmp.gff3
cat sort/${sample}_H2.TE.tmp.gff3|./gff2out.pl > sort/${sample}_H2.TE.tmp.out
perl ~/miniconda3/envs/EDTA/share/EDTA/util/buildSummary.pl sort/${sample}_H2.TE.tmp.out > sort/${sample}_H2.TE.tbl

rm sort/${sample}_H2.TE.tmp.gff3 sort/${sample}_H2.TE.tmp.out

```



### Minigraph-Cactus

```bash
# Minigraph-Cactus
cactus-bin-v2.5.2/cactus_env/bin/cactus-pangenome js genome.pos.txt --outDir out --outName potato_pg --reference DM --vcf --giraffe --gfa --gbz --workDir $PWD --restart

# Stats


# Covered regions

zcat potato_pg.gfa.gz|grep "^W"|cut -f1-6 > potato_pg.MC.keep.bed

perl scripts/stat_MC.pl potato_pg.MC.keep.bed group.xls > potato_pg.MC.keep.perchr.stats.tsv

```

