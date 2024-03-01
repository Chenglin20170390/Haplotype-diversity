## Variation calling (minimap2,samtools and syri)
<img width="300" alt="image" src="https://github.com/Chenglin20170390/Haplotype-diversity/assets/33062118/8000eb0b-bcec-4495-961d-f0c8a6517d04">

- Alignment（Genomes and Reads） for BAM file
```
for sample in $(cat list);do
    ref=DM.reference.fa
    for hap in H1 H2;do
    asm=$sample.$hap.assembly.fa  ##haplotypes
    ccs=$sample.$hap.hifi.fa   ##hifi reads
    th=10 ##threads for running
    minimap2 -ax asm5 -t $th --eqx $ref $asm | samtools sort -O BAM - > DM.$sample.$hap.bam
    samtools index DM.$sample.$hap.bam
    minimap2 -ax asm5 -t $th --eqx $ref $ccs | samtools sort -O BAM - > DM.$sample.$hap.ccs.bam
    samtools index DM.$sample.$hap.ccs.bam
    done
done
```
- Variation calling from syri based on BAM from assembly
```
for sample in $(cat list);do
    for hap in H1 H2;do
    ref=DM.reference.fa
    py38=/home/softwares/miniconda3/envs/py38/bin/python
    syri=/home/softwares/syri1.5/bin/syri
    plotsr=/home/softwares/syri/bin/plotsr
    $py38 $syri -c DM.$sample.$hap.bam -r $ref -q $asm -F B --prefix DM.$sample.$hap.
    $py38 $plotsr DM.$sample.$hap.syri.out $ref $asm -H 8 -W 5
    done
done    
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

