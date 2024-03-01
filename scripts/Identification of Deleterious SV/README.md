# Identification of Deleterious SVs (dSVs)

<img width="300" alt="image" src="https://github.com/Chenglin20170390/Haplotype-diversity/assets/33062118/c9798f27-336d-4bb6-bc0d-ad357b2dca8e">

- Frequency of SV no more than 5%(three)
```
cat 60.hap.sv.fre.bed |awk '$6<=3' > 01_het_fre/sv.all.3.bed
```

- We need two files for identify dSV 1.[The CDS of Reference genome file](https://academic.oup.com/gigascience/article/9/9/giaa100/5910251) (DMv6.gff file) 2.[Constraint regions](https://doi.org/10.1016/j.cell.2023.04.008)（GERP>=2）
```
Here refrence file named  DM.gff3 ; constraint file named Sol_msa_AllChrs_GERP_withDepth.bed_Conserved2.info

##CDS文件
grep -w 'CDS' DM.gff3 |awk 'OFS="\t"{print $1,$4,$5,$3,$9}'| sort -k1,1 -k2,2n > DM.CDS.bed 

##Intron
grep -w mRNA DM.gff3 | awk 'OFS="\t"{print $1,$4,$5,$3,$9}'| sort -k1,1 -k2,2n > DM.mRNA.bed 
bedtools subtract -a DM.mRNA.bed -b DM.CDS.bed > DM.intron.bed  

##UTR
grep 'UTR' DM.gff3 |awk 'OFS="\t"{print $1,$4,$5,$3,$9}' > DM.UTR.bed 

##DM.UpDown5K.bed 
sed 's/"/\t/g'  DM.gff3 | awk 'BEGIN{OFS=FS="\t"}{if($3=="mRNA") {if($7=="+") {start=$4-5000; end=$4;} if(start<0) start=0; print $1,start,end,$7,$9;}}' |sort -k1,1 -k2,2n >DM.UpDown5K.bed
sed 's/"/\t/g'  DM.gff3 | awk 'BEGIN{OFS=FS="\t"}{if($3=="mRNA") {if($7=="-") {start=$5; end=$5+5000;} if(start<0) start=0; print $1,start,end,$7,$9;}}' |sort -k1,1 -k2,2n >>DM.UpDown5K.bed
sort -k1,1 -k2,2n DM.UpDown5K.bed  |uniq > DM.UpDown5K.sort.bed 

##Intergenic
genome=DM_v6.1_all_chr.fa
faidx ${genome} -i chromsizes |grep 'chr'> DM.genome.size
cat DM.mRNA.bed  DM.UpDown5K.sort.bed | sort -k1,1 -k2,2n |grep "chr"> DM.mRNA.UpDown5K.bed 
bedtools complement -i DM.mRNA.UpDown5K.bed -g DM.genome.size > DM.Intergenic.bed 
```

- Identification of dSV
```
cat DM.CDS.bed  Sol_msa_AllChrs_GERP_withDepth.bed_Conserved2.info|grep 'chr'|sed 's/\s\+/\t/g'|sort -k1,1 -k2,2n |awk 'OFS="\t"{print $1,$2,$3}' > CDS.constraint.bed
done

bedtools coverage -a 01_het_fre/sv.all.3.bed  -b CDS.constraint.bed  |sort |uniq > 02_bed/sv.all.3.CDS.constraint.bed 
```

- Identification of dSNP
```
for sample in $(cat list30);do
    for hap in H1 H2;do
    nohup comm -12 $sample.$hap.SNP.bed  CDS.constraint.bed   > $sample.$hap.dSNP.bed  &
    done
done
```
