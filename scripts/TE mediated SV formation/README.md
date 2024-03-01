# TE mediated SV formation
<img width="600" alt="image" src="https://github.com/Chenglin20170390/Haplotype-diversity/assets/33062118/928373ad-bc5d-414f-8595-0200cf8c84b2">

Reference:Balachandran, P., Walawalkar, I.A., Flores, J.I. et al. Transposable element-mediated rearrangements are prevalent in human genomes. Nat Commun 13, 7115 (2022).

## Identification of repeats elements
- Transposable elements (TEs)
```
for sample in $(cat list30);do
    for hap in H1 H2;do
    th=36
    asm=$sample.$hap.fa  #assembly genome file
    species=others #others if your species neither rice or maize
    sensitive=0 #1 for running repeatmodeler/ 0 for not run
    EDTA.pl  --genome $asm --species others  --sensitive $sensitive  --anno 1  --threads $th  > edta.log 2>&1
    done
done

##Then convert gff to bed file
grep -v "Parent" $sample.$hap.chr.fa.mod.EDTA.TEanno.gff3 |awk '{OFS="\t"}{print $1,$4,$5,$3}' > $sample.$hap.te.bed

```
- Tandem Repeats(TRs)
```
for sample in $(cat ist30);do
    for hap in H1 H2;do
    nohup trf $sample.$hap.chr.fa 2 6 6 80 10 50 2000 -h  &
    sleep 1s
    done
done
```
  Then you can use [TRF2GFF](https://github.com/Adamtaranto/TRF2GFF) to convert TRF to GFF and bed format. 
```
for sample in $(cat list30);do
    for hap in H1 H2;do
    python  TRF2GFF.py -d $sample.$hap.chr.fa.2.6.6.80.10.50.2000.dat \
    -o $sample.$hap.tr.gff
    grep -v '^#' $sample.$hap.tr.gff|sed '/^$/d' | awk '{OFS="\t"}{print $1,$4,$5,$9}' > $sample.$hap.tr.bed

    done
done
```

- Segmental duplications (SDs)
```
for sample in $(cat list30);do
    for hap in H1 H2;do
    th=20
    export PATH=/home/chenglin/softwares/asgart_dir:$PATH
    asgart 01_soft_mask/$sample.$hap.chr.softmask.fa --threads $th -CRSv  --out $sample.$hap
    asgart-plot 02_sd_res/$sample.$hap.json  chord  --min-length 5000 --out=$sample.$hap.svg
    awk '{OFS="\t"}{print $1,$4,$5,$9}'  03_sd_gff/$sample.$hap.gff3 > $sample.$hap.sd.bed
    done
done
```
- SV bed file (from syri) for overlap with Repeat elements(TEs,SDs and TRs)
```
##Only  INS DEL INV DUP
for svtype in INS DEL INV DUP;do
mkdir -p 03_coverage/$svtype
    for sample in $(cat list30);do
        for hap in H1 H2;do
            for bin in  100 ;do
                for hap in H1 H2;do
                    for type in TE LTR TIR HEL SD TR ALL;do
                    bedtools coverage -a 02_svind_bed/$sample.$hap.$svtype.up.$bin.bed -b 01_bk_feature/$sample.$hap.$type.bed -wo > 03_coverage/$svtype/$sample.$hap.up.$bin.$type.bed
                    bedtools coverage -a 02_svind_bed/$sample.$hap.$svtype.dn.$bin.bed -b 01_bk_feature/$sample.$hap.$type.bed -wo > 03_coverage/$svtype/$sample.$hap.dn.$bin.$type.bed
                    paste  03_coverage/$svtype/$sample.$hap.up.$bin.$type.bed  03_coverage/$svtype/$sample.$hap.dn.$bin.$type.bed  > 03_coverage/$svtype/$sample.$hap.all.$bin.$type.bed
                    done
                done
            done
        done
    done
    ' > $svtype.sh
    nohup  sh $svtype.sh &
done


```
- Simulation SVs and overlap with Repeat elements. 
```
##simualtion SV
bin=100
for svtype in INS DEL INV DUP;do
    for sample in $(cat list30);do
        for hap in H1 H2;do
        awk '{OFS="\t"}{print $1,$2}' $sample.$hap.chr.fa.fai >01_genome/$sample.$hap.genome
        ####make comlement for repeat region (to speed up)
        bedtools shuffle -i  ../02_svind_bed/$sample.$hap.$svtype.up.100.bed  -g 01_genome/$sample.$hap.genome -excl 02_svind_bed/$sample.$hap.svall.bed > 02_svind_bed/$sample.$hap.$svtype.up.100.bed
        bedtools shuffle -i  ../02_svind_bed/$sample.$hap.$svtype.dn.100.bed  -g 01_genome/$sample.$hap.genome -excl 02_svind_bed/$sample.$hap.svall.bed > 02_svind_bed/$sample.$hap.$svtype.dn.100.bed
        done
    done
done

##Overlap with Repeat elements
for svtype in INS DEL INV DUP;do
    for sample in $(cat list30);do
        for hap in H1 H2;do
            for bin in  100 ;do
                for type in TE LTR TIR HEL SD TR ALL;do
                bedtools coverage -a 02_svind_bed/$sample.$hap.$svtype.up.$bin.bed -b 01_bk_feature/$sample.$hap.$type.bed -wo > 03_coverage/$svtype/$sample.$hap.up.$bin.$type.bed
                bedtools coverage -a 02_svind_bed/$sample.$hap.$svtype.dn.$bin.bed -b 01_bk_feature/$sample.$hap.$type.bed -wo > 03_coverage/$svtype/$sample.$hap.dn.$bin.$type.bed
                paste  03_coverage/$svtype/$sample.$hap.up.$bin.$type.bed  03_coverage/$svtype/$sample.$hap.dn.$bin.$type.bed  > 03_coverage/$svtype/$sample.$hap.all.$bin.$type.bed
                done
            done
        done
    done
    nohup  sh $svtype.sh &
done
```


  
