# Correlation between dSV and dSNP

- Simulation non-deleterious SV(nSV) file
```
awk '{OFS="\t"}{print $1,$2}' DM.chr.fa.fai  > DM.genome

for sample in $(cat list30);do
    for hap in H1 H2;do
        for bin in 500 1000 5000 10000 50000 100000 500000 1000000 5000000;do 
        ####make comlement for repeat region (to speed up)
        bedtools shuffle -i  01_dbin/$sample.$hap.dSV.$bin.st.bed  -g DM.genome | sort -k1,1 -k2,2n > 01_dbin/$sample.$hap.shuf.$bin.st.bed
        done
    done
done
```
- Bin file for dSV and nSV
```
########################### dSV bin ###########################
mkdir 01_dbin
for sample in $(cat /home/chenglin/work/01_sv/01_genome/02_hap/list30);do
    for hap in H1 H2;do
        for bin in 500 1000 5000 10000 50000 100000 500000 1000000 5000000;do 
        awk 'OFS="\t"{print $1,$2-'$bin',$2,$4,$6,"Up"}' dSV/$sample.$hap.dSV.bed | sed -r 's/-[0-9]+/0/g' > 01_dbin/$sample.$hap.dSV.$bin.bed
        awk 'OFS="\t"{print $1,$3,$3+'$bin',$4,$6,"Dn"}' dSV/$sample.$hap.dSV.bed >> 01_dbin/$sample.$hap.dSV.$bin.bed
        sort -k1,1 -k2,2n   01_dbin/$sample.$hap.dSV.$bin.bed >  01_dbin/$sample.$hap.dSV.$bin.st.bed
        rm 01_dbin/$sample.$hap.dSV.$bin.bed
        done
    done
done
########################### nSV bin ###########################
for sample in $(cat /home/chenglin/work/01_sv/01_genome/02_hap/list30);do
    for hap in H1 H2;do
        for bin in 500 1000 5000 10000 50000 100000 500000 1000000;do 
        awk 'OFS="\t"{print $1,$2-'$bin',$2,$4,$6,"Up"}' dSV/$sample.$hap.nSV.bed | sed -r 's/-[0-9]+/0/g' > 01_dbin/$sample.$hap.nSV.$bin.bed
        awk 'OFS="\t"{print $1,$3,$3+'$bin',$4,$6,"Dn"}' dSV/$sample.$hap.nSV.bed >> 01_dbin/$sample.$hap.nSV.$bin.bed
        sort -k1,1 -k2,2n   01_dbin/$sample.$hap.nSV.$bin.bed >  01_dbin/$sample.$hap.nSV.$bin.st.bed
        rm 01_dbin/$sample.$hap.nSV.$bin.bed
        done
    done
done
```
- Overlap dSV nSV with dSNP and SNP
```
###bedtools
mkdir 02_cov
for sample in $(cat list30);do
    for hap in H1 H2;do
        for bin in 500 1000 5000 10000 50000 100000 500000 1000000;do 
            for type in dSV nSV ;do
            bedtools coverage -a  01_dbin/$sample.$hap.$type.$bin.st.bed -b dSNP/$sample.$hap.SNP.bed -wo -nonamecheck > 02_cov/$sample.$hap.$bin.$type.SNP.bed
            bedtools coverage -a  01_dbin/$sample.$hap.$type.$bin.st.bed -b dSNP/$sample.$hap.dSNP.bed -wo -nonamecheck > 02_cov/$sample.$hap.$bin.$type.dSNP.bed
            done
        done
    done
done
```

## For dSV enrichment analysis  of dSNP ( up tp 1Mb)
- Bin file for enrichment analysis
```
mkdir 01_dbin
for sample in $(cat /home/chenglin/work/01_sv/01_genome/02_hap/list30);do
    for hap in H1 H2;do
        for bin in $(seq 100000 100000 1000000) ;do
        awk 'OFS="\t"{print $1,$2-'$bin',$2-'$bin'+100000,$4,$6,"Up"}' dSV/$sample.$hap.dSV.bed | sed -r 's/-[0-9]+/0/g'|sort > 01_dbin/$sample.$hap.dSV.$bin.Up.bed
        awk 'OFS="\t"{print $1,$3+'$bin'-100000,$3+'$bin',$4,$6,"Dn"}' dSV/$sample.$hap.dSV.bed|sort > 01_dbin/$sample.$hap.dSV.$bin.Dn.bed
        done
    done
done
```
- Overlap Bin file with dSNPs
```
mkdir 02_cov
for sample in $(cat /home/chenglin/work/01_sv/01_genome/02_hap/list30);do
    for hap in H1 H2;do
        for bin in $(seq 100000 100000 1000000) ;do
            for type in dSV ;do
            bedtools coverage -a  01_dbin/$sample.$hap.dSV.$bin.Up.bed -b dSNP/$sample.$hap.dSNP.bed -wo -nonamecheck > 02_cov/$sample.$hap.$bin.$type.dSNP.Up.bed
            bedtools coverage -a  01_dbin/$sample.$hap.dSV.$bin.Dn.bed -b dSNP/$sample.$hap.dSNP.bed -wo -nonamecheck > 02_cov/$sample.$hap.$bin.$type.dSNP.Dn.bed
            done
        done
    done
done

```





