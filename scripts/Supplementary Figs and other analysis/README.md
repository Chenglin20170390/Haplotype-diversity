### - Comparing haplotype divergence (heterozygosity) based on SNPs, k-mer and genome variants by sequence length(GVSL).

- VCF files with SNPs and SVs (insertion, deletion, inversion and duplication) from syri.output

```
for sample in $(cat /home/chenglin/work/01_sv/01_genome/02_hap/list30);do
echo $sample 
grep -E '#|SNP|INS|DEL|INV|DUP|' 02_syri/$sample.H1.$sample.H2/$sample.H1.$sample.H2.syri.out | grep -v -E '<HDR>|<INVTR>|<INVTRAL>|<DUPAL>|<INVAL>|INVDPAL' |grep -v -w -E "INVTR|INVTRAL|INVAL|DUPAL" > 03_hapout/$sample.H1.$sample.H2.ft.out 
bedtools merge -i <(sort -k1,1 -k2,2n 03_hapout/$sample.H1.$sample.H2.ft.out)  > 04_haprmdup/$sample.H1.$sample.H2.rmdup.bed 
done
```
- Statistic for comparison of haplotype divergence.
```
for sample in $(cat /home/chenglin/work/01_sv/01_genome/02_hap/list30);do
snp_size=`grep -w 'SNP' 03_hapout/$sample.H1.$sample.H2.ft.out |wc -l `
indels_s1=`grep -w 'INS' 03_hapout/$sample.H1.$sample.H2.ft.out |awk '$8-$7<50'|awk '{print $8-$7}' |sed 's/-//g'  |awk '{sum+=$1}END{print sum}'  `
indels_s2=`grep -w 'DEL' 03_hapout/$sample.H1.$sample.H2.ft.out | awk '$3-$2<50' |awk '{sum+=$3-$2} END{print sum}' `
ins_size=`grep -w 'INS' 03_hapout/$sample.H1.$sample.H2.ft.out |awk '$8-$7>49'|awk '{print $8-$7}' |sed 's/-//g'  |awk '{sum+=$1}END{print sum}'  `
del_size=`grep -w 'DEL' 03_hapout/$sample.H1.$sample.H2.ft.out | awk '$3-$2>49'|awk '{sum+=$3-$2}END{print sum}' `
inv_size=`grep -w 'INV' 03_hapout/$sample.H1.$sample.H2.ft.out | awk '{sum+=$3-$2}END{print sum}' `
dup_size=`grep -w 'DUP' 03_hapout/$sample.H1.$sample.H2.ft.out | awk '{sum+=$3-$2}END{print sum}' `
indels=`expr $indels_s1 + $indels_s2`
v_size=`expr $snp_size + $ins_size + $del_size + $inv_size + $dup_size`
asm_size=`awk '{sum+=$2}END{print sum/2}' $asmdir/$sample.H1.chr.fa.fai $asmdir/$sample.H2.chr.fa.fai`
asm_per=`awk 'BEGIN{printf "%.2f%\n",('$v_size'/'$asm_size')*100}'`
echo $sample $snp_size $indels $ins_size  $del_size  $inv_size  $dup_size $asm_size $asm_per |sed 's/\s\+/\t/g' >> divergence.all.txt
done

sed  -i '1iID\tSNP\tInDels\tINS\tDEL\tINV\tDUP\tGenome_size\tVariant_percent' divergence.all.txt
```

### - ABBA-BABA statistics for introgression
- The software ABBABABAwindows.py from https://github.com/simonhmartin/genomics_general was used for introgesssion.
```
python ABBABABAwindows.py -g parse.snp.gt -f phased  -o $sample.abba.csv  -P1 PHU  \
-P2 STN   -P3 GON  -O CND  --popsFile pop.C027.txt   \
--minData 0.5 -w 1000000 -s 500000 -m 2  -T 30
```
