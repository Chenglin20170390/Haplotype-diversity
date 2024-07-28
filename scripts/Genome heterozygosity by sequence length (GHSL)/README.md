## - Genome heterozygosity by sequence length(GHSL) 
<img width="730" alt="image" src="https://github.com/user-attachments/assets/a37c1e9a-7c6b-424d-919e-4fe3d659846d">


- Step
  1) The phased genome files, for example: chr1_Hap1, chr1_Hap2.
  2) Calling variantions (SNPs, SVs)
  3) GHSL caculation (the length ofdivergent sequence between two haplotypes)


- For example: VCF calling from Syri with SNPs and SVs (insertion, deletion, inversion and duplication)

```
## filtering variations
for sample in $(cat list30);do
echo $sample 
grep -E '#|SNP|INS|DEL|INV|DUP|' $sample.H1.$sample.H2.syri.out | grep -v -E '<HDR>|<INVTR>|<INVTRAL>|<DUPAL>|<INVAL>|INVDPAL' |grep -v -w -E "INVTR|INVTRAL|INVAL|DUPAL" > 03_hapout/$sample.H1.$sample.H2.ft.out
done
```
- Statistic for GHSL
```
for sample in $(cat list30);do
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

##Note:: for Non-redundant lenght you may used 'bedtools merge' commands.##
```

