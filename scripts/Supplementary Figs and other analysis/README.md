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

### - haplotype consensus tree and density tree
- haplotype consensus tree (-t: threads; -M output multiple sequence alignment format for orthologous genes)
```
##prepare input sequences for iqtree
orthofinder -t 46 -fg /home/chenglin/work/01_sv/01_genome/04_analysis/08_orthfinder/03_pep/OrthoFinder/Results_Oct12_1 -M msa 

##running tree with best module (if no option -m, software will estimated suitable module for your tree)
iqtree2 -s $dir/MultipleSequenceAlignments/SpeciesTreeAlignment.fa  -m  JTT+F+R5 -o Outgroup -T 48 -B 1000
```

- split chromosome tree
```
ls ../Single_Copy_Orthologue_Sequences |wc -l ## check how many synteny genes from orthfinder

##make single copy gene list 
ls ../Single_Copy_Orthologue_Sequences/ | sed 's/\.fa//g' | while read line;do
# ../Orthogroups/Orthogroups.tsv 
#echo $line 
grep -w $line ../Orthogroups/Orthogroups.tsv  >> single.copy.gene.list
done

##make a reference gene ID(DM)for determination of chromosome level
mkdir 01_chr
awk '{print $1,$58}' single.copy.gene.list > single.copy.gene.DM.list
for i in {01..12};do
grep 'DM.'$i''  single.copy.gene.DM.list > 01_chr/DM.chr$i.list
done

###make tree for each chromosome
mkdir 02_og_rename_se
for i in {01..12};do
  echo '
  i='$i'
  rm 02_og_rename_se/chr$i.rename.list
  cat 01_chr/DM.chr$i.list| while read line;do
  og=`echo $line |awk '"'"'{print $1}'"'"'`
  echo 02_og_rename_se/$og.rename.fa >> 02_og_rename_se/chr$i.rename.list
  #cat MultipleSequenceAlignments/$og.fa | sed '"'"'s/_.\+//g'"'"' > 02_og_rename_se/$og.rename.fa
  done
  ' > chr$i.sh
  nohup sh chr$i.sh &
done

###merge multiple fasta with same in  one file
mkdir 03_chr_fa
for i in {01..12};do
nohup seqkit concat --infile-list 02_og_rename_se/chr$i.rename.list -o 03_chr_fa/chr$i.pep.aligned.fa &
done

##running iqtree for each tree 
for i in {02..12};do
nohup iqtree2 -s  03_chr_fa/chr$i.pep.aligned.fa  -m  JTT+F+R5 -o Outgroup -T 4 -B 1000 > chr$i.log &
done
i=01
nohup iqtree2 -s  03_chr_fa/chr$i.pep.aligned.fa  -m  JTT+F+R5 -o Outgroup -T 10 -B 1000 > chr$i.log  &

##then you can visulization by the function of ggdensitree from ggtree(https://yulab-smu.top/treedata-book/chapter4.html)
```
