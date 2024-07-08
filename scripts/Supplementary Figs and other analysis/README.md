
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

```
- then you can visulization by the function of ggdensitree from ggtree(https://yulab-smu.top/treedata-book/chapter4.html)

### - Principal component analysis(PCA) for haplotype diversity
```
plink --vcf all.vcf.gz --recode --out 60hap --allow-extra-chr --pca  5
```

