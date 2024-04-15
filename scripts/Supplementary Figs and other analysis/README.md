## Comparing heterozygosity based on SNPs, k-mer and genome variants by sequence length.

- 1VCF files with SNPs and SVs (insertion, deletion, inversion and duplication) from syri.output

```
for sample in $(cat /home/chenglin/work/01_sv/01_genome/02_hap/list30);do
echo $sample 
grep -E '#|SNP|INS|DEL|INV|DUP|' 02_syri/$sample.H1.$sample.H2/$sample.H1.$sample.H2.syri.out | grep -v -E '<HDR>|<INVTR>|<INVTRAL>|<DUPAL>|<INVAL>|INVDPAL' |grep -v -w -E "INVTR|INVTRAL|INVAL|DUPAL" > 03_hapout/$sample.H1.$sample.H2.ft.out 
bedtools merge -i <(sort -k1,1 -k2,2n 03_hapout/$sample.H1.$sample.H2.ft.out)  > 04_haprmdup/$sample.H1.$sample.H2.rmdup.bed 
done
```

