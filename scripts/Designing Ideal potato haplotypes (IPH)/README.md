
<img width="600" alt="image" src="https://github.com/Chenglin20170390/Haplotype-diversity/assets/33062118/873b0f8b-86d9-4a47-9ff4-269ca17ad16c">

- support from  @Nan Wang
---
# The general workflow for designing ideal potato haplotypes (IPHs).
- Our aim is to 1) identify best recombination for selecting haplotype blocks containing no deleterious SVs and minimal dSNPs 2) explore the heterosis based on two heterotic groups. Here, we set the priority condition of purging dSVs, and the most ideal haplotype of heterotic A and E groups was the chimeric haplotype with the minimum value of IPHs burden. For our purpose, we dependently the ideal recombination and combination of haplotype blocks for each chromosome in two heterotic groups (hereby A and E groups).
---
## The number of deleterious variants for each chromosome was estimated using the following formula:
<img width="400" alt="image" src="https://github.com/Chenglin20170390/Haplotype-diversity/assets/52204680/9eccfbc8-3c0d-4a07-85e7-588b16b0825a">

Where n represents the recombination fragments.
The parameter j indicates recombination events in each chromosome (50 bins). 
The P indicates the number of dSVs (dSV burden) and p indicates the number of dSNPs (dSNP burden). 
The most ideal haplotype of heterotic A and E groups was the chimeric haplotype with the minimum value of IPHs burden.

- Here are two important parameters should be considered in this study:
1) j indicates recombination events in each haplotype (here are 4)
2) the number of bins for each chromosome (here are 50)
---
Preparation deleterious variation map from all samples (for example C005 H1)
A bed formation of dSNPs file, split by tab
```
chr01	10007440	10007441
chr01	10009811	10009812
chr01	10432127	10432128
chr01	10453615	10453616
chr01	10476270	10476271
chr01	10545151	10545152
chr01	10820244	10820245
chr01	10820812	10820813
chr01	10820819	10820820
chr01	10820839	10820840
chr01	10821616	10821617
```
## IPH calculate process
- 1. distigushed the heterotic groups based on haplotype phylogeny
```
mashtree --mindepth 0 --numcpus 12 *.fasta > potato_mashtree.nwk
```
- 2. for example, design the minimal dSVs in heterotic A group
```
Python Design_less_dSVs_in_A_line_IPH.py -chr chr01 -o combinations chr01_A_minimal_dSVs.txt
```
The output like below, while the potential minimal dSVs recombination were collected for further dSNPs estimation.
```
dSVs	breakpoint1 	breakpoint2 	breakpoint3 	breakpoint4 	fragment1 	fragment2 	fragment3 	fragment4 	fragment5 	
0	49611344.16	54926845.32	65557847.64	77960683.68	C151.H1	C005.H1_C025.H1_C056.H1_C091.H1_C093.H1_C151.H1_C098.H1_C121.H1_C174.H1_C091.H2_C093.H2_C151.H2_C121.H2_C174.H2	C005.H1_C091.H1_C121.H1_C025.H2_C098.H2	AE.H1	C174.H1_C121.H2
0	46067676.72	53155011.6	63786013.92	76188849.96	C151.H1	C005.H1_C056.H1_C091.H1_C093.H1_C151.H1_C098.H1_C121.H1_C174.H1_C025.H2_C093.H2_C151.H2_C121.H2	C005.H1_C025.H1_C091.H1_C093.H1_C121.H1_C091.H2_C121.H2	AE.H1	C174.H1
4	1771833.72	3543667.44	5315501.16	7087334.88	C025.H1_C056.H1_C091.H1_C151.H1_C121.H1_C025.H2_C091.H2_C151.H2_C098.H2_C174.H2_AE.H1	C056.H1_C091.H1_C151.H1_C121.H1_C056.H2_C091.H2_C098.H2_C121.H2_C174.H2	C091.H1_C151.H1_C121.H1_C025.H2	C005.H1_C091.H1_C151.H1_C121.H1_C174.H1_C025.H2_C056.H2_C093.H2_C174.H2	C151.H1_C121.H2
4	1771833.72	3543667.44	5315501.16	8859168.6	C025.H1_C056.H1_C091.H1_C151.H1_C121.H1_C025.H2_C091.H2_C151.H2_C098.H2_C174.H2_AE.H1	C056.H1_C091.H1_C151.H1_C121.H1_C056.H2_C091.H2_C098.H2_C121.H2_C174.H2	C091.H1_C151.H1_C121.H1_C025.H2	C005.H1_C091.H1_C151.H1_C174.H1_C025.H2_C174.H2	C151.H1_C121.H2
```
- 3. design the minimal dSNPs in heterotic E group, this script will calculate the lowest dSNPs based on the combinations of recombination minimal dSVs datasets. Our design considers the low recombination rate of inversions (INVs) and centromeres for each chromosome.
```
Python Design_less_dSNPs_in_two_IPH_lines.py -chr chr01 -file combinations chr01_A_minimald_SVs.txt -o combinations chr01_A_minimal_dSVs_dSNPs.txt
```
The output like below, while the potential minimal dSVs and dSNPs recombination were identified as the most ideal potato haplotype for chromosome 1 of heterotic A group.
```
CHR	dSVs	dSNPs	breakpoint1 	breakpoint2 	breakpoint3 	breakpoint4 	fragment1 	fragment2 	fragment3 	fragment4 	fragment5 	
chr01	0	5529	44295843.0	65557847.64	69101515.08	76188849.96	C151.H1	C121.H1	C025.H1_C091.H2	AE.H1	C174.H1
chr01	0	5529	42524009.28	65557847.64	69101515.08	76188849.96	C151.H1	C121.H1	C025.H1_C091.H2	AE.H1	C174.H1
chr01	0	5536	31893006.96	65557847.64	69101515.08	76188849.96	C151.H1	C121.H1	C025.H1_C091.H2	AE.H1	C174.H1
chr01	0	5539	46067676.72	65557847.64	69101515.08	76188849.96	C151.H1	C121.H1	C025.H1_C091.H2	AE.H1	C174.H1
chr01	0	5543	42524009.28	65557847.64	69101515.08	77960683.67999999	C151.H1	C121.H1	C025.H1_C091.H2	AE.H1	C174.H1
chr01	0	5543	44295843.0	65557847.64	69101515.08	77960683.67999999	C151.H1	C121.H1	C025.H1_C091.H2	AE.H1	C174.H1
chr01	0	5545	47839510.44	65557847.64	69101515.08	76188849.96	C151.H1	C121.H1	C025.H1_C091.H2	AE.H1	C174.H1
chr01	0	5547	33664840.68	65557847.64	69101515.08	76188849.96	C151.H1	C121.H1	C025.H1_C091.H2	AE.H1	C174.H1
chr01	0	5550	40752175.56	65557847.64	69101515.08	76188849.96	C151.H1	C121.H1	C025.H1_C091.H2	AE.H1	C174.H1
chr01	0	5550	37208508.12	65557847.64	69101515.08	76188849.96	C151.H1	C121.H1	C025.H1_C091.H2	AE.H1	C174.H1
```
