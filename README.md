# Ideal Potato Haplotyes (IPHs v1.0)![github](https://img.shields.io/badge/3C-Certification-red)        
![github](https://img.shields.io/badge/Potato-square-hex)        ![github](https://img.shields.io/badge/Haplotype--resolved-green)         ![github](https://img.shields.io/badge/Deleterious-SV-red)       


IPHs v1.0  was designed form two heterotic groups with minimized dSVs and dSNPs, which will guide recurrent selection of elite inbred lines (The Upotato Plan:https://agis.caas.cn/en/research/principalinvestigator/237658.htm).

In the future, we will develop IPH v2.0, incorporating of beneficial alleles associated with vigorous agronomic traits and further reducing deleterious variations through synthetic biology and genome editing. Although this step is currently challenging, we believe it will be possible to artificially synthesize potato chromosomes, even genomes, to create desirable ideal potato varieties within the next few decades.
<img width="790" alt="image" src="https://github.com/Chenglin20170390/Haplotype-diversity/assets/33062118/9fc55715-d1a4-4a34-8c2e-1344a3913151">

# Workflow Overview
<img width="600" alt="image" src="https://github.com/Chenglin20170390/Haplotype-diversity/assets/33062118/97935282-4ee3-4084-a70c-0ef8f6f94134">

The scripts expect to be run in roughly this order:

- [1. Genome survey, assembly, accessment and annotation.](https://github.com/Chenglin20170390/Haplotype-diversity/tree/main/scripts/Genome%20assembly%20and%20annotation)
- [2. Graph construction by PGGB and MC](https://github.com/Chenglin20170390/Haplotype-diversity/tree/main/scripts/Graph%20construction%20by%20PGGB%20and%20MC)
- [3. SV calling based on genome alignments](https://github.com/Chenglin20170390/Haplotype-diversity/tree/main/scripts/SV%20calling%20based%20on%20genome%20alignments)
- [4. TE mediated SV formation](https://github.com/Chenglin20170390/Haplotype-diversity/tree/main/scripts/TE%20mediated%20SV%20formation)
- [5. Identification of Deleterious SV](https://github.com/Chenglin20170390/Haplotype-diversity/tree/main/scripts/Identification%20of%20Deleterious%20SV)
- [6. Correlation between dSV and dSNP](https://github.com/Chenglin20170390/Haplotype-diversity/tree/main/scripts/Correlation%20between%20dSV%20and%20dSNP)
- [7. Designing Ideal potato haplotypes (IPH)](https://github.com/Chenglin20170390/Haplotype-diversity/tree/main/scripts/Designing%20Ideal%20potato%20haplotypes%20(IPH))![github](https://img.shields.io/badge/3C-Certification-red)  


# Citations
Please cite our paper if you find scripts useful:
***
AA, BB, CC. et al. Leveraging a phased pangenome to design ideal haplotypes for hybrid potato breeding. article 1–7 (2023/2024).



# Acknowledgements
I want to thank [Nan Wang](https://github.com/wangnan9394), [Zhigui Bao](https://github.com/baozg) for contribution of the script for this analysis. I also want to thank [Peter L. Morrell](https://github.com/pmorrell), [ZhiYang Zhang](https://github.com/zhangzhiyangcs) and many more others for suggesting, testing, debugging, and improving the pipeline.





