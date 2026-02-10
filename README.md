## Analyzing RNA-Seq data from TCGA

Cholangiocarcinoma (CCA) is a rare and aggressive malignancy arising from the bile duct epithelium, characterized by poor prognosis and limited therapeutic options. This project aims to explore the transcriptional landscape of CCA using RNA-seq data from the TCGA-CHOL cohort [1, 2].
As an initial analysis and training exercise, we focus on paired tumor–normal samples to identify tumor-specific gene expression changes while controlling for inter-patient variability. Differential expression analysis is performed using a paired experimental design, followed by Gene Ontology (GO) enrichment analysis to characterize biological processes associated with tumor activation and suppression.
The workflow and statistical framework are adapted from limma-based vignette [3].

### Summary of Results
For intitial practice, difference between paired tumor and normal tissue's gene expression profile will be analyzed by following this vignette:
[https://cran.r-project.org/web/packages/easybio/vignettes/example_limma.html](https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html)

<img width="877" height="442" alt="dotplotGO" src="https://github.com/user-attachments/assets/ee7ed080-3344-4946-92fb-769f80cf7cde" />

GO enrichment analysis of differentially expressed genes highlights distinct biological programs altered in CCA tumors relative to matched normal tissues. 
Activated GO terms were enriched in cell cycle regulation, chromosomal segregation, and mitotic processes, which are hallmarks of cancer in general.
Meanwhile, suppressed GO terms were enriched in metabolism and specifically catabolism processes, corresponding to how cancer cells switch from catabolism to prioritize anabolism to promote growth.


### Reference
[1] The Cancer Genome Atlas Research Network. Comprehensive molecular characterization of human cancers. Nature 490, 61–70 (2012). https://doi.org/10.1038/nature11213

[2] Farshidfar, F. et al. Integrative genomic analysis of cholangiocarcinoma identifies distinct IDH-mutant molecular profiles. Cell Reports 18, 2780–2794 (2017). https://doi.org/10.1016/j.celrep.2017.02.033

[3] UC Davis Bioinformatics Core. (2018). Differential Expression with Limma-Voom. In 2018 June RNA-Seq Workshop (pp. RNA-Seq DE tutorial). UC Davis Bioinformatics Training. Retrieved [date you accessed], from https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html
