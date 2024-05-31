# Alternative splicing coupled to nonsense-mediated decay coordinates downregulation of non-neuronal genes in developing mouse neurons

This repository contains scripts used to generate the data
presented in our manuscript. The codes for processing RNA-sequencing data
for gene expression and alternative splicing analyses can be found in the `scripts`
directory. A `.qmd` file with codes to reproduce relevant figures from the manuscript
can be found in the `vignette` directory.

## Abstract

**Background**  
The functional coupling between alternative pre-mRNA splicing (AS) and the mRNA quality control mechanism called nonsense-mediated decay (NMD) can modulate transcript abundance. Previous studies have identified several examples of such a regulation in developing neurons. However, the systems-level effects of AS-NMD in this context are poorly understood.

**Results **  
We developed an R package, factR2, which offers a comprehensive suite of AS-NMD analysis functions. Using this tool, we conducted a longitudinal analysis of gene expression in pluripotent stem cells undergoing induced neuronal differentiation. Our analysis uncovered hundreds of AS-NMD events with significant potential to regulate gene expression. Notably, this regulation was significantly overrepresented in specific functional groups of developmentally downregulated genes. Particularly strong association with gene downregulation was detected for alternative cassette exons stimulating NMD (NS-CEs) upon their inclusion into mature mRNA. By combining bioinformatics analyses with CRISPR/Cas9 genome editing and other experimental approaches we show that NS-CEs regulated by the RNA-binding protein PTBP1 dampen the expression of their genes in developing neurons. We also provide evidence that the NS-CE activity is temporally coordinated with NMD-independent gene repression mechanisms.

**Conclusions **  
Our study provides an accessible workflow for the discovery and prioritization of AS-NMD targets. It further argues that the AS-NMD pathway plays a widespread role in developing neurons by facilitating the downregulation of functionally related non-neuronal genes. 

Please cite:
{Insert manuscript/publication citation here}
