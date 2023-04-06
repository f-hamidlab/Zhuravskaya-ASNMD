## Load dependencies for this workflow
library(tidyverse)
library(DESeq2)
library(reticulate)

# set variables
wd <- "/media/cdn-bc/RAID/Projects/FHyyy_Anya_ES/gitrepo/data"
pco <- 0.05

## Set working directory
setwd(wd)

## load DMSO-only DDS object
dds.lrt <- readRDS("rds/DDS_DMSO_only.rds")

## Stabilize variance of counts
vsd <- assay(vst(dds.lrt, blind=FALSE))
colnames(vsd) <- paste("VST", dds.lrt$day, dds.lrt$batch, sep = "_")


## Prepare matrix for DP_GP_cluster.py

### Get CDS-containing genes with significant LRT change
results.lrt <- read_tsv("tables/lrt_CDSs_DMSO_only.tsv")
sig.lrt <- results.lrt %>% 
    as.data.frame() %>% 
    filter(padj < pco) %>% 
    pull(gene_id)

# get vst-transformed expression of significant LRT genes
vsd.sig <- vsd[sig.lrt,]

# Average expression from replicates
vsd.sig.averaged <- t(rowsum(t(vsd.sig), rep(0:4, each=2)))/2

# output dataframe for DP_GP_cluster.py
vsd.sig.averaged %>% 
    as.data.frame() %>% 
    rownames_to_column("gene_id") %>% 
    write_tsv("matrices/vst_sig_cds_dmso_averaged.tsv")


use_condaenv("py2.7")

fn <- "DP_GP_cluster.py -i matrices/vst_sig_cds_dmso_averaged.tsv
-o DP_GP/vst_sig_cds --fast -n 100 --max_iters 100 
--check_burnin_convergence --check_convergence 
--cluster_uncertainty_estimate --plot -p pdf"


system2(str_remove_all(fn, "\n"))