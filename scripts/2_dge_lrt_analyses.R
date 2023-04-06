## Load dependencies for this workflow
library(tidyverse)
library(tximport)
library(DESeq2)

## set variables

### Working directory containing input data
wd <- "/media/cdn-bc/RAID/Projects/FHyyy_Anya_ES/gitrepo/data" 
pco <- 0.05

## Set working directory
setwd(wd)

## Load transcript metadata
tx.meta <- read_tsv("metadata/txs_metadata.tsv.gz")


## Create vector containing path to kallisto files
sample_id <-  dir("kallisto")
kal_files <-  file.path("kallisto", sample_id, "abundance.h5")
names(kal_files) <- sample_id

## Load sample metadata
s2c_master <- read_tsv("metadata/samples.txt")
s2c_master <- s2c_master %>%
    as.data.frame() %>% 
    mutate(group = paste0(condition,day)) %>%
    mutate(day = paste0("Day", day)) %>% 
    mutate(across(-sample, .fns = as.factor))
rownames(s2c_master) <- s2c_master$sample


## Import Kallisto transcript expression and collapse to gene-level expression
kal_files <- kal_files[s2c_master$sample]
txi.kallisto <- tximport(kal_files, type = "kallisto", 
                         tx2gene = tx.meta[,c("transcript_id", "gene_id")], 
                         ignoreAfterBar = F)

# create DESeq2DataSet object
dds.gn  <- DESeqDataSetFromTximport(txi.kallisto, s2c_master, ~group)


## Subsetting dds object for CHX vs DMSO Wald-test DESeq2 analysis
### - Iterate each timepoint
### - Relevel groups and days
### - Reset design

timepoints <- levels(dds.gn$day)
names(timepoints) <- paste0("Day",timepoints)

ddslist.gn <- lapply(timepoints, function(t){
    dds.tmp <- dds.gn[,dds.gn$day == t]
    dds.tmp$group <- droplevels(dds.tmp$group)
    dds.tmp$condition <- relevel(dds.tmp$condition, ref = "DMSO")
    design(dds.tmp) = ~condition
    dds.tmp
})


## Filter genes by timepoint
### So that at least 2 samples have count>=5

ddslist.gn <- lapply(ddslist.gn, function(dds){
    dds <- estimateSizeFactors(dds)
    keep <- rowSums(counts(dds, normalized = TRUE) >=5 ) >= 2
    dds[keep,]
})

## Run DESeq for each dds object
ddslist.gn <- lapply(ddslist.gn, DESeq, fitType = "local")

## Get results from CHX vs DMSO test and compile results
deg.compiled <- do.call(bind_rows,lapply(ddslist.gn, function(dds){
    results(dds) %>% 
        as.data.frame() %>% 
        rownames_to_column("gene_id") %>% 
        arrange(desc(padj)) %>% 
        mutate(Timepoint = unique(dds$day))
}))

## Run LRT test on DMSO samples only
### subset dds for DMSO samples only
### Also changing the design formula to ~ day
### Note that I am not doing ~ batch + day since
### "batches" are just technical replicates and are very tight
dds.lrt <- dds.gn[,dds.gn$condition == "DMSO"]
dds.lrt$condition = droplevels(dds.lrt$condition)
design(dds.lrt) = ~ day 
dds.lrt <- estimateSizeFactors(dds.lrt)

## Filter genes so that at least 2 samples have count>=5
idx     = rowSums( counts(dds.lrt, normalized=TRUE) >= 5 ) >= 2 
dds.lrt = dds.lrt[idx,]

## Run DESeq and get LRT results
dds.lrt <- DESeq(dds.lrt, test="LRT", reduced =  ~ 1)
lrt <- results(dds.lrt) %>% 
    as.data.frame() %>% 
    rownames_to_column("gene_id")


## Subset for CDS-containing genes
cds.gns <- tx.meta %>% 
    filter(!is.na(is_NMD)) %>% 
    pull(gene_id) %>% 
    unique()

## Filter for CDS genes
deg.cds.compiled <- deg.compiled %>% 
    filter(gene_id %in% cds.gns)
lrt.cds <- lrt %>% 
    filter(gene_id %in% cds.gns)


## Save DDS objects
saveRDS(dds.gn, "rds/DDS_all_samples.rds")
saveRDS(ddslist.gn, "rds/DDS_list_by_timepoint.rds")
saveRDS(dds.lrt, "rds/DDS_DMSO_only.rds")

## Save results of comparative analyses
write_tsv(deg.cds.compiled, "tables/dge_CDSs_CHX_vs_DMSO_bytimepoint.tsv")
write_tsv(lrt.cds, "tables/lrt_CDSs_DMSO_only.tsv")

## Prepare matrix of variance-stabilised counts

### Stabilize variance of counts
vsd <- assay(vst(dds.lrt, blind=FALSE))
colnames(vsd) <- paste("VST", dds.lrt$day, dds.lrt$batch, sep = "_")


### Get CDS-containing genes with significant LRT change
sig.lrt <- lrt %>% 
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






































