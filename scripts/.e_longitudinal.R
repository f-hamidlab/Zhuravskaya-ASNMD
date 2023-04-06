library(tidyverse)
library(DESeq2)

# set variables
wd <- "/media/cdn-bc/RAID/Projects/FHyyy_Anya_ES/gitrepo/data"
pco <-  1e-3

## Set working directory
setwd(wd)

# Load Annotations/Metadata
tx2gene_meta <- read.delim("tables/transcripts_info_final_transcriptome.tsv.gz", 
                           sep="\t", header=T)
# Droping unneeded columns and renaming remaining columns:
tx2gene_meta <- tx2gene_meta %>% 
    select(chr=seqnames, beginning = start, ending = end, nt = width,
           str = strand, source, gene_id, gene_name,gene_type, level = match_level,
           mgi_id, havana_gene, transcript_id)


nmd <- read.delim("tables/transcripts_final_transcriptome_factR2.tsv.gz", 
                  sep="\t", header=T) %>% 
    dplyr::select(transcript_id, stop_to_lastEJ, num_of_downEJs, 
                  UTR3_length = X3.UTR_length, is_NMD)

# Combine txinfo with NMD info
tx2gene_meta = left_join(tx2gene_meta, nmd, by = "transcript_id")

# load DESeq2DataSet object containing gene expression data of all samples
dds <- readRDS("rds/DDS_all_samples.rds")



# subset dds for DMSO samples only
# Also changing the design formula to ~ day
# Note that I am not doing ~ batch + day since
# "batches" are just technical replicates and are very tight
dds.lrt <- dds[,dds$condition == "DMSO"]
dds.lrt$condition = droplevels(dds.lrt$condition)
design(dds.lrt) = ~ day 


dds.lrt = estimateSizeFactors(dds.lrt)
idx     = rowSums( counts(dds.lrt, normalized=TRUE) >= 5 ) >= 2 # Filter so that at least 2 samples have count>=5
dds.lrt = dds.lrt[idx,]


dds.lrt <- DESeq(dds.lrt, test="LRT", reduced =  ~ 1)
results.lrt <- results(dds.lrt)


# Make time-ordered VST-transformed gene expession table containing DMSO data
vsd <- assay(vst(dds.lrt, blind=FALSE))
colnames(vsd) <- paste("VST", dds.lrt$day, dds.lrt$batch, sep = "_")

KeT = apply(vsd, 1, function(x) cor.test(rep(1:5, each=2), x, method ="kendall")$estimate[[1]])
KeP = apply(vsd, 1, function(x) cor.test(rep(1:5, each=2), x, method ="kendall")$p.value[[1]])

# prepare data frame containing coding genes showing significant LRT for DP_GP_cluster.py

cdss <- tx2gene_meta %>% 
    filter(!is.na(is_NMD)) %>% 
    pull(gene_id) %>% 
    unique()

sig.cdss <- results.lrt %>% 
    as.data.frame() %>% 
    filter(padj < pco) %>% 
    rownames_to_column("Gene") %>% 
    filter(Gene %in% cdss) %>% 
    pull(Gene)

# get vst-transformed expression of significant LRT genes
vsd.sig <- vsd[sig.cdss,]

# Average expression from replicates
vsd.sig.averaged <- t(rowsum(t(vsd.sig), rep(0:4, each=2)))/2

# output dataframe for DP_GP_cluster.py
vsd.sig.averaged %>% 
    as.data.frame() %>% 
    rownames_to_column("gene_id") %>% 
    write_tsv("tables/DP_GP_ALLsig.vst.tsv")

# save DDS object
saveRDS(dds.lrt, "rds/DDS_DMSOonly.rds")

results.lrt$KeT <- KeT
results.lrt$KeP <- KeP
results.lrt[sig.cdss,] %>% 
    as.data.frame() %>% 
    rownames_to_column("gene_id") %>% 
    write_tsv("tables/ALLsig_LRT.tsv")































