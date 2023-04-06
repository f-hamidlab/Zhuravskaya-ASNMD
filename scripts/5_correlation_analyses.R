## Load dependencies for this workflow
library(tidyverse)
library(DESeq2)
library(vroom)
devtools::load_all("/home/cdn-bc/Github_repo/factR2/")

## Set variables
### Working directory containing input data
wd <- "/media/cdn-bc/RAID/Projects/FHyyy_Anya_ES/gitrepo/data" 
co_KeT <- 0.75
co_KeP <-  0.05


## Set working directory
setwd(wd)

## Load cds_event_time metadata
cds.event.time <- read_tsv("metadata/cds_event_time_meta.tsv")

## Correlation 1: Kendall correlation of VST-transformed GEX per DPGP cluster

### load vst-transformed, LRT-significant genes
vst.lrt <- read_tsv("matrices/vst_sig_cds_dmso_averaged.tsv") %>% 
    column_to_rownames("gene_id") %>% 
    as.matrix()

### load DPGP output
gn.clusters <- read_tsv("/media/cdn-bc/RAID/Projects/FHyyy_Anya_ES/originals/Anya_AS_NMD/DESeq2_results/GP_DPGP_ALLsig.VST_output_1k/DP_GP_ALLsig.VST_optimal_clustering.txt")
gn.clusters <- gn.clusters %>% 
    as.data.frame() %>% 
    column_to_rownames("gene")

### Scale matrix
vst.lrt <-  vst.lrt - apply(vst.lrt, 1, mean) # center
vst.lrt <-  vst.lrt / apply(vst.lrt, 1, sd)   # scale
vst.lrt.clust <- vst.lrt %>% 
    as.data.frame() %>% 
    mutate(cluster = gn.clusters[rownames(.),]$cluster) %>% 
    filter(!is.na(cluster))

### Calculate Kendall TAUs for clusters
vst.lrt.clust.kendall <- do.call(bind_rows,lapply(unique(vst.lrt.clust$cluster), 
                                                  function(clust){
                                                      clust.genes <- rownames(gn.clusters[gn.clusters$cluster==clust,])
                                                      dat <- vst.lrt[clust.genes,] %>% 
                                                          as.data.frame() %>% 
                                                          rownames_to_column("gene") %>% 
                                                          pivot_longer(cols = -gene, values_to = "exp", names_to = "time") %>% 
                                                          mutate(time = as.numeric(time))
                                                      kendall.out <- cor.test(dat$exp, dat$time, method = "kendall")
                                                      data.frame(cluster = clust,
                                                                 n_genes = length(clust.genes),
                                                                 cluster_KeT = kendall.out$estimate,
                                                                 cluster_KeP = kendall.out$p.value)
                                                  }))
rownames(vst.lrt.clust.kendall) <- NULL

### Rank clusters by increasing KeT and by decreasing KeP
vst.lrt.clust.kendall.ranked <- vst.lrt.clust.kendall %>% 
    mutate(rank_by_KeT = rank(cluster_KeT),
           rank_by_KeP = rank(cluster_KeP, ties.method = "random")) %>% 
    arrange(rank_by_KeT) %>% 
    mutate(trend = "CX") %>% 
    mutate(trend = ifelse(cluster_KeT < -0.5 & cluster_KeP < 0.05, "DN", trend)) %>% 
    mutate(trend = ifelse(cluster_KeT > 0.5 & cluster_KeP < 0.05, "UP", trend)) 


### Annotated rank_by_KeT clusters into cds_event_time metadata
gn.clusters <- gn.clusters %>% 
    rownames_to_column("gene_id") %>% 
    left_join(vst.lrt.clust.kendall.ranked)

cds.event.time <- cds.event.time %>% 
    left_join(gn.clusters %>% dplyr::select(gene_id, rank_by_KeT))




### Output ranked_clusters 
write_tsv(vst.lrt.clust.kendall.ranked, 
          "tables/DP_GP_clusters_stats_ranked.tsv")

## Correlation 2: Kendall's correlation of VST-transformed GEX across timepoint

### Load GEX data
dds <- readRDS("rds/DDS_DMSO_only.rds")
dds <- estimateSizeFactors(dds)

### VST-transform counts
vst <- assay(vst(dds))

### Run kendall correlation
timepoint <- rep(1:5, each = 2)
kendall.vst <- do.call(bind_rows,apply(vst, 1, function(x){
    test <- cor.test(timepoint, x,  method = "kendall")
    data.frame(gex.KeT = test$estimate, gex.KeP = test$p.value)
}))
rownames(kendall.vst) <- rownames(vst)
kendall.vst %>% 
    mutate(gene_id = rownames(vst)) %>% 
    write_tsv("tables/gex_kendall.tsv.gz")

### Annotate direction of kendall trend into cds_event_time metadata
kendall.vst <- kendall.vst %>% 
    rownames_to_column("gene_id") %>% 
    mutate(gex.Kend.dir = "Not-regulated") %>% 
    mutate(gex.Kend.dir = ifelse(gex.KeT < -co_KeT & gex.KeP < co_KeP, "Down", gex.Kend.dir)) %>% 
    mutate(gex.Kend.dir = ifelse(gex.KeT > co_KeT & gex.KeP < co_KeP, "Up", gex.Kend.dir)) 

cds.event.time <- cds.event.time %>% 
    left_join(kendall.vst %>% dplyr::select(gene_id, gex.Kend.dir))


## Correlation 3: Kendall's correlation of VST-transformed PSI across timepoint

### Import samples metadata
s2c_master <- read.table("metadata/samples.txt", header = TRUE)
s2c_dmso <- s2c_master %>% filter(condition == "DMSO")

### Import Whippet PSI data
whippet.outs <- file.path("whippet/PSI", paste0(s2c_dmso$sample,".psi.gz") )
whippet.psis <- vroom(whippet.outs, 
                      id = "Sample")
whippet.psis$Sample <- str_remove(whippet.psis$Sample,
                                  "whippet/PSI/")
whippet.psis$Sample <- str_remove(whippet.psis$Sample,
                                  ".psi.gz")

## Rename sample to be descriptive
s2c_dmso <- s2c_dmso %>% 
    mutate(name = str_glue("dmso_D{str_pad(day,2,pad=0)}_{batch}"))
whippet.psis <- whippet.psis %>% 
    left_join(s2c_dmso %>% dplyr::select(Sample=sample,Name=name))

whippet.psis <- whippet.psis %>% 
    mutate(event.id = paste0(Coord, Strand, Gene))

## import AS-NMD data
asnmd <- read_tsv("tables/ds_iNs_CHX_vs_DMSO_bytimepoint.tsv")

## Simplify asnmd dataframe and append PSI values from each sample
asnmd <- asnmd %>% 
    dplyr::select(Gene, event.id, Type, ASNMDtype, ASNMD.in.cds) %>% 
    filter(!is.na(event.id))
whippet.psi.pivot <- whippet.psis %>% 
    dplyr::select(event.id, Name, Psi) %>% 
    filter(event.id %in% asnmd$event.id) %>% 
    pivot_wider(names_from = Name, values_from = Psi) %>% 
    column_to_rownames("event.id")

### Run Kendall correlation on PSI
asinTransform <-  function(x) { 2 * asin(sqrt(x))/pi }
whippet.vst.pivot <- asinTransform(as.matrix(whippet.psi.pivot))
whippet.vst.pivot <- whippet.vst.pivot[rowSums(is.na(whippet.vst.pivot)) <= 5, ]

timepoint <- rep(1:5, each = 2)
kendall.psi <- do.call(bind_rows,apply(whippet.vst.pivot, 1, function(x){
    test <- cor.test( timepoint[!is.na(x)],x[!is.na(x)], method = "kendall")
    data.frame(psi.KeT = test$estimate, psi.KeP = test$p.value)
}))

rownames(kendall.psi) <- rownames(whippet.vst.pivot)
kendall.psi %>% 
    rownames_to_column("event.id") %>% 
    write_tsv("tables/exons_kendall.tsv.gz")

### Annotate direction of kendall trend into cds_event_time metadata
kendall.psi <- kendall.psi %>% 
    rownames_to_column("event.id") %>% 
    mutate(psi.Kend.dir = "Not-regulated") %>% 
    mutate(psi.Kend.dir = ifelse(psi.KeT < -co_KeT & psi.KeP < co_KeP, "Down", psi.Kend.dir)) %>% 
    mutate(psi.Kend.dir = ifelse(psi.KeT > co_KeT & psi.KeP < co_KeP, "Up", psi.Kend.dir)) 

cds.event.time <- cds.event.time %>% 
    left_join(kendall.psi %>% dplyr::select(event.id, psi.Kend.dir))


## Correlation 4: Pearson's correlation of PSI vs GEX (vst-transformed)
### Load factR object
factR.obj <- readRDS("rds/factR_final_transcriptome.rds")

# prep ASE to match selected as.nmd events
ase.factR <- ase(factR.obj) %>% 
    mutate(event.id = paste0(coord, strand, gene_id)) %>% 
    filter(event.id %in% rownames(whippet.psi.pivot)) %>% 
    filter(!is.na(ASNMDtype))
rownames(ase.factR) <- ase.factR$event.id
rownames(whippet.psi.pivot) <- ase.factR[rownames(whippet.psi.pivot),]$AS_id

rownames(ase.factR) <- ase.factR$AS_id
factR.obj@sets$AS@rowData <- ase.factR
factR.obj@sets$AS@data <- as.matrix(whippet.psi.pivot)


# prep GEX 
count <- DESeq2::counts(dds)
normcount <- DESeq2::counts(dds, normalized=TRUE)
colnames(count) <- colnames(whippet.psi.pivot)
colnames(normcount) <- colnames(whippet.psi.pivot)
factR.obj@sets$gene@counts <- count
factR.obj@sets$gene@data <- normcount
rownames(s2c_dmso) <- s2c_dmso$name
factR.obj@colData <- s2c_dmso
factR.obj@sets$AS@rowData <- ase.factR[ase.factR$gene_id %in% rownames(count),]

# run correlation
factR.obj <- testGeneCorr(factR.obj, alternative = "greater")

# output factR object
saveRDS(factR.obj, "rds/factR_final_transcriptome_gex_psi.rds")

# output correlation scores
ase(factR.obj) %>% 
    dplyr::select(AS_id, event.id, psi.gex.PeT = gene.cor.estimate, 
                  psi.gex.PeP = gene.cor.pval) %>% 
    write_tsv("tables/exons_gex_correlation.tsv.gz")


### Annotate Pearson's correlation significance into cds_event_time metadata
ase <- ase(factR.obj) %>% 
    mutate(psi.gex.sig = ifelse(gene.cor.estimate > 0 & gene.cor.pval < 0.05, TRUE, FALSE)) 
cds.event.time <- cds.event.time %>% 
    left_join(ase %>% dplyr::select(event.id, psi.gex.sig))


## Correlation 5: Pearson's correlation of PSI vs PTBP1 expression
### Load factR object


ptbp1 <- vst["ENSMUSG00000006498.17",]
pearson.psi.ptb <- do.call(bind_rows,apply(whippet.vst.pivot, 1, function(x){
    test <- cor.test( ptbp1[!is.na(x)],x[!is.na(x)], method = "pearson")
    data.frame(psi.ptb.PeR = test$estimate, psi.ptb.PeP = test$p.value)
}))
rownames(pearson.psi.ptb) <- rownames(whippet.vst.pivot)

pearson.psi.ptb %>% rownames_to_column("event.id") %>% 
    write_tsv("tables/event_PTBP_corr.tsv.gz")


### Output cds_event_time metadata
write_tsv(cds.event.time, "metadata/cds_event_time_meta.tsv")






























