library(tidyverse)
library(DESeq2)
devtools::load_all("/home/cdn-bc/Github_repo/factR2/")

# set variables
wd <- "/media/cdn-bc/RAID/Projects/FHyyy_Anya_ES/gitrepo/data"

## Set working directory
setwd(wd)

# Import samples metadata
s2c_master <- read.table("metadata/samples.txt", header = TRUE)
s2c_dmso <- s2c_master %>% filter(condition == "DMSO")

## Import Whippet PSI data
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
    mutate(id = paste0(Coord, Strand, Gene))

## import AS-NMD data
asnmd <- read_tsv("tables/deg_AS_CDSs_CHX_vs_DMSO_bytimepoint.tsv")

## Simplify asnmd dataframe and append PSI values from each sample
asnmd <- asnmd %>% 
    dplyr::select(gene_id, id, Type, ASNMDtype, ASNMD.in.cds) %>% 
    filter(!is.na(id))
whippet.psi.pivot <- whippet.psis %>% 
    dplyr::select(id, Name, Psi) %>% 
    filter(id %in% asnmd$id) %>% 
    pivot_wider(names_from = Name, values_from = Psi) %>% 
    column_to_rownames("id")

# Load factR object
factR.obj <- readRDS("rds/factR_final_transcriptome.rds")

# prep ASE to match selected as.nmd events
ase.factR <- ase(factR.obj) %>% 
    mutate(id = paste0(coord, strand, gene_id)) %>% 
    filter(id %in% rownames(whippet.psi.pivot)) %>% 
    filter(!is.na(ASNMDtype))
rownames(ase.factR) <- ase.factR$id
rownames(whippet.psi.pivot) <- ase.factR[rownames(whippet.psi.pivot),]$AS_id

rownames(ase.factR) <- ase.factR$AS_id
factR.obj@sets$AS@rowData <- ase.factR
factR.obj@sets$AS@data <- as.matrix(whippet.psi.pivot)

# prep GEX 
dds <- readRDS("rds/DDS_all_samples.rds")
dds <- dds[,s2c_dmso$sample]
dds <- estimateSizeFactors(dds)
count <- DESeq2::counts(dds)
normcount <- DESeq2::counts(dds, normalized=TRUE)
colnames(count) <- colnames(whippet.psi.pivot)
colnames(normcount) <- colnames(whippet.psi.pivot)
factR.obj@sets$gene@data <- normcount
factR.obj@sets$gene@counts <- count

rownames(s2c_dmso) <- s2c_dmso$name
factR.obj@colData <- s2c_dmso


# run correlation
factR.obj <- testGeneCorr(factR.obj, alternative = "greater")

# output factR object
saveRDS(factR.obj, "rds/factR_final_transcriptome_gex_psi.rds")

# output correlation scores
ase(factR.obj) %>% 
    dplyr::select(AS_id, id,gene.cor.estimate, gene.cor.pval) %>% 
    write_tsv("tables/exons_gex_correlation.tsv.gz")



# run Kendall correlation on psi
asinTransform = function(x) { 2 * asin(sqrt(x))/pi }
whippet.psi.pivot <- asinTransform(as.matrix(factR.obj2@sets$AS@data ))
whippet.psi.pivot <- whippet.psi.pivot[rowSums(is.na(whippet.psi.pivot)) <= 5, ]

timepoint <- rep(1:5, each = 2)
kendall.psi <- do.call(bind_rows,apply(whippet.psi.pivot, 1, function(x){
    test <- cor.test( timepoint[!is.na(x)],x[!is.na(x)], method = "kendall")
    data.frame(KeT = test$estimate, KeP = test$p.value)
}))

rownames(kendall.psi) <- rownames(whippet.psi.pivot)
kendall.psi$id <- ase.factR[rownames(kendall.psi),]$id
kendall.psi %>% 
    write_tsv("tables/exons_kendall.tsv.gz")


deseq.counts <- factR.obj2@sets$gene@counts
deseq.vst <- varianceStabilizingTransformation(deseq.counts)
kendall.vst <- do.call(bind_rows,apply(deseq.vst, 1, function(x){
    test <- cor.test(timepoint[!is.na(x)],x[!is.na(x)],  method = "kendall")
    data.frame(KeT = test$estimate, KeP = test$p.value)
}))
rownames(kendall.vst) <- rownames(deseq.vst)
kendall.vst %>% 
    mutate(gene_id = rownames(deseq.vst)) %>% 
    write_tsv("tables/gex_kendall.tsv.gz")


































