### load libraries### 
library(Seurat)
library(tidyverse)
library(rtracklayer)
library(Biobase)
library(reshape2)
library(ggpubr)
library(MuSiC)

### get gene names from latest gencode annotation ### 
gtf <- import("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M29/gencode.vM29.annotation.gtf.gz")
genenames <- gtf %>% 
    as.data.frame() %>% 
    distinct(gene_id, gene_name) %>% 
    mutate(gene_id_short = str_remove(gene_id, ".[0-9]+$"))

### prepare mESC reference### 
# retrieve mESC single-cell dataset
download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE79578&format=file",
              "GSE79578_RAW.tar")
untar("GSE79578_RAW.tar")

# create count matrix
esc.counts <- read_tsv("GSM2098554_smartseq_2i.txt.gz", col_names = F)
esc.counts <- rowsum(esc.counts[,-1], esc.counts$X1, reorder = T)  # combine counts from same features
esc.counts <- esc.counts[,!is.na(colSums(esc.counts))] # remove cells with NA feature counts
esc.counts <- esc.counts[(rowSums(esc.counts)>0) ,] # remove not-detected features
esc.counts <- esc.counts[rownames(esc.counts) %in% genenames$gene_name,] # remove inconsistent named features

# create seurat object
esc.seurat <- CreateSeuratObject(esc.counts, project = "ESC")
esc.seurat$subclass_label <- "ESC"
esc.seurat$class_label <- "ESC"
esc.seurat$region_label <- "ESC"


### prepare NPC reference ### 
# retrieve gene expression matrix
download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE67833&format=file&file=GSE67833%5FGene%5Fexpression%5Fmatrix%5FGSM1684656%2D704%2Ecsv%2Egz",
              "GSE67833_Gene_expression_matrix_GSM1684656-704.csv.gz")
npc.counts <- read_csv("GSE67833_Gene_expression_matrix_GSM1684656-704.csv.gz", 
                       col_names = T) %>% 
    column_to_rownames("...1")

# create counts
npc.genes <- genenames %>% 
    filter(gene_id_short %in% rownames(npc.counts))
npc.counts <- rowsum(npc.counts[npc.genes$gene_id_short,], npc.genes$gene_name)
npc.counts <- npc.counts[,!is.na(colSums(npc.counts))] # remove cells with NA feature counts
npc.counts <- npc.counts[(rowSums(npc.counts)>0) ,] # remove not-detected features

# create metadata
npc.seurat <- CreateSeuratObject(npc.counts, project = "NPC")
npc.seurat$subclass_label <- "NPC"
npc.seurat$class_label <- "NPC"
npc.seurat$region_label <- "NPC"


### prepare adult cortex reference ### 
# retrieve Seurat object
options(timeout=10000)
download.file("https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_mouse_ctx-hpf_smart-seq/Seurat.ss.rda",
              "Seurat.ss.rda")
load("Seurat.ss.rda")

# sample cell clusters from 5 brain regions
sampled.ss.df <- ss.seurat@meta.data %>% 
    rownames_to_column("cell_id") %>% 
    filter(region_label %in% c("VISp", "ALM", "HIP", "SSp", "MOp")) %>% 
    mutate(subclass_region = paste0(region_label,"_", subclass_label)) %>% 
    group_by(subclass_region) %>% 
    filter(n() > 20) %>% 
    sample_n(ifelse(n() > 150, 150, n()))
ss.seurat.filtered <- ss.seurat[,colnames(ss.seurat) %in% sampled.ss.df$cell_id] 


### create pooled reference seurat object ###
ref.seurat <- merge(ss.seurat.filtered, y = c(esc.seurat, npc.seurat))

# clean up data
rm(esc.counts, npc.counts, ss.seurat.filtered)

### prepare induced neuron dataset (bulk RNAseq) ### 
# load and process induced neuron dataset
dds <- readRDS("/media/cdn-bc/RAID/Projects/FHyyy_Anya_ES/gitrepo/data/rds/DDS_DMSO_only.rds")
gene.meta <- read_tsv("/media/cdn-bc/RAID/Projects/FHyyy_Anya_ES/gitrepo/data/metadata/genes_metadata.tsv.gz")
gene.meta <- gene.meta %>%
    filter(gene_id %in% rownames(dds))
iNs.bulk.counts <- counts(dds)
iNs.bulk.counts <- rowsum(iNs.bulk.counts[gene.meta$gene_id,], gene.meta$gene_name)
iNs.bulk.counts <- iNs.bulk.counts[!is.na(rownames(iNs.bulk.counts)),]
colnames(iNs.bulk.counts) <- paste(dds$day, dds$batch, sep = "_")


### Prepare expressiondataset for deconvolution ### 
# prepare test dataset
bulk.metadata <- data.frame(row.names = colnames(iNs.bulk.counts),
                       group = str_sub(colnames(iNs.bulk.counts), end = -3))
bulk.metalabels <- data.frame(labelDescription=c("group"))
bulk.metadata <- AnnotatedDataFrame(data=bulk.metadata, varMetadata=bulk.metalabels)
bulk.eset <- ExpressionSet(as.matrix(iNs.bulk.counts), phenoData = bulk.metadata)

# prepare reference dataset
ref.metadata <- ref.seurat@meta.data %>% as.data.frame() %>% 
    mutate(group = paste0(region_label,"_", subclass_label)) %>% 
    dplyr::select(group)
ref.metadata$sample_id <- rownames(ref.metadata)
ref.metalabels <- data.frame(labelDescription=colnames(ref.metadata))
ref.metadata <- AnnotatedDataFrame(data=ref.metadata, varMetadata=ref.metalabels)
ref.eset <- ExpressionSet(as.matrix(ref.seurat@assays$RNA@counts), phenoData = ref.metadata)

### Run deconvolution using MuSiC ### 
music.out <- music_prop(bulk.eset = bulk.eset, sc.eset = ref.eset, clusters = 'group',
                                      samples = 'sample_id', verbose = T)
### Clean up output and visualise ### 
# clean up MuSiC output
m.music.out <-  rbind(melt(music.out$Est.prop.weighted), 
                      melt(music.out$Est.prop.allgene))
colnames(m.music.out) = c('Sub', 'CellType', 'Prop')
m.music.out$Method = factor(rep(c('MuSiC', 'NNLS'), each = nrow(m.music.out)/2), levels = c('MuSiC', 'NNLS'))

# append class information and sum proportions for each class
class_categories <- ref.seurat@meta.data %>% 
    distinct(subclass_label, class_label, region_label) %>% 
    mutate(region_subclass = paste0(region_label,"_", subclass_label)) %>% 
    dplyr::select("class_label",  CellType = "region_subclass")
m.music.out.sum <- m.music.out %>% left_join(class_categories) %>% 
    mutate(group = str_remove(as.character(Sub), "_[A,B]$")) %>% 
    mutate(group = str_replace(group, "CNT_", "Day")) %>% 
    mutate(group = factor(group, levels = c("Day0", "Day3", "Day6", "Day12", "Day24"))) %>% 
    group_by(Sub, group, class_label, Method) %>% 
    summarise(Prop = sum(Prop)) %>%
    filter(Method == "MuSiC")

write_tsv(m.music.out.sum, "/media/cdn-bc/RAID/Projects/FHyyy_Anya_ES/gitrepo/data/tables/MuSiC_DMSO_proportions.tsv")






















