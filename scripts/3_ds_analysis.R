## Load dependencies for this workflow
library(tidyverse)
library(vroom)

# set variables
wd <- "/media/cdn-bc/RAID/Projects/FHyyy_Anya_ES/gitrepo/data"
dco <-  0.1      # DeltaPsi cutoff
yco <-  0.9       # Probability cutoff
pco <- 1e-3         # P-value cutoff; applied to FDR
lco <-  log2(1)     # Log2 FC cutoff

## Set working directory
setwd(wd)

## Import NMD.exons data from factR 
NMD.exons <- read_tsv("metadata/events_ASNMD_metadata.tsv.gz")


## Import CHX-regulated exons data from Whippet.delta
whippet.outs <- list.files("whippet/iNs", full.names = TRUE)
CHX.exons.df <- vroom(whippet.outs, 
                      id = "Timepoint")
CHX.exons.df$Timepoint <- str_remove(CHX.exons.df$Timepoint,
                                     "whippet/iNs/")
CHX.exons.df$Timepoint <- str_remove(CHX.exons.df$Timepoint,
                                     "_CHX_v_DMSO.diff.gz")
CHX.exons.df$Timepoint <- factor(CHX.exons.df$Timepoint, 
                                 c("Day0","Day3","Day6","Day12","Day24"))

## Add gene names
tx.meta <- read_tsv("metadata/txs_metadata.tsv.gz")
gn_name_id <- tx.meta %>% 
    dplyr::select(gene_id, gene_name) %>% 
    distinct(gene_id, .keep_all = TRUE)

CHX.exons.df <- CHX.exons.df %>% 
    left_join(gn_name_id, by = c("Gene" = "gene_id"))

## Generate coordinate-strand-gene id columns:
NMD.exons <- NMD.exons %>% 
    mutate(event.id = paste0(coord, strand, gene_id))
CHX.exons.df <- CHX.exons.df %>% 
    mutate(event.id = paste0(Coord, Strand, Gene))


## Subset all Whippet-detected events that may cause NMD according to factR:
CHX.NMD.df  <-  CHX.exons.df %>% 
    filter(event.id %in% NMD.exons$event.id) %>% 
    left_join(NMD.exons %>% dplyr::select(event.id, ASNMDtype, ASNMD.in.cds))


## Filter significantly regulated by CHX Whippet events that may cause NMD
### Note that I am using conditional DeltaPsi filter to accommodate both poison (NMD-stimulating) and ORF-maintaining (NMD-repressing) events:
### I am also including cases where Probability <= yco but (Psi_A+Psi_B)/2 > 1-dco or (Psi_A+Psi_B)/2 < dco, accordingly,
### i.e. where there is no "room" for CHX to induce a splicing cha.nge of desired magnitude.
CHX.NMD.df <-  CHX.NMD.df  %>% 
    mutate(pass_positive = (Probability > yco & DeltaPsi > dco) | 
               (Probability <= yco & (Psi_A+Psi_B)/2 > 1- dco)) %>% 
    mutate(pass_negative = (Probability > yco & DeltaPsi < -dco) | 
               (Probability <= yco & (Psi_A+Psi_B)/2 < dco)) %>% 
    mutate(ds.CHX.sig = ifelse((ASNMDtype == "Stimulating" & pass_positive) | 
                                       (ASNMDtype == "Repressing" & pass_negative),
                                  TRUE, FALSE)) %>% 
    dplyr::select(-pass_positive, -pass_negative)
        

## Save results as TSV file
write_tsv(CHX.NMD.df, "tables/ds_iNs_CHX_vs_DMSO_bytimepoint.tsv")


## Create cds_event_time metadata
dge.cds <- read_tsv("tables/dge_CDSs_CHX_vs_DMSO_bytimepoint.tsv")
CHX.NMD.df <- mutate(CHX.NMD.df, factR.event = TRUE)
cds.event.time.meta <- dge.cds %>% 
    dplyr::select(gene_id,Timepoint) %>% 
    left_join(CHX.NMD.df %>% dplyr::select(gene_id=Gene, Timepoint, 
                                           event.id, Type, factR.event,
                                           ASNMDtype,ds.CHX.sig)) %>% 
    replace_na(list(factR.event = FALSE))

## Update gene-level upregulation to CHX treatment
dge.cds <- dge.cds %>% 
    mutate(dge.CHX.up = ifelse(log2FoldChange > lco & padj < pco, TRUE, FALSE))

cds.event.time.meta <- cds.event.time.meta %>% 
    left_join(dge.cds %>% dplyr::select(gene_id, Timepoint, dge.CHX.up))
    

write_tsv(cds.event.time.meta, "metadata/cds_event_time_meta.tsv")

## Repeat for primary neurons dataset

## Import CHX-regulated exons data from Whippet.delta
whippet.outs <- list.files("whippet/primary", full.names = TRUE)
CHX.exons.df <- vroom(whippet.outs, 
                      id = "File")
CHX.exons.df$File <- str_remove(CHX.exons.df$File,
                                "whippet/primary/")
CHX.exons.df$File <- str_remove(CHX.exons.df$File,
                                ".diff.gz")
CHX.exons.df <- CHX.exons.df %>% 
    mutate(Sample = str_split_i(File, "_", -1)) %>% 
    mutate(Comparison = str_remove(File, paste0("_",Sample))) %>% 
    dplyr::select(File, Sample, Comparison, Gene:Complexity)

CHX.exons.df <- CHX.exons.df %>% 
    left_join(gn_name_id, by = c("Gene" = "gene_id"))

CHX.exons.df <- CHX.exons.df %>% 
    mutate(event.id = paste0(Coord, Strand, Gene))

CHX.NMD.df  <-  CHX.exons.df %>% 
    filter(event.id %in% NMD.exons$event.id) %>% 
    left_join(NMD.exons %>% dplyr::select(event.id, ASNMDtype, ASNMD.in.cds))

CHX.NMD.df <-  CHX.NMD.df  %>% 
    mutate(pass_positive = (Probability > yco & DeltaPsi > dco) | 
               (Probability <= yco & (Psi_A+Psi_B)/2 > 1- dco)) %>% 
    mutate(pass_negative = (Probability > yco & DeltaPsi < -dco) | 
               (Probability <= yco & (Psi_A+Psi_B)/2 < dco)) %>% 
    mutate(ds.CHX.sig = ifelse((ASNMDtype == "Stimulating" & pass_positive) | 
                                       (ASNMDtype == "Repressing" & pass_negative),
                                   TRUE, FALSE)) %>% 
    dplyr::select(-pass_positive, -pass_negative)


## Save results as TSV file
write_tsv(CHX.NMD.df, "tables/ds_primary_CHX_vs_DMSO_bytimepoint.tsv")













































