library(tidyverse)
library(vroom)
library(ggpubr)

# set variables
wd <- "/media/cdn-bc/RAID/Projects/FHyyy_Anya_ES/gitrepo/data"
dco <-  0.1      # DeltaPsi cutoff
yco <-  0.9       # Probability cutoff
pco <- 1e-3         # P-value cutoff; applied to FDR
lco <-  log2(1)     # Log2 FC cutoff

## Set working directory
setwd(wd)

# Import NMD.exons data from factR (some of the events are actually introns but never mind)
AS.exons <- read.delim("tables/exons_final_transcriptome_factR2.NMD.tsv.gz", 
                       sep="\t", header=T)
NMD.exons <- AS.exons %>% 
    filter(!is.na(ASNMDtype))

# Import CHX-regulated exons data from Whippet.delta
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



# Select for CE, AD, AA, and RI events
NMD.exons <- filter(NMD.exons, AStype %in% c("CE","AD","AA","RI"))
CHX.exons.df <- filter(CHX.exons.df, Type %in% c("CE","AD","AA","RI"))

# Add gene names
tx2gene_meta <- read.delim("tables/transcripts_info_final_transcriptome.tsv.gz", 
                           sep="\t", header=T)
gn_name_id <- tx2gene_meta %>% 
    dplyr::select(gene_id, gene_name) %>% 
    distinct(gene_id, .keep_all = TRUE)

CHX.exons.df <- CHX.exons.df %>% 
    left_join(gn_name_id, by = c("Gene" = "gene_id"))

# Generate coordinate-strand-gene id columns:
NMD.exons <- NMD.exons %>% 
    mutate(id = paste0(coord, strand, gene_id))
CHX.exons.df <- CHX.exons.df %>% 
    mutate(id = paste0(Coord, Strand, Gene))


# Subset all Whippet-detected events that may cause NMD according to factR:
CHX.NMD.df  <-  CHX.exons.df %>% 
    filter(id %in% NMD.exons$id) %>% 
    left_join(NMD.exons %>% dplyr::select(id, ASNMDtype, ASNMD.in.cds))

# Filter significantly regulated by CHX Whippet events that may cause NMD
# Note that I am using conditional DeltaPsi filter to accommodate both poison (NMD-stimulating) and ORF-maintaining (NMD-repressing) events:
# I am also including cases where Probability <= yco but (Psi_A+Psi_B)/2 > 1-dco or (Psi_A+Psi_B)/2 < dco, accordingly,
# i.e. where there is no "room" for CHX to induce a splicing cha.nge of desired magnitude.
CHX.NMD.df <-  CHX.NMD.df  %>% 
    mutate(pass_positive = (Probability > yco & DeltaPsi > dco) | (Probability <= yco & (Psi_A+Psi_B)/2 > 1- dco)) %>% 
    mutate(pass_negative = (Probability > yco & DeltaPsi < -dco) | (Probability <= yco & (Psi_A+Psi_B)/2 < dco)) %>% 
    mutate(CHX.responsive = ifelse((ASNMDtype == "Stimulating" & pass_positive) | (ASNMDtype == "Repressing" & pass_negative),
                                   TRUE, FALSE)) %>% 
    dplyr::select(-pass_positive, -pass_negative)


CHX.NMD.sig.df <- filter(CHX.NMD.df, CHX.responsive)
CHX.NMD.sig.df %>% 
    group_by(File) %>% 
    tally()

# output dataframe
write_tsv(CHX.NMD.df, "tables/AS_primary_data.tsv")





































