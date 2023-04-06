## Load dependencies for this workflow
library(tidyverse)
library(factR2)
library(rtracklayer)
devtools::load_all("/home/cdn-bc/Github_repo/factR2/")

## set variables

### Working directory containing input data
wd <- "/media/cdn-bc/RAID/Projects/FHyyy_Anya_ES/gitrepo/data" 

## Set working directory
setwd(wd)

## Import custom transcriptome containing de novo transcripts
custom.gtf <- import("annotations/gffcmp.annotated.gtf")

## Retrieve reference transcriptome
ref.url <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/
release_M26/gencode.vM26.primary_assembly.annotation.gtf.gz"
ref.gtf <- import(str_remove(ref.url, "\n"))

## Unify custom and reference transcriptomes
union.seqlvls <- unique(c(seqlevels(custom.gtf), seqlevels(ref.gtf)))
seqlevels(custom.gtf) <- union.seqlvls
seqlevels(ref.gtf) <- union.seqlvls
union.gtf <- c(custom.gtf, ref.gtf)

## Create factR2 object
factR.obj <- createfactRObject(gtf = union.gtf,
                               reference = "vM26")

## Run factR pipeline
factR.obj <- runfactR(factR.obj)

## Compute conservation scores of intronic regions flanking exons 
factR.obj <- getAScons(factR.obj, type="flanks", padding=100)
factR.obj <- getAScons(factR.obj, type="upstream", padding=200)
factR.obj <- getAScons(factR.obj, type="downstream", padding=200)
factR.obj <- getAScons(factR.obj, type="upstream", padding=-100)
factR.obj <- getAScons(factR.obj, type="downstream", padding=-100)

## Retain AS-NMD events of CE,AA,AD,RI types
factR.obj@sets$AS@rowData <- factR.obj@sets$AS@rowData %>% 
    filter(!is.na(ASNMDtype), AStype %in% c("CE","AD","AA","RI"))


## Refine intron conservation scores based on exon type
### CE: 100bp upstream and 100bp downstream
### AD: 200bp downstream
### AA: 200bp upstream
### RI: First 100bp and last 100bp of intron
factR.obj@sets$AS@rowData <- factR.obj@sets$AS@rowData %>% 
    mutate(PhastCons.score = Cons.flanks.pad100) %>% 
    mutate(PhastCons.score = ifelse(as.character(AStype=="AD"), 
                                       Cons.downstream.pad200, 
                                    PhastCons.score)) %>% 
    mutate(PhastCons.score = ifelse(as.character(AStype=="AA"), 
                                       Cons.upstream.pad200, 
                                    PhastCons.score)) %>% 
    mutate(PhastCons.score = ifelse(as.character(AStype=="RI"), 
                                       (`Cons.downstream.pad-100`+
                                            `Cons.downstream.pad-100`)/2, 
                                    PhastCons.score)) %>% 
    dplyr::select(AS_id:ASNMD.in.cds, PhastCons.score)

## Output GTF
export(factR.obj@transcriptome, "annotations/final_transcriptome.gtf.gz")

## Output metadata at gene, transcript and event levels
genes(factR.obj) %>% 
    write_tsv("metadata/genes_metadata.tsv.gz")
txs(factR.obj) %>% 
    write_tsv("metadata/txs_metadata.tsv.gz")
ase(factR.obj) %>% 
    write_tsv("metadata/events_ASNMD_metadata.tsv.gz")

## Save factR object
saveRDS(factR.obj, "rds/factR_final_transcriptome.rds")



######## POTENTIALLY REMOVE

# output transcripts info and coordinates
factr.obj@transcriptome %>% 
    as.data.frame() %>% 
    filter(type == "transcript") %>% 
    write.table("tables/transcripts_info_final_transcriptome.tsv.gz",
                sep = "\t", quote = F, col.names = T, row.names = F)


















