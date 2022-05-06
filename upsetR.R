rm(list=ls()) #reset working space
graphics.off() #closing current or all graphical windows

library(ggplot2)
library(UpSetR)
library(grid)

setwd("~/MA_Bioinformatics/Thesis_Project/04_overlap_analysis")

common_originals <- read.table(file = 'bed_outputs/deletions/common.with.originals.bed', 
                               sep = '\t', header = FALSE)
unique_common <- read.table(file = 'bed_outputs/deletions/sv_80_overlap.bed', sep = '\t', header = FALSE)

#filter for scaffolds
common_originals <- common_originals[!grepl("caffold", common_originals$V1),]
unique_common <- unique_common[!grepl("caffold", unique_common$V1),]

length(unique(common_originals$V2))
length(unique(unique_common$V2))

cnvnator_bed <- read.table(file = 'bed_outputs/deletions/cnvnator_del.bed', sep = '\t', header = FALSE)
delly_bed <- read.table(file = 'bed_outputs/deletions/delly_del.bed', sep = '\t', header = FALSE)
dysgu_bed <- read.table(file = 'bed_outputs/deletions/dysgu_del.bed', sep = '\t', header = FALSE)
manta_bed <- read.table(file = 'bed_outputs/deletions/manta_del.bed', sep = '\t', header = FALSE)
matchclips_bed <- read.table(file = 'bed_outputs/deletions/matchclips_del.bed', sep = '\t', header = FALSE)
tiddit_bed <- read.table(file = 'bed_outputs/deletions/tiddit_del.bed', sep = '\t', header = FALSE)

length(c(delly_IDs, setdiff(delly_bed$V4, common_originals$V16)))
#cnvnator_id<- length(df[grepl("delly", df$V16),][ , "V16"])


plot_overlap <- function(sv_size_range){
  #sv_size_range <- 'S'
  sv_sizes <- list(SS=c(50,100), S=c(100,1000), M=c(1000,100000), L=c(100000,1000000))
  sv_size = get(sv_size_range, sv_sizes)
  #filter for sv size range
  common_originals <- common_originals[(common_originals$V3 - common_originals$V2) > sv_size[1] &
                                         (common_originals$V3 - common_originals$V2) < sv_size[2], ]
  unique_common <- unique_common[(unique_common$V3 - unique_common$V2) > sv_size[1] &
                                   (unique_common$V3 - unique_common$V2) < sv_size[2], ]
  
  cnvnator_bed <- cnvnator_bed[(cnvnator_bed$V3 - cnvnator_bed$V2) > sv_size[1] &
                           (cnvnator_bed$V3 - cnvnator_bed$V2) < sv_size[2], ]
  delly_bed <- delly_bed[(delly_bed$V3 - delly_bed$V2) > sv_size[1] &
                           (delly_bed$V3 - delly_bed$V2) < sv_size[2], ]
  dysgu_bed <- dysgu_bed[(dysgu_bed$V3 - dysgu_bed$V2) > sv_size[1] &
                           (dysgu_bed$V3 - dysgu_bed$V2) < sv_size[2], ]
  manta_bed <- manta_bed[(manta_bed$V3 - manta_bed$V2) > sv_size[1] &
                           (manta_bed$V3 - manta_bed$V2) < sv_size[2], ]
  matchclips_bed <- matchclips_bed[(matchclips_bed$V3 - matchclips_bed$V2) > sv_size[1] &
                           (matchclips_bed$V3 - matchclips_bed$V2) < sv_size[2], ]
  tiddit_bed <- tiddit_bed[(tiddit_bed$V3 - tiddit_bed$V2) > sv_size[1] &
                           (tiddit_bed$V3 - tiddit_bed$V2) < sv_size[2], ]
  
  unique_common$uniq_id <- 1:nrow(unique_common) 
  #print(nrow(unique_common))

  #making data set for UpSetR
  data <- unique_common[, c("V2", "V3", "uniq_id")]
  df <- merge(common_originals, data, by=c("V2","V3"))
  
  cnvnator_IDs<- df[grepl("cnvnator", df$V16),][ , "uniq_id"] #these are full of non unique SVs, but upsetr only looks at unique svs occuring between software
  delly_IDs <- df[grepl("delly", df$V16),][ , "uniq_id"]
  dysgu_IDs <- df[grepl("dysgu", df$V16),][ , "uniq_id"]
  manta_IDs <- df[grepl("manta", df$V16),][ , "uniq_id"]
  matchclips_IDs <- df[grepl("matchclips", df$V16),][ , "uniq_id"]
  tiddit_IDs <- df[grepl("tiddit", df$V16),][ , "uniq_id"]
  
  cnvnator_IDs <- c(cnvnator_IDs, setdiff(cnvnator_bed$V4, df$V16))
  delly_IDs <- c(delly_IDs, setdiff(delly_bed$V4, df$V16))
  dysgu_IDs <- c(dysgu_IDs, setdiff(dysgu_bed$V4, df$V16))
  manta_IDs <- c(manta_IDs, setdiff(manta_bed$V4, df$V16))
  matchclips_IDs <- c(matchclips_IDs, setdiff(matchclips_bed$V4, df$V16))
  tiddit_IDs <- c(tiddit_IDs, setdiff(tiddit_bed$V4, df$V16))
  
  listInput <- list(cnvnator = cnvnator_IDs, delly = delly_IDs, dysgu = dysgu_IDs,
                    manta = manta_IDs, matchclips = matchclips_IDs, tiddit = tiddit_IDs)
  Upset <- upset(fromList(listInput), order.by = "freq" , nsets = 6)
  print(Upset)
  grid.text(sv_size_range,x = 0.65, y=0.95, gp=gpar(fontsize=12))
  #  print(Upset)
}

plot_overlap('M')



