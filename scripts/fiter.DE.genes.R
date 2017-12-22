library(dplyr)
library(tidyr)
library(tibble)
library(broom)

setwd("~/Box\ Sync/RNA-Seq/Sleep.Dep.dev/")

#####################################################
################### GENE LEVEL ######################
#####################################################
dedata <- read.table("RESULTS/DE.w.limmavoom/gene_level/DE_SS.vs.SD_sense.txt", header = TRUE)
qvalue.cutoff = 0.1
log2FC.cutoff = 0.5 #equals FC=1.5

#change the order of the samples by condition and timepoint
order.dedata = dedata[,c("id",
                         "B2.0hr","B11.0hr","C10.0hr","D9.0hr","E8.0hr","F7.0hr",
                         "B3.SS3hr","C2.SS3hr","C11.SS3hr","D10.SS3hr","E9.SS3hr","F8.SS3hr",
                         "B4.SS6hr","C3.SS6hr","D2.SS6hr","D11.SS6hr","E10.SS6hr","F9.SS6hr",
                         "B5.SS9hr","C4.SS9hr","D3.SS9hr","E2.SS9hr","F10.SS9hr",
                         "B6.SS12hr","C5.SS12hr","D4.SS12hr","E3.SS12hr","F2.SS12hr","F11.SS12hr",
                         "B7.SD3hr","C6.SD3hr","D5.SD3hr","E4.SD3hr","F3.SD3hr","G2.SD3hr",
                         "B8.SD6hr","C7.SD6hr","D6.SD6hr","E5.SD6hr","F4.SD6hr","G3.SD6hr",
                         "B9.SD9hr","C8.SD9hr","D7.SD9hr","E6.SD9hr","F5.SD9hr","G4.SD9hr",
                         "B10.SD12hr","C9.SD12hr","E7.SD12hr","F6.SD12hr","G5.SD12hr",
                         "geneCoordinate","geneSymbol","logFC","LimmaVoom.pvalue","BH.qvalue")]

# order.dedata %>%
  # filter(geneSymbol == "Arc")

#filter the data by q-value cutoff<=0.1
filter.dedata = order.dedata %>%
  filter(BH.qvalue <= qvalue.cutoff & abs(logFC) >= log2FC.cutoff)

#export tab delimited file with filtered genes
write.table(filter.dedata, file = "RESULTS/DE.w.limmavoom/gene_level/DE_SS.vs.SD_sense_filtered.txt", sep = "\t", row.names = FALSE)



#####################################################
################### GENE ANTISENSE ##################
#####################################################
dedata <- read.table("RESULTS/DE.w.limmavoom/gene_level/DE_SS.vs.SD_antisense.txt", header = TRUE)
qvalue.cutoff = 0.1
log2FC.cutoff = 0.5 #equals FC=1.5

#filter the data by q-value cutoff<=0.1
filter.dedata = dedata %>%
  filter(BH.qvalue <= qvalue.cutoff & abs(logFC) >= log2FC.cutoff)

#export tab delimited file with filtered genes
write.table(filter.dedata, file = "RESULTS/DE.w.limmavoom/gene_level/DE_SS.vs.SD_antisense_filtered.txt", sep = "\t", row.names = FALSE)



#####################################################
################### INTRON LEVEL ######################
#####################################################

dedata <- read.table("RESULTS/DE.w.limmavoom/intron_level/DE_SS.vs.SD_intron_sense.txt", header = TRUE, sep = "\t")
qvalue.cutoff = 0.1
log2FC.cutoff = 0.5 #equals FC=1.5

select.dedata = dedata[,!colnames(dedata) %in% c("NovelIntron", "highExp")]

#filter the data by q-value cutoff<=0.1
filter.dedata = select.dedata %>%
  filter(BH.qvalue <= qvalue.cutoff & abs(logFC) >= log2FC.cutoff)

#export tab delimited file with filtered genes
write.table(filter.dedata, file = "RESULTS/DE.w.limmavoom/intron_level/DE_SS.vs.SD_intron.sense_filtered.txt", sep = "\t", row.names = FALSE)




#####################################################
################### INTRON ANTISENSE ################
#####################################################

dedata <- read.table("RESULTS/DE.w.limmavoom/intron_level/DE_SS.vs.SD_intron_antisense.txt", header = TRUE, sep = "\t")
qvalue.cutoff = 0.1
log2FC.cutoff = 0.5 #equals FC=1.5

select.dedata = dedata[,!colnames(dedata) %in% c("NovelIntron", "highExp")]

#filter the data by q-value cutoff<=0.1
filter.dedata = select.dedata %>%
  filter(BH.qvalue <= qvalue.cutoff & abs(logFC) >= log2FC.cutoff)

#export tab delimited file with filtered genes
write.table(filter.dedata, file = "RESULTS/DE.w.limmavoom/intron_level/DE_SS.vs.SD_intron.antisense_filtered.txt", sep = "\t", row.names = FALSE)
