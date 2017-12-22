library(DESeq2)
library(dplyr)

setwd("~/Box Sync/RNA-Seq/Sleep.Dep.dev/")
counts = read.delim("master_list_of_gene_counts_MIN.sense.UNNORMQUANT.REORDERED.txt")
row.names(counts) = counts[,1]

counts.select = counts %>%
  select(B2.CTL0hr:G5.SD12hr) %>%
  as.matrix()

coldata = as.data.frame(colnames(counts.select))
coldata$timepoint = as.factor(gsub("^[B-G][[:digit:]]+.(SS|SD|CTL)([[:digit:]]+)hr","\\2", coldata[,1]))
coldata$condition = as.factor(gsub("^[B-G][[:digit:]]+.(SS|SD|CTL)([[:digit:]]+)hr","\\1", coldata[,1]))

row.names(coldata) = coldata[,1]
coldata.select = coldata %>%
  select(-`colnames(counts.select)`)

counts.deseq.dataset = DESeqDataSetFromMatrix(countData = counts.select,
                                              colData = coldata.select,
                                              design = ~ condition + timepoint)# + condition:timepoint)


