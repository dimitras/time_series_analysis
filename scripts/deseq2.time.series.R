library(DESeq2)
library(dplyr)
library(ggplot2)

setwd("~/Box Sync/RNA-Seq/Sleep.Dep.dev/")
counts = read.delim("master_list_of_gene_counts_MIN.sense.UNNORMQUANT.REORDERED.txt", row.names = "id")

counts.select = counts %>%
  select(B2.CTL0hr:G5.SD12hr) %>%
  # select(-B2.CTL0hr, -B11.CTL0hr, -C10.CTL0hr, -E8.CTL0hr, -F7.CTL0hr, -D9.CTL0hr, -geneCoordinate, -geneSymbol) %>%
  # select(-B2.CTL0hr, -B11.CTL0hr, -C10.CTL0hr, -E8.CTL0hr, -F7.CTL0hr, -D9.CTL0hr, -D10.SS3hr, -C3.SS6hr, -B6.SS12hr, -E4.SD3hr, -F4.SD6hr, -G4.SD9hr, -geneCoordinate, -geneSymbol) %>%
  # select(-D9.CTL0hr, -D10.SS3hr, -C3.SS6hr, -B6.SS12hr, -E4.SD3hr, -F4.SD6hr, -G4.SD9hr, -geneCoordinate, -geneSymbol) %>%
  # filter(!D9.CTL0hr, !D10.SS3hr, !C3.SS6hr, !B6.SS12hr, !E4.SD3hr, !F4.SD6hr, !G4.SD9hr, !geneCoordinate, !geneSymbol) %>%
  as.matrix()

coldata = data.frame(row.names = colnames(counts.select))
coldata$timepoint = as.factor(gsub("^[B-G][[:digit:]]+.(SS|SD|CTL)([[:digit:]]+)hr","\\2", row.names(coldata)))
coldata$condition = as.factor(gsub("^[B-G][[:digit:]]+.(SS|SD|CTL)([[:digit:]]+)hr","\\1", row.names(coldata)))
coldata$timepoint = factor(coldata$timepoint, levels = c(0, 3, 6, 9, 12))
# coldata$rep = as.factor(gsub("^([B-G][[:digit:]]+).(SS|SD|CTL)([[:digit:]]+)hr","\\1", row.names(coldata)))

# coldata = transform(coldata, nest = as.numeric(interaction(timepoint, condition, drop=TRUE)))

# coldata$group = as.factor(gsub("^[B-G][[:digit:]]+.((SS|SD|CTL)[[:digit:]]+)hr","\\1", row.names(coldata)))

#flag the nested pairs, remove the zeros, rerun DESeq w the new model and betaPrior=FALSE)
# mm = model.matrix( ~ condition + timepoint + condition:nest + condition:timepoint, coldata)
# res = DESeq(counts.deseq.dataset, full=mm, betaPrior=FALSE)

#same error
#https://support.bioconductor.org/p/95602/#95633
mm = model.matrix( ~ timepoint + condition:timepoint, coldata)
mm.full = as.matrix(mm[,-c(6, 11)])
mm.reduced = as.matrix(mm.full[,-c(2:5)]) 

dds = DESeqDataSetFromMatrix(countData = counts.select,
                             colData = coldata,
                             design = ~ condition)

dds = DESeq(dds, full=mm.full, reduced=mm.reduced, test="LRT")
res <- results(dds)

#group the condition and time
counts.deseq.dataset = DESeqDataSetFromMatrix(countData = counts.select,
                                              colData = coldata,
                                              design = ~ group)
counts.deseq.dataset = estimateSizeFactors(counts.deseq.dataset)
counts.deseq.dataset = estimateDispersions(counts.deseq.dataset)
counts.deseq.dataset = nbinomWaldTest(counts.deseq.dataset)
res = DESeq(counts.deseq.dataset)
out = results(res)
resultsNames(res)
summary(out)


##### time series according to manual
counts.deseq.dataset = DESeqDataSetFromMatrix(countData = counts.select,
                                              colData = coldata,
                                              design = ~ condition + timepoint + condition:timepoint)
                                              #design = ~ condition + timepoint + rep + condition:timepoint + timepoint:rep + condition:rep)

colData(counts.deseq.dataset)
# dsq = DESeq(counts.deseq.dataset)
dsq = DESeq(counts.deseq.dataset, test="LRT", reduced = ~ condition + timepoint)
# dsq = DESeq(counts.deseq.dataset, test="LRT", full = ~ condition + timepoint + condition:timepoint, reduced = ~ condition + timepoint)
res = results(dsq)
# res = results(dsq, contrast=c("condition", "SS", "SD"))
summary(res)
resultsNames(dsq)
plotDispEsts(dsq)
plotMA(res)
plotCounts(dsq, gene=which.min(res$padj), intgroup="condition")
plotCounts(dsq, gene=which.min(res$padj), intgroup="timepoint")

##### Plot with ggplot
#plot the counts for the groups over time using ggplot2, 
#for the gene with the smallest adjusted p value, testing for condition-dependent time profile 
#and accounting for differences at time 0 (figure below). 
#Keep in mind that the interaction terms are the difference between the two groups at a given time 
#after accounting for the difference at time 0.
fiss <- plotCounts(dsq, which.min(res$padj), 
                   intgroup = c("condition","timepoint"), returnData = TRUE)
ggplot(fiss, aes(x = timepoint, y = count, color = condition, group = condition)) + 
  geom_point() + geom_smooth(se = FALSE, method = "loess") + scale_y_log10() + geom_text(aes(label=rownames(fiss)),hjust=0, vjust=0, size = 2.2)


#####write results to a spreadsheet
res = results(dsq, contrast=c("condition", "SS", "SD"))
df = cbind(key = rownames(res), as.data.frame(res))
write.table(df, file="SS.vs.SD.w.deseq2.results.txt", sep = "\t", quote = FALSE, row.names = FALSE)


#more plotting?
counts.deseq.dataset = estimateSizeFactors(counts.deseq.dataset)
counts.deseq.dataset = estimateDispersions(counts.deseq.dataset)
counts.deseq.dataset = nbinomWaldTest(counts.deseq.dataset)

design(counts.deseq.dataset) = formula(~ condition)# + timepoint + condition:timepoint)

res = DESeq(counts.deseq.dataset, test="LRT", reduced = ~ condition + timepoint + condition:timepoint)
out = results(res)
resultsNames(res)
summary(out)
plotDispEsts(res)
plotMA(out)
plotCounts(res, gene=which.min(out$padj), intgroup="condition")
ntd = normTransform(res)
library("vsn")
meanSdPlot(assay(ntd))

counts.deseq.dataset$group <- factor(paste0(counts.deseq.dataset$timepoint, counts.deseq.dataset$condition))
design(counts.deseq.dataset) = ~ group
res = DESeq(res)
out = results(res)
plotMA(out)


# check out solutions:
# https://support.bioconductor.org/p/64480/
# https://support.bioconductor.org/p/93936/
# http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#model-matrix-not-full-rank
# http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#multi-factor-designs
# http://www.bioconductor.org/help/workflows/rnaseqGene/#time

# https://support.bioconductor.org/p/95929/
# https://support.bioconductor.org/p/95602/#95633

  