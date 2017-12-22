library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(broom)
library(genefilter)

setwd("~/Box\ Sync/RNA-Seq/Sleep.Dep.dev/")

sel <- read.table("RESULTS/DE.w.limmavoom/gene_level/DE_SS.vs.SD_sense_filtered.txt", header = TRUE, sep = "\t")

sel$mean.SS0hr = sel %>%
  select(B2.0hr:F7.0hr) %>%
  rowMeans()

sel$mean.SS3hr = sel %>%
  select(B3.SS3hr:F8.SS3hr) %>%
  rowMeans()

sel$mean.SS6hr = sel %>%
  select(B4.SS6hr:F9.SS6hr) %>%
  rowMeans() 

sel$mean.SS9hr = sel %>%
  select(B5.SS9hr:F10.SS9hr) %>%
  rowMeans() 

sel$mean.SS12hr = sel %>%
  select(B6.SS12hr:F11.SS12hr) %>%
  rowMeans() 

sel$mean.SD0hr = sel$mean.SS0hr

sel$mean.SD3hr = sel %>%
  select(B7.SD3hr:G2.SD3hr) %>%
  rowMeans()

sel$mean.SD6hr = sel %>%
  select(B8.SD6hr:G3.SD6hr) %>%
  rowMeans() 

sel$mean.SD9hr = sel %>%
  select(B9.SD9hr:G4.SD9hr) %>%
  rowMeans() 

sel$mean.SD12hr = sel %>%
  select(B10.SD12hr:G5.SD12hr) %>%
  rowMeans() 


sel.m = sel %>%
  select(id, geneSymbol, mean.SS0hr:mean.SD12hr) %>%
  gather(Sample_id,MeanExpression,mean.SS0hr:mean.SD12hr) %>%
  mutate(Condition=gsub("^mean.(SS|SD)([[:digit:]]+hr)","\\1",Sample_id),
         TP=gsub("^mean.(SS|SD)([[:digit:]]+)hr","\\2",Sample_id),
         TP=factor(TP, levels = c("0","3","6","9","12")),
         Rep=gsub("^mean.(SS|SD\\[[:digit:]]+hr)","\\1",Sample_id))
  
sel.m %>%  
  ggplot(aes(x=TP, y=MeanExpression, group=Condition, color=Condition)) + 
  geom_path(size=2) +
  facet_wrap(~geneSymbol, ncol=10, scales = "free") +
  theme_bw(base_size = 30) +
  theme(legend.position = c(0.94,0.001), legend.title = element_text(size = 30),
        legend.text = element_text(size = 30), legend.key.height = unit(2,"cm")) +
  labs(y="Mean expression", x="Timepoint (in hours)")
  
ggsave("RESULTS/Longtitudinal/longtitudinal.sense.genes.filtered.w.free.scale.pdf", height = 120, width = 120, units ="cm")


########### ANTISENSE ########### 

sel <- read.table("RESULTS/DE.w.limmavoom/gene_level/DE_SS.vs.SD_antisense_filtered.txt", header = TRUE, sep = "\t")

sel$mean.SS0hr = sel %>%
  select(B2.0hr:F7.0hr) %>%
  rowMeans()

sel$mean.SS3hr = sel %>%
  select(B3.SS3hr:F8.SS3hr) %>%
  rowMeans()

sel$mean.SS6hr = sel %>%
  select(B4.SS6hr:F9.SS6hr) %>%
  rowMeans() 

sel$mean.SS9hr = sel %>%
  select(B5.SS9hr:F10.SS9hr) %>%
  rowMeans() 

sel$mean.SS12hr = sel %>%
  select(B6.SS12hr:F11.SS12hr) %>%
  rowMeans() 

sel$mean.SD0hr = sel$mean.SS0hr

sel$mean.SD3hr = sel %>%
  select(B7.SD3hr:G2.SD3hr) %>%
  rowMeans()

sel$mean.SD6hr = sel %>%
  select(B8.SD6hr:G3.SD6hr) %>%
  rowMeans() 

sel$mean.SD9hr = sel %>%
  select(B9.SD9hr:G4.SD9hr) %>%
  rowMeans() 

sel$mean.SD12hr = sel %>%
  select(B10.SD12hr:G5.SD12hr) %>%
  rowMeans() 

sel.m = sel %>%
  select(id, geneSymbol, mean.SS0hr:mean.SD12hr) %>%
  gather(Sample_id,MeanExpression,mean.SS0hr:mean.SD12hr) %>%
  mutate(Condition=gsub("^mean.(SS|SD)([[:digit:]]+hr)","\\1",Sample_id),
         TP=gsub("^mean.(SS|SD)([[:digit:]]+)hr","\\2",Sample_id),
         TP=factor(TP, levels = c("0","3","6","9","12")))

sel.m %>%  
  ggplot(aes(x=TP, y=MeanExpression, group=Condition, color=Condition)) + 
  geom_path(size=2) +
  facet_wrap(~geneSymbol, ncol=4, scales = "free") +
  theme_bw(base_size = 40) +
  theme(legend.position = c(0.94,0.25), legend.title = element_text(size = 40),
        legend.text = element_text(size = 40), legend.key.height = unit(2,"cm")) +
  labs(y="Mean expression", x="Timepoint (in hours)")

ggsave("RESULTS/Longtitudinal/longtitudinal.antisense.genes.filtered.w.free.scale.pdf", height = 50, width = 100, units ="cm")


############ INTRON SENSE ############ 

sel <- read.table("RESULTS/DE.w.limmavoom/intron_level/DE_SS.vs.SD_intron.sense_filtered.txt", header = TRUE, sep = "\t")

sel$mean.SS0hr = sel %>%
  select(B2.0hr:F7.0hr) %>%
  rowMeans()

sel$mean.SS3hr = sel %>%
  select(B3.SS3hr:F8.SS3hr) %>%
  rowMeans()

sel$mean.SS6hr = sel %>%
  select(B4.SS6hr:F9.SS6hr) %>%
  rowMeans() 

sel$mean.SS9hr = sel %>%
  select(B5.SS9hr:F10.SS9hr) %>%
  rowMeans() 

sel$mean.SS12hr = sel %>%
  select(B6.SS12hr:F11.SS12hr) %>%
  rowMeans() 

sel$mean.SD0hr = sel$mean.SS0hr

sel$mean.SD3hr = sel %>%
  select(B7.SD3hr:G2.SD3hr) %>%
  rowMeans()

sel$mean.SD6hr = sel %>%
  select(B8.SD6hr:G3.SD6hr) %>%
  rowMeans() 

sel$mean.SD9hr = sel %>%
  select(B9.SD9hr:G4.SD9hr) %>%
  rowMeans() 

sel$mean.SD12hr = sel %>%
  select(B10.SD12hr:G5.SD12hr) %>%
  rowMeans() 

#Find the duplicated gene symboks and append a suffix
# n_occur = data.frame(table(sel$geneSymbol))
# n_occur[n_occur$Freq>1,]
# sel[sel$geneSymbol %in% n_occur$Var1[n_occur$Freq>1],]

# genesymbol.rle = rle(as.character(sel$geneSymbol))
# sel$geneSymbol.new = paste0(rep(genesymbol.rle$values, times = genesymbol.rle$lengths),"#",unlist(lapply(genesymbol.rle$lengths, seq_len))) 

# transform(sel, geneSymbol.new = ifelse(duplicated(geneSymbol) | duplicated(geneSymbol, fromLast = TRUE),
#           paste(geneSymbol, ave(geneSymbol, geneSymbol, FUN = seq_along), sep = "#"), geneSymbol))

# sel %>%
#   filter(geneSymbol=="Amigo2")%>%
#   group_by(geneSymbol)%>%
#   mutate(geneSymbol.new=paste0(geneSymbol,"#",row_number()))%>%
#   ungroup() %>%
#   as.data.frame()%>%
#   head

sel.m = sel %>%
  group_by(geneSymbol) %>%
  mutate(geneSymbol.new= paste0(geneSymbol,"#",row_number())) %>% 
  ungroup() %>%
  select(id, geneSymbol.new, mean.SS0hr:mean.SD12hr) %>%
  gather(Sample_id,MeanExpression,mean.SS0hr:mean.SD12hr) %>%
  mutate(Condition=gsub("^mean.(SS|SD)([[:digit:]]+hr)","\\1",Sample_id),
         TP=gsub("^mean.(SS|SD)([[:digit:]]+)hr","\\2",Sample_id),
         TP=factor(TP, levels = c("0","3","6","9","12")))

sel.m %>%  
  ggplot(aes(x=TP, y=MeanExpression, group=Condition, color=Condition)) + 
  geom_path(size=2) +
  facet_wrap(~geneSymbol.new, ncol=10, scales = "free") +
  theme_bw(base_size = 30) +
  theme(legend.position = c(0.94,0.001), legend.title = element_text(size = 30),
        legend.text = element_text(size = 30), legend.key.height = unit(2,"cm")) +
  labs(y="Mean expression", x="Timepoint (in hours)")

ggsave("RESULTS/Longtitudinal/longtitudinal.sense.introns.filtered.w.free.scale.pdf", height = 120, width = 120, units ="cm")



############ INTRON ANTISENSE ############ 

sel <- read.table("RESULTS/DE.w.limmavoom/intron_level/DE_SS.vs.SD_intron.antisense_filtered.txt", header = TRUE, sep = "\t")

sel$mean.SS0hr = sel %>%
  select(B2.0hr:F7.0hr) %>%
  rowMeans()

sel$mean.SS3hr = sel %>%
  select(B3.SS3hr:F8.SS3hr) %>%
  rowMeans()

sel$mean.SS6hr = sel %>%
  select(B4.SS6hr:F9.SS6hr) %>%
  rowMeans() 

sel$mean.SS9hr = sel %>%
  select(B5.SS9hr:F10.SS9hr) %>%
  rowMeans() 

sel$mean.SS12hr = sel %>%
  select(B6.SS12hr:F11.SS12hr) %>%
  rowMeans() 

sel$mean.SD0hr = sel$mean.SS0hr

sel$mean.SD3hr = sel %>%
  select(B7.SD3hr:G2.SD3hr) %>%
  rowMeans()

sel$mean.SD6hr = sel %>%
  select(B8.SD6hr:G3.SD6hr) %>%
  rowMeans() 

sel$mean.SD9hr = sel %>%
  select(B9.SD9hr:G4.SD9hr) %>%
  rowMeans() 

sel$mean.SD12hr = sel %>%
  select(B10.SD12hr:G5.SD12hr) %>%
  rowMeans() 

sel.m = sel %>%
  select(id, geneSymbol, mean.SS0hr:mean.SD12hr) %>%
  gather(Sample_id,MeanExpression,mean.SS0hr:mean.SD12hr) %>%
  mutate(Condition=gsub("^mean.(SS|SD)([[:digit:]]+hr)","\\1",Sample_id),
         TP=gsub("^mean.(SS|SD)([[:digit:]]+)hr","\\2",Sample_id),
         TP=factor(TP, levels = c("0","3","6","9","12")))

sel.m %>%  
  ggplot(aes(x=TP, y=MeanExpression, group=Condition, color=Condition)) + 
  geom_path(size=2) +
  facet_wrap(~geneSymbol, ncol=4, scales = "free") +
  theme_bw(base_size = 40) +
  theme(legend.position = c(0.9,0.2), legend.title = element_text(size = 40),
        legend.text = element_text(size = 40), legend.key.height = unit(2,"cm")) +
  labs(y="Mean expression", x="Timepoint (in hours)")

ggsave("RESULTS/Longtitudinal/longtitudinal.antisense.introns.filtered.w.free.scale.pdf", height = 50, width = 100, units ="cm")
