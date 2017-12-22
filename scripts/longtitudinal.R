library(ggplot2)
library(dplyr)
library(tidyr)

setwd("~/Box\ Sync/RNA-Seq/Sleep.Dep.dev/")

########### SENSE GENES #############
sel <- read.table("RESULTS/DE.w.limmavoom/selected.genes/sense.genes.selected.txt", header = TRUE, sep = "\t")

order.sel = sel[,c("id",
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

sel = order.sel

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
  facet_wrap(~geneSymbol, ncol=3, scales = "free") +
  theme_bw(base_size = 50) +
  theme(legend.position = c(0.94,0.3), legend.title = element_text(size = 40),
        legend.text = element_text(size = 40), legend.key.height = unit(2,"cm")) +
labs(y="Mean expression", x="Timepoint (in hours)")

ggsave("RESULTS/Longtitudinal/longtitudinal.sense.genes.selected.pdf", height = 30, width = 100, units ="cm")



########### ANTISENSE GENES #############

sel <- read.table("RESULTS/DE.w.limmavoom/selected.genes/antisense.genes.selected.txt", header = TRUE, sep = "\t")

order.sel = sel[,c("id",
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

sel = order.sel

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
  facet_wrap(~geneSymbol, ncol=3, scales = "free") +
  theme_bw(base_size = 50) +
  theme(legend.position = c(0.96,0.85), legend.title = element_text(size = 40),
        legend.text = element_text(size = 40), legend.key.height = unit(2,"cm")) +
  labs(y="Mean expression", x="Timepoint (in hours)")

ggsave("RESULTS/Longtitudinal/longtitudinal.antisense.genes.selected.pdf", height = 30, width = 100, units ="cm")



############ INTRONS SENSE ##############

sel <- read.delim("RESULTS/DE.w.limmavoom/selected.genes/sense.introns.selected.txt", header = TRUE, sep = "\t")

order.sel = sel[,c("id",
                   "B2.0hr","B11.0hr","C10.0hr","D9.0hr","E8.0hr","F7.0hr",
                   "B3.SS3hr","C2.SS3hr","C11.SS3hr","D10.SS3hr","E9.SS3hr","F8.SS3hr",
                   "B4.SS6hr","C3.SS6hr","D2.SS6hr","D11.SS6hr","E10.SS6hr","F9.SS6hr",
                   "B5.SS9hr","C4.SS9hr","D3.SS9hr","E2.SS9hr","F10.SS9hr",
                   "B6.SS12hr","C5.SS12hr","D4.SS12hr","E3.SS12hr","F2.SS12hr","F11.SS12hr",
                   "B7.SD3hr","C6.SD3hr","D5.SD3hr","E4.SD3hr","F3.SD3hr","G2.SD3hr",
                   "B8.SD6hr","C7.SD6hr","D6.SD6hr","E5.SD6hr","F4.SD6hr","G3.SD6hr",
                   "B9.SD9hr","C8.SD9hr","D7.SD9hr","E6.SD9hr","F5.SD9hr","G4.SD9hr",
                   "B10.SD12hr","C9.SD12hr","E7.SD12hr","F6.SD12hr","G5.SD12hr",
                   "geneid","geneSymbol","logFC","LimmaVoom.pvalue","BH.qvalue")]

sel = order.sel

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
  facet_wrap(~geneSymbol, ncol=5, scales = "free") +
  theme_bw(base_size = 50) +
  theme(legend.position = c(0.7,0.1), legend.title = element_text(size = 40),
        legend.text = element_text(size = 40), legend.key.height = unit(2,"cm")) +
  labs(y="Mean expression", x="Timepoint (in hours)")

ggsave("RESULTS/Longtitudinal/longtitudinal.sense.introns.selected.pdf", height = 100, width = 100, units ="cm")




############ ANTISENSE INTRONS #################
sel <- read.delim("RESULTS/DE.w.limmavoom/selected.genes/antisense.introns.selected.txt", header = TRUE, sep = "\t")

order.sel = sel[,c("id",
                   "B2.0hr","B11.0hr","C10.0hr","D9.0hr","E8.0hr","F7.0hr",
                   "B3.SS3hr","C2.SS3hr","C11.SS3hr","D10.SS3hr","E9.SS3hr","F8.SS3hr",
                   "B4.SS6hr","C3.SS6hr","D2.SS6hr","D11.SS6hr","E10.SS6hr","F9.SS6hr",
                   "B5.SS9hr","C4.SS9hr","D3.SS9hr","E2.SS9hr","F10.SS9hr",
                   "B6.SS12hr","C5.SS12hr","D4.SS12hr","E3.SS12hr","F2.SS12hr","F11.SS12hr",
                   "B7.SD3hr","C6.SD3hr","D5.SD3hr","E4.SD3hr","F3.SD3hr","G2.SD3hr",
                   "B8.SD6hr","C7.SD6hr","D6.SD6hr","E5.SD6hr","F4.SD6hr","G3.SD6hr",
                   "B9.SD9hr","C8.SD9hr","D7.SD9hr","E6.SD9hr","F5.SD9hr","G4.SD9hr",
                   "B10.SD12hr","C9.SD12hr","E7.SD12hr","F6.SD12hr","G5.SD12hr",
                   "geneid","geneSymbol","logFC","LimmaVoom.pvalue","BH.qvalue")]

sel = order.sel

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
  facet_wrap(~geneSymbol, ncol=5, scales = "free") +
  theme_bw(base_size = 50) +
  theme(legend.position = c(0.9,0.1), legend.title = element_text(size = 40),
        legend.text = element_text(size = 40), legend.key.height = unit(2,"cm")) +
  labs(y="Mean expression", x="Timepoint (in hours)")

ggsave("RESULTS/Longtitudinal/longtitudinal.antisense.introns.selected.pdf", height = 100, width = 100, units ="cm")
