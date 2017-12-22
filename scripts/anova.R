library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(broom)
library(genefilter)

setwd("~/Box\ Sync/RNA-Seq/Sleep.Dep.dev/")

#GENE
sel <- read.table("RESULTS/DE.w.limmavoom/gene_level/DE_SS.vs.SD_sense_filtered.w.genetype.txt", sep = "\t", header = TRUE, check.names=FALSE)

sel.m = sel %>%
  select(id, geneSymbol, B3.SS3hr:G5.SD12hr) %>%
  gather(Sample_id, Expression, B3.SS3hr:G5.SD12hr) %>%
  mutate(Condition=gsub("^(B|C|D|E|F|G)([[:digit:]]+).(SS|SD)([[:digit:]]+hr)","\\3",Sample_id),
         TP=gsub("^(B|C|D|E|F|G)([[:digit:]]+).(SS|SD)([[:digit:]]+)hr","\\4",Sample_id),
         TP=factor(TP, levels = c("3","6","9","12")),
         Rep=gsub("^(B|C|D|E|F|G)([[:digit:]]+).(SS|SD\\[[:digit:]]+hr)","\\3",Sample_id))

sel.aov = aov(Expression~Condition*TP, sel.m)
summary(sel.aov)

#INTRON
sel <- read.table("RESULTS/DE.w.limmavoom/intron_level/DE_SS.vs.SD_intron.sense_filtered.txt", sep = "\t", header = TRUE, check.names=FALSE)

sel.m = sel %>%
  select(id, geneSymbol, B3.SS3hr:G5.SD12hr) %>%
  gather(Sample_id, Expression, B3.SS3hr:G5.SD12hr) %>%
  mutate(Condition=gsub("^(B|C|D|E|F|G)([[:digit:]]+).(SS|SD)([[:digit:]]+hr)","\\3",Sample_id),
         TP=gsub("^(B|C|D|E|F|G)([[:digit:]]+).(SS|SD)([[:digit:]]+hr)","\\4",Sample_id),
         TP=factor(TP, levels = c("3hr","6hr","9hr","12hr")),
         Rep=gsub("^(B|C|D|E|F|G)([[:digit:]]+).(SS|SD\\[[:digit:]]+hr)","\\3",Sample_id))

sel.aov = aov(Expression~Condition*TP, sel.m)
summary(sel.aov)



#GENE with 0hr
sel <- read.table("RESULTS/DE.w.limmavoom/gene_level/DE_SS.vs.SD_sense_filtered.w.genetype.RENAMED.txt", sep = "\t", header = TRUE, check.names=FALSE)

sel.m = sel %>%
  select(id, geneSymbol, B2.CTL0hr:G5.SD12hr) %>%
  gather(Sample_id, Expression, B2.CTL0hr:G5.SD12hr) %>%
  mutate(Condition=gsub("^(B|C|D|E|F|G)([[:digit:]]+).(SS|SD|CTL)([[:digit:]]+hr)","\\3",Sample_id),
         TP=gsub("^(B|C|D|E|F|G)([[:digit:]]+).(SS|SD|CTL)([[:digit:]]+)hr","\\4",Sample_id),
         TP=factor(TP, levels = c("0","3","6","9","12")))
  # mutate(Condition=gsub("^(B|C|D|E|F|G)([[:digit:]]+).(SS|SD|.{0})([[:digit:]]+hr)","\\3",Sample_id),
  #      TP=gsub("^(B|C|D|E|F|G)([[:digit:]]+).(SS|SD|.{0})([[:digit:]]+)hr","\\4",Sample_id),
  #      TP=factor(TP, levels = c("0","3","6","9","12")))

sel.aov = aov(Expression~Condition*TP, sel.m)
summary(sel.aov)



#INTRON with 0hr
sel <- read.table("RESULTS/DE.w.limmavoom/intron_level/DE_SS.vs.SD_intron_sense.RENAMED.txt", sep = "\t", header = TRUE, check.names=FALSE)

sel.m = sel %>%
  select(id, geneSymbol, B2.CTL0hr:G5.SD12hr) %>%
  gather(Sample_id, Expression, B2.CTL0hr:G5.SD12hr) %>%
  mutate(Condition=gsub("^(B|C|D|E|F|G)([[:digit:]]+).(SS|SD|CTL)([[:digit:]]+hr)","\\3",Sample_id),
         TP=gsub("^(B|C|D|E|F|G)([[:digit:]]+).(SS|SD|CTL)([[:digit:]]+)hr","\\4",Sample_id),
         TP=factor(TP, levels = c("0","3","6","9","12")))

sel.aov = aov(Expression~Condition*TP, sel.m)
summary(sel.aov)
