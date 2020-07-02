#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)

## Arguments:
# 1. Home directory
# 2. Dataset

library(alphaOutlier)
library(plyr)
library(readr)
library("GenBinomApps")

HOMEDIR <- args[1]
# Read in significant filtered associations data: selected.sign.anovas.df and selected.all.IC50s.df
load(paste0(HOMEDIR, '/results/', args[2], '/sign_anova_out.RData'))
load(paste0(HOMEDIR, '/results/', args[2], '/high_var_results.RData'))

#### Neyman-Pearson Application
assoc.numbers.all <- selected.all.IC50s.df %>%
  mutate(tissue.assoc = paste0(tissue, '-', association)) %>%
  select(tissue.assoc, log10.dr,mut.status) %>%
  group_by(tissue.assoc) %>%
  arrange(-log10.dr) %>%
  dplyr::mutate(index=1:n()) %>%
  arrange(tissue.assoc)

## Format data for fast sampling
selected.all.IC50s.to.sample.all <- selected.all.IC50s.df %>%
  mutate(tissue.assoc = paste0(tissue, '-', association)) %>%
  group_by(tissue.assoc) %>%
  dplyr::mutate(nr.mut=sum(mut.status == 1),
         total.cell = length(mut.status)) %>%
  ungroup() %>%
  select(tissue.assoc, log10.dr, nr.mut,total.cell,mut.status) %>%
  nest(nest.log10dr = c(log10.dr),
       nest.mut = c(mut.status))

# significant level
alpha <- 0.15

# outliers from resistance.analysis
high.var.associations.filter <- high.var.associations[high.var.associations$outlier.q.value<alpha,]
outliers.resistance.analysis <- unique(high.var.associations.filter %>%
                                         select(association,tissue,outlier.p.value,outlier.q.value))


# 2: define normal distribution values: only IC50s of mut=1
n.cells <- nrow(assoc.numbers.all[assoc.numbers.all$mut.status==1,])
sd.cells <- sd(assoc.numbers.all[assoc.numbers.all$mut.status==1,]$log10.dr)
mean.cells <- mean(assoc.numbers.all[assoc.numbers.all$mut.status==1,]$log10.dr)
outlier.binomial <- matrix(nrow = nrow(selected.all.IC50s.to.sample.all),
                           ncol = 2)
rownames(outlier.binomial) <- selected.all.IC50s.to.sample.all$tissue.assoc

for (i in 1:nrow(selected.all.IC50s.to.sample.all)) {
  up <- qnorm(1-alpha,mean.cells,sd.cells)
  mean.mutated.i <- (data.frame(selected.all.IC50s.to.sample.all[i,]$nest.log10dr)[
    unlist(selected.all.IC50s.to.sample.all[i,]$nest.mut)==1,])
  outlier.binomial.aux<- (mean.mutated.i>=up)
  outlier.binomial[i,1] <- any(outlier.binomial.aux)
  outlier.binomial[i,2] <- length(which(outlier.binomial.aux))
}

high.var.associations.new <- map_dfr(rownames(outlier.binomial), function(x) {
  attrib <- strsplit(x, '-')[[1]]
  selected.sign.anovas.df[selected.sign.anovas.df$association == paste(attrib[-1], collapse = '-') &
                            selected.sign.anovas.df$tissue ==attrib[1],]
})
high.var.associations.new$outlier <- outlier.binomial[,1]
high.var.associations.new$nr.outlier <- outlier.binomial[,2]
outliers.Neyman.Pearson <- unique(high.var.associations.new %>%
                           filter(outlier == TRUE) %>%
                           mutate(tissue.assoc = paste0(tissue, '-', association)) %>%
                           select(tissue.assoc,nr.outlier))
