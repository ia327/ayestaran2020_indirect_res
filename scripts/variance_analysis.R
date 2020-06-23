#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)

# args should be a vector of length 4-5:
# 1. Home directory of the repository
# 2. "gdsc" or "ctrp"
# 3. ANOVA filter result directory
# 4. Output directory
# 5. (optional) Permutation ID, also deactivates plots


#### This script takes as input significant drug-CFE associations and
#### it performs a analysis of change in variance to detect outliers.
#### The results of associations in which there are outliers between the mutant
#### cell lines are saved.

suppressMessages(library(tidyverse))
library(parallel)

HOMEDIR <- args[1]
INDIR <- args[3]
OUTDIR <- args[4]

######################### ANALYSIS #########################
stopifnot(file.exists(paste0(INDIR, 'sign_anova_out.RData')))

# Read in significant filtered associations data: selected.sign.anovas.df and selected.all.IC50s.df
load(paste0(INDIR, 'sign_anova_out.RData'))

stopifnot(nrow(selected.sign.anovas.df) > 0)

#### Variance change analysis

# Select mut cell lines, and sort and add index by decreasing IC50
assoc.numbers <- selected.all.IC50s.df %>%
  filter(mut.status == 1) %>%
  mutate(tissue.assoc = paste0(tissue, '-', association)) %>%
  select(tissue.assoc, log10.dr) %>%
  group_by(tissue.assoc) %>%
  arrange(-log10.dr) %>%
  mutate(index=1:n()) %>%
  arrange(tissue.assoc)

# Observed variances will be stored here
variances <- matrix(nrow = length(unique(assoc.numbers$tissue.assoc)),
                    # ncol =  max(assoc.numbers$index))
                    ncol = 5+1) # Max n outliers to consider + 1
rownames(variances) <- unique(assoc.numbers$tissue.assoc)

# First column: initial SD of the entire mut subpopulation (sigma_0)
variances[,1] <- assoc.numbers %>%
  summarise(var1=sd(log10.dr)) %>%
  select(var1) %>%
  unlist() %>% unname()

all.assocs.ref <- ungroup(unique(assoc.numbers['tissue.assoc']))

# For subsequent columns, remove the top i IC50s and get SD
# variances[,2:ncol(variances)] <- map_dfc(2:ncol(variances), ~assoc.numbers %>%
#                                            filter(index >= .x) %>%
#                                            summarise(var1=sd(log10.dr)) %>%
#                                            right_join(all.assocs.ref,
#                                                       by = "tissue.assoc") %>% # Makes sure that the length is correct
#                                            select(var1) %>%
#                                            unlist() %>% unname()) %>%
#   as.matrix()

variances[,2:ncol(variances)] <- map_dfc(2:ncol(variances), ~assoc.numbers %>%
                                           ungroup() %>%
                                           filter(index >= .x) %>%
                                           split(.$tissue.assoc) %>%
                                           map_dbl(~sd(.x$log10.dr)) %>%
                                           data.frame(var1=.) %>%
                                           rownames_to_column('tissue.assoc') %>%
                                           right_join(all.assocs.ref,
                                                      by = "tissue.assoc") %>% # Makes sure that the length is correct
                                           select(var1) %>%
                                           unlist() %>% unname()) %>%
  as.matrix()


# Calculate difference with respect to sigma_0 
diff.variances <- (variances[,2:ncol(variances)] - variances[,1])
dim(diff.variances) <- dim(variances) + c(0, -1)
rownames(diff.variances) <- rownames(variances)

## Format data for fast sampling
selected.all.IC50s.to.sample <- selected.all.IC50s.df %>%
  mutate(tissue.assoc = paste0(tissue, '-', association)) %>%
  group_by(tissue.assoc) %>%
  mutate(to.sample=sum(mut.status == 1)) %>%
  ungroup() %>%
  select(tissue.assoc, log10.dr, to.sample) %>%
  nest(nest.log10dr = c(log10.dr))

######## Null bootstrap model for change in variance
# Get random values for IC50s  SAMPLING WITH REPLACEMENT
null.model.var.change <- function(j){
  if (j %% 50 == 0){
    cat('\tRep#', j, '\033[0K\r')
  }
  
  null.assoc.numbers <- selected.all.IC50s.to.sample %>%
    mutate(nest.log10dr = map2(nest.log10dr, to.sample, sample_n, replace=TRUE)) %>%
    unnest(cols = nest.log10dr) %>%
    group_by(tissue.assoc) %>%
    arrange(-log10.dr) %>%
    select(-to.sample) %>%
    mutate(index=1:n()) %>%
    arrange(tissue.assoc)
  
  null.variances <- matrix(nrow = length(unique(null.assoc.numbers$tissue.assoc)),
                           # ncol =  max(null.assoc.numbers$index))
                           ncol=5+1)
  rownames(null.variances) <- unique(null.assoc.numbers$tissue.assoc)
  
  null.variances[,1] <- null.assoc.numbers %>%
    summarise(var1=sd(log10.dr)) %>%
    select(var1) %>%
    unlist() %>% unname()
  
  # null.variances[,2:ncol(null.variances)] <- map_dfc(2:ncol(null.variances),
  #                                                    ~null.assoc.numbers %>%
  #                                                      filter(index >= .x) %>%
  #                                                      summarise(var1=sd(log10.dr)) %>%
  #                                                      right_join(all.assocs.ref,
  #                                                                 by = "tissue.assoc") %>% # Makes sure that the length is correct
  #                                                      select(var1) %>%
  #                                                      unlist() %>% unname()) %>%
  #   as.matrix()
  
  null.variances[,2:ncol(null.variances)] <- map_dfc(2:ncol(null.variances),
                                                     ~null.assoc.numbers %>%
                                                       ungroup() %>%
                                                       filter(index >= .x) %>%
                                                       split(.$tissue.assoc) %>%
                                                       map_dbl(~sd(.x$log10.dr)) %>%
                                                       data.frame(var1=.) %>%
                                                       rownames_to_column('tissue.assoc') %>%
                                                       right_join(all.assocs.ref,
                                                                  by = "tissue.assoc") %>% # Makes sure that the length is correct
                                                       select(var1) %>%
                                                       unlist() %>% unname()) %>%
    as.matrix()
  
  
  
  null.diff.variances <- (null.variances[,2:ncol(null.variances)] - null.variances[,1])
  dim(null.diff.variances) <- dim(null.variances) + c(0, -1)
  rownames(null.diff.variances) <- rownames(null.variances)
  
  
  return(null.diff.variances)
}

## Replicate the null model many times and get its statistics
reps <- 10000
cat('Running', reps, ' bootstrap iterations: \n')
null.model.res <- mclapply(1:reps, function (ii) null.model.var.change(ii),
                           mc.cores = 8)
cat('\nBootstrap done\n')

################################################################

pval.null.model <- matrix(nrow = nrow(diff.variances),
                          ncol = ncol(diff.variances))
rownames(pval.null.model) <- rownames(diff.variances)

expected.var.diff <- matrix(nrow = nrow(diff.variances),
                            ncol = ncol(diff.variances))
rownames(expected.var.diff) <- rownames(diff.variances)

# pdf(paste0(HOMEDIR, '/plots/', args[2], '/all_sub_hist.pdf'))
# par(mfrow=c(2,2))
for (i in 1:nrow(pval.null.model)){
  this.assoc <- rownames(pval.null.model)[i]
  size.this.assoc <- sum(assoc.numbers$tissue.assoc == this.assoc)
  for (j in 1:ncol(pval.null.model)){
    if (j <= round(size.this.assoc / 2)){
      sub <- map_dbl(null.model.res, ~.x[i,j])
      if (sum(is.na(sub)) != length(sub)){
        # The distribution of sub follows a normal distribution, approx
        # hist(sub, main = paste0(rownames(diff.variances)[i],' - ', j, ' out'),
        # breaks = 50, xlab = bquote(paste(Delta, sigma[.(j)]^2)))
        # abline(v = diff.variances[i, j], col = 'red', lwd = 2)
        # abline(v = median(sub), col = 'blue', lwd = 2)
        pval.null.model[i, j] <- (sum(sub < diff.variances[i, j])+1) / (length(sub)+1)
        expected.var.diff[i, j] <- median(sub)
      }
    }
  }
}
# invisible(dev.off())

# qval.null.model <- t(apply(pval.null.model, 1, p.adjust, method = 'BH'))
# qval.null.model <- pval.null.model
# qval.null.model[,] <- p.adjust(pval.null.model, method = 'BH')

## Hierarchical FDR control
# level 0: check which cells with sensitive biomarker p.value<alpha
alpha <- 0.15

# level 1: define child hypothesis based on significant parent hypothesis 
## Select outliers when the change in variance is significant
high.var.change <- which(pval.null.model <= 1, # & pval.null.model < 0.05,
                         arr.ind = TRUE)

high.var.associations <- map_dfr(rownames(high.var.change), function(x) {
  attrib <- strsplit(x, '-')[[1]]
  selected.sign.anovas.df[selected.sign.anovas.df$association == paste(attrib[-1], collapse = '-') &
                            selected.sign.anovas.df$tissue ==attrib[1],]
}
)

qval.corrected <- p.adjust(pval.null.model[high.var.change], method = "BH")
sign.children.hypothesis <- qval.corrected < alpha

# calculate FDR of full tree 
n.families.tested <- sum(selected.sign.anovas.df$p.value<alpha)
n.children.sign <- sum(qval.corrected < alpha)
n.tree.discoveries <- n.families.tested + n.children.sign
fdr.tree <- min(1,(n.tree.discoveries+n.families.tested)/(n.tree.discoveries+1)*alpha)

save(list = 'fdr.tree',
     file = paste0(OUTDIR, 'fdr.tree.value.RData'))

high.var.associations$n.outliers <- high.var.change[,2]
high.var.associations$outlier.p.value <- pval.null.model[high.var.change]
high.var.associations$outlier.q.value <- qval.corrected
high.var.associations$initial.var <- variances[high.var.change[,1], 1]
high.var.associations$var.diff <- diff.variances[high.var.change]
high.var.associations$expected.var.diff <- expected.var.diff[high.var.change]
high.var.associations <- mutate(high.var.associations,
                                effect.size=-(var.diff-expected.var.diff)/(initial.var-expected.var.diff))
rownames(high.var.associations) <- paste0(rownames(high.var.change), '_',
                                          high.var.associations$n.outliers, '_out')


## Some sorting
high.var.associations <- high.var.associations[order(high.var.associations$drug),]
high.var.associations <- high.var.associations[order(high.var.associations$alteration),]
high.var.associations <- high.var.associations[order(high.var.associations$tissue),]


## Save the associations with 1 or more outliers as an RData object
save(list = 'high.var.associations',
     file = paste0(OUTDIR, 'high_var_results.RData'))
