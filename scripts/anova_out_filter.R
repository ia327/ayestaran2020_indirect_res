#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)

## Arguments:
# 1. Home directory
# 2. Dataset
# 3. ANOVA result directory
# 4. Output directory
# 5. (optional) Do permutation? also deactivates plots

#### This script takes as input the results of tissue-specific ANOVA and
#### saves significant and effective enough drug-gene associations (sign_anova_out)

## Filters used:
## 1: P value < 0.001, FDR < 25%, effect size > 1 towards sensitivity
## 2: Associations contain at least one driver gene (according to COSMIC Gene census,
##    previously according to GDSC driver gene list)
## 3: At least 50% (rounded) of mutant cell lines must have an IC50 below the
##    maximum concentration of the drug
## 4: A minimum number of mutant cell lines (n = 4)

## Volcano plots and barcharts are drawn for the initial ANOVA results
## and after the filtering process
## Stepwise volcano plot visulazations are commented out

### Data used:
## DRUG_NAME, DRUG_ID, TCGA_DESC, COSMIC_ID,
## PUTATIVE_TARGET, LN_IC50, MAX_CONC

suppressMessages(library(tidyverse))
library(ggbeeswarm)
library(ggrepel)
suppressMessages(library(scales))

HOMEDIR <- args[1]
DATADIR <- paste0(HOMEDIR, '/data/combined_for_pipeline/')
ANOVAOUTDIR <- args[3]
OUTDIR <- args[4]
do_permutation <- !is.na(args[5]) & args[5] == 'permute'

load(paste0(DATADIR, args[2], '_dr_data.RData'))
load(paste0(DATADIR, 'utils.RData'))
load(paste0(DATADIR, 'bems_tidy.RData'))
load(paste0(DATADIR, 'cn_annotation.RData'))
load(paste0(DATADIR, 'cancer_gene_lists.RData'))


dataset <- switch(args[2],
                  gdsc = all.gdsc,
                  ctrp = all.ctrp)

all.bems.tidy <- switch(args[2],
                        gdsc = all.bems.tidy.gdsc,
                        ctrp = all.bems.tidy.ctrp)
#### Plotting functions
# function for transforming y-axis
reverselog_trans <- function (base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^-x
  trans_new(paste0("log-", format(base)), trans, inv, log_breaks(base = base),
            domain = c(1e-100, Inf))
}

# function for relabeling the y-axis
scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scientific_format()(x)))
}



plot.volcano.gg <- function(data, label.thres = 0.01, main = '', fixed.x = FALSE,
                            filter.sig=FALSE){
  
  drug.label <- sapply(data$drug,
                       function(x) unique(dataset$DRUG_NAME[which(dataset$DRUG_ID== x)]))
  labels <- paste0(drug.label, '\n', data$alteration)
  labels[data$corrected.p.value...q.value. > label.thres] <- ''
  col.vector <- tissue.cols[data$tissue]
  
  all <- ggplot(data=data,
                aes(x=effect, y=p.value,
                    colour=tissue, size = n.samples)) +
    scale_y_continuous(trans=reverselog_trans(10),
                       label=scientific_10,
                       limits=c(0.3, min(data$p.value)/10)) +
    geom_point(aes(size=(n.samples+10)/10)) +
    geom_vline(xintercept=0, lwd=0.3) +
    geom_vline(xintercept=c(-1,1), lwd=0.3, linetype=2) +
    geom_hline(yintercept=max.pval, lwd=0.3, linetype=2) +
    annotate("text", x=-5.8 + 0.5, y=max.pval-0.0003, size=3.5,
             label=paste('P-value = ', max.pval, sep=""), colour="black") +
    annotate("text", x=min.effect.size + 1.5,
             y=min(data$p.value)/9, size=3.5,
             label=paste('Effect size of ', min.effect.size, sep=""), colour="black") +
    xlab("Drug effect") +
    ylab("P-value") +
    ggtitle(main) +
    theme_bw() +
    theme(plot.title=element_text(size=14, face="bold", hjust = 0.5)) +
    scale_colour_manual('Tissue', values = tissue.cols,
                        breaks = names(tissue.cols)[-length(tissue.cols)]) +
    scale_size('N. altered samples',
               breaks = c((40+10)/10, (20+10)/10, (3+10)/10),
               labels = c('40', '20', '3'))
  if (label.thres > 0){
    all + geom_text_repel(mapping=aes(x=effect,
                                      y=p.value,
                                      label=labels), colour="gray40", size=3,
                          box.padding = unit(0.5, 'lines'))
    
  } else {
    if (fixed.x){
      all +  scale_x_continuous(breaks =c(-6, -4, -2, 0),
                                limits=c(-6, 0))
    } else {
      all + scale_x_continuous(breaks =seq(floor(min(data$effect)),
                                           ceiling(max(data$effect)+2), 2),
                               limits=c(floor(min(data$effect)),
                                        ceiling(max(data$effect)+2)))
    }
  }
}

barchart.volcano <- function(data, fixed.y = TRUE){
  tissue.table <- as.data.frame(table(data$tissue))
  all <- ggplot(data = tissue.table, aes(x = Var1, y=Freq, fill = Var1)) +
    geom_bar(colour = 'black', stat = 'identity') +
    xlab('') +
    ylab('Number of associations') +
    theme_bw() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x = element_text(angle = 60, vjust = 0.5)) +
    scale_fill_manual('Tissue', values = tissue.cols,
                      breaks = names(tissue.cols),
                      guide = FALSE)
  if (fixed.y){
    all #+ ylim(0,15)
  } else {
    all
  }
}


######################### ANALYSIS #########################

# Read in result data from ANOVA
# Use only tissues with at least 15 cell lines
used.tissues <- dataset %>%
  select(TCGA_DESC, COSMIC_ID) %>%
  unique() %>%
  count(TCGA_DESC) %>%
  filter(n >= 15 & TCGA_DESC %in% all.tissues) %>%
  select(TCGA_DESC) %>% unlist() %>% unname() %>% as.character()

all.anovas.df <- map_dfr(paste0(ANOVAOUTDIR, used.tissues, '_anova_out.csv'),
                         read_csv, col_type = cols(),
                         .id = 'id') %>%
  rename_all(make.names) %>%
  mutate(drug=map_chr(str_split(association, ':'), 1),
         alteration=str_replace(association, paste0(drug,':'), ''),
         tissue=used.tissues[as.numeric(id)]) %>%
  select(-id) %>%
  as.data.frame()

# Select only significant results with a big enough effect
max.pval <- 0.001
min.effect.size <- -1 # Min in absolute terms

sign.anovas.df <- all.anovas.df %>%
  filter(p.value < max.pval &
           effect < min.effect.size)
# sign.anovas.df <- all.anovas.df %>%
#   filter(p.value < max.pval &
#            effect < min.effect.size)

# Filter out sensitivity markers that do not contain any known driver gene (GDSC)
sign.anovas.genes <- lapply(1:nrow(sign.anovas.df), function(x){
  alt <- sign.anovas.df$alteration[x]
  alt.type <- ifelse(grepl('_mut', alt), 'mut', 'cna')
  if (alt.type == 'mut'){
    res <- strsplit(alt, '_')[[1]][1]
  }
  if (alt.type == 'cna'){
    region <- strsplit(alt, ' ')[[1]][1]
    region <- strsplit(region, ':')[[1]][2]
    res <- as.character(RACS_to_gene[RACS_to_gene$Identifier == region, 'Contained genes'])
    res <- strsplit(res, ',')[[1]]
  }
  return(res)
}
)



driver.mask <- sapply(1:nrow(sign.anovas.df), function(x) {
  is.driver <- sum(sign.anovas.genes[[x]] %in% cancer.genes.gdsc$Gene) >= 1
  gene.tissues <- lapply(sign.anovas.genes[[x]],
                         function(y) unlist(cancer.genes.gdsc[cancer.genes.gdsc$Gene == y,
                                                              'Cancer Type']))
  correct.tissue <- sapply(gene.tissues, function(y) sign.anovas.df$tissue[x] %in% y)
  return(is.driver & any(correct.tissue))
}
)

selected.sign.anovas.df <- sign.anovas.df[driver.mask,]

stopifnot(nrow(selected.sign.anovas.df) > 0)
## We want to see the statistical properties of the IC50 values
## in the mutant populations

all.IC50s <- apply(selected.sign.anovas.df, 1, function(x) {
  simple.bem <- all.bems.tidy %>%
    filter(Tissue == as.character(x['tissue'])) %>%
    select(-Tissue)
  
  # Get Drug response and mutation status data
  selected.dr <- dataset %>%
    filter(DRUG_ID == as.character(x['drug'])) %>%
    right_join(simple.bem, by = 'COSMIC_ID') %>%
    filter(alteration == as.character(x['alteration'])) %>%
    arrange(COSMIC_ID) %>%
    unique() %>%
    mutate(log10.dr = log10(exp(LN_IC50)))
  
  # Get max conc value (if there is a single different value, treat it as and error)
  max.conc <- log10(as.numeric(names(sort(table(selected.dr$MAX_CONC),
                                          decreasing = TRUE)))[1])
  
  res <- selected.dr %>%
    mutate(tissue=as.character(x['tissue']), 
           log10.max.conc=max.conc,
           gene=as.character(x['alteration'])) %>%
    select(tissue, COSMIC_ID, DRUG_ID, PUTATIVE_TARGET,
           gene, mut.status, log10.dr, log10.max.conc) %>%
    mutate_if(is.character, as.factor) %>%
    mutate(association=as.character(x['association']),
           COSMIC_ID = as.character(COSMIC_ID))
  
  if (do_permutation){
    perm.dr <- read.csv(paste0(ANOVAOUTDIR, x['tissue'],'_dr.csv'),
                        row.names = 1, check.names = FALSE)
    drug_target <- unique(na.omit(res$PUTATIVE_TARGET))
    drug_target <- ifelse(length(drug_target) == 0, NA,
                          drug_target)
    res <- tibble(LN_IC50_perm=perm.dr[, as.character(x['drug'])],
                  COSMIC_ID=rownames(perm.dr)) %>%
      full_join(res, by = 'COSMIC_ID') %>%
      mutate(log10.dr=log10(exp(LN_IC50_perm)),
             DRUG_ID = as.character(x['drug']),
             PUTATIVE_TARGET = drug_target) %>%
      select(-LN_IC50_perm)
  }
  return(res[!is.na(res$log10.dr),])
})

all.IC50s.df <- do.call('rbind', all.IC50s)


## More filters

# Select only associations with at least 4 cell lines
min.n.cl <- 4

## Only keep associations with at least 50% of the mutant cell lines
## with a IC50 below the maximum concentration

max.frac.below <- 0.5
filters.to.keep <- all.IC50s.df %>%
  filter(mut.status == 1) %>%
  group_by(association, tissue) %>%
  summarise(n.cell.lines= n(),
            below.max = sum(log10.dr <= log10.max.conc),
            frac.below.max = below.max/n.cell.lines) %>%
  ungroup() %>%
  filter(frac.below.max >= max.frac.below &
           n.cell.lines >= min.n.cl) %>%
  select(tissue, association) %>%
  mutate_at(c('tissue', 'association'), as.character)


selected.sign.anovas.df <- selected.sign.anovas.df %>%
  mutate_at('association', as.character) %>%
  right_join(filters.to.keep, by = c("association", "tissue"))

selected.all.IC50s.df <- all.IC50s.df %>%
  mutate_at('tissue', as.character) %>%
  right_join(filters.to.keep, by = c("association", "tissue"))


if (!do_permutation) {
  pdf(paste0(HOMEDIR, '/plots/', args[2], '/volcano_filtered.pdf'), width = 8, height = 6)
  plot.volcano.gg(selected.sign.anovas.df, label.thres = 0,
                  main = 'Filtered', fixed.x = TRUE)
  invisible(dev.off())
  pdf(paste0(HOMEDIR, '/plots/', args[2], '/barchart_filtered.pdf'), width = 3.6, height = 2.4)
  barchart.volcano(selected.sign.anovas.df)
  invisible(dev.off())
}

save(list = c('selected.sign.anovas.df', 'selected.all.IC50s.df'),
     file = paste0(OUTDIR, 'sign_anova_out.RData'))




#### Export the resulting list #####
# Auxiliary function to rename CNAs
rename.cna <- function(x) {
  if (grepl('mutation', x)){
    return(x)
  }
  region <- strsplit(x, ':')[[1]][2]
  region2 <- strsplit(region, ' ')[[1]][1]
  return(unique(paste0(RACS_to_gene$locus[RACS_to_gene$Identifier == region2], ' ',
                       RACS_in_cell_lines$Alteration.Type[RACS_in_cell_lines$Region.identifier == region],
                       ' ', strsplit(region, ' ')[[1]][2])))
}

save.df <- selected.sign.anovas.df[order(selected.sign.anovas.df$corrected.p.value...q.value.),]
save.df$DRUG_ID <- as.character(save.df$drug)
save.df <- dataset %>%
  filter(DRUG_ID %in% save.df$DRUG_ID) %>%
  select(c(DRUG_ID, DRUG_NAME, PUTATIVE_TARGET)) %>%
  unique() %>%
  right_join(save.df, by = "DRUG_ID") %>%
  mutate_at('alteration', function(x) gsub('_mut', ' mutation', x)) %>%
  mutate_at('alteration', function(x) sapply(x, rename.cna)) %>%
  select(c(tissue, DRUG_NAME, PUTATIVE_TARGET, alteration,
           n.samples, p.value, corrected.p.value...q.value.,
           effect))

names(save.df) <- c('Tissue', 'Drug', 'Drug Target', 'Sensitivity marker',
                    'N. of sensitive samples', 'P-value', 'FDR',
                    'Effect size')
write.csv(save.df, file = paste0(OUTDIR, 'sign_anova_out.csv'),
          row.names = FALSE)


if (!do_permutation) {
  ########################## OTHER PLOTS ##########################
  ### Volcano plot for sign.anovas.df
  all.anovas.df.to.plot <- all.anovas.df
  all.anovas.df.to.plot$tissue[all.anovas.df.to.plot$p.value > max.pval |
                                 # all.anovas.df.to.plot$corrected.p.value...q.value. > max.fdr |
                                 abs(all.anovas.df.to.plot$effect) < abs(min.effect.size)] <- 'Not significant'
  all.anovas.unfiltered <- all.anovas.df.to.plot
  all.anovas.df.to.plot <- all.anovas.df.to.plot[all.anovas.df.to.plot$tissue != 'Not significant',]
  all.anovas.df.to.plot$tissue <- as.factor(all.anovas.df.to.plot$tissue)
  
  
  pdf(paste0(HOMEDIR, '/plots/',args[2],'/volcano_sign.pdf'), width = 10, height = 7)
  all.anovas.df.to.plot <- dataset %>%
    select(DRUG_ID, DRUG_NAME) %>%
    unique() %>%
    mutate(DRUG_ID=as.character(DRUG_ID)) %>%
    right_join(all.anovas.df.to.plot, by = c('DRUG_ID' = 'drug')) %>%
    arrange(corrected.p.value...q.value.) %>%
    mutate(labels=ifelse(corrected.p.value...q.value. < 0.02,
                         paste0(DRUG_NAME, '\n', alteration),
                         ''))
  # Upper limit of labels
  all.anovas.df.to.plot$labels[30:nrow(all.anovas.df.to.plot)] <- ''
  
  ggplot(all.anovas.df.to.plot,
         aes(x=effect, y=p.value,
             colour=tissue, size = n.samples)) +
    scale_y_continuous(trans=reverselog_trans(10),
                       label=scientific_10,
                       limits=c(0.05, min(all.anovas.df$p.value)/10)) +
    scale_x_continuous(breaks =seq(floor(min(all.anovas.df$effect)),
                                   ceiling(max(all.anovas.df$effect)+2), 2),
                       limits=c(floor(min(all.anovas.df$effect)),
                                ceiling(max(all.anovas.df$effect)+2))) +
    geom_point(aes(size=(n.samples+10)/10)) +
    geom_vline(xintercept=0, lwd=0.3) +
    geom_vline(xintercept=c(-1,1), lwd=0.3, linetype=2) +
    geom_hline(yintercept=max.pval, lwd=0.3, linetype=2) +
    annotate("text", x=min(all.anovas.df$effect) + 0.65, y=max.pval-0.0003, size=3.5,
             label=paste('P-value = ', max.pval, sep=""), colour="black") +
    annotate("text", x=min.effect.size + 1.5,
             y=min(all.anovas.df$p.value)/9, size=3.5,
             label=paste('Effect size of ', 1, sep=""), colour="black") +
    xlab("Drug effect") +
    ylab("P-value") +
    theme_bw() +
    theme(plot.title=element_text(size=14, face="bold", hjust = 0.5)) +
    scale_colour_manual('Tissue', values = tissue.cols,
                        breaks = names(tissue.cols)) +
    geom_text_repel(mapping=aes(x=effect,
                                y=p.value,
                                label=labels), colour="gray40", size=2,
                    box.padding = unit(0.2, 'lines'),
                    segment.alpha = 0.4) +
    scale_size('N. altered samples',
               breaks = c((40+10)/10, (20+10)/10, (3+10)/10),
               labels = c('40', '20', '3'))
  invisible(dev.off())
  
  pdf(paste0(HOMEDIR, '/plots/',args[2],'/barchart_sign.pdf'), width = 3.6, height = 2.4)
  barchart.volcano(all.anovas.df, fixed.y = FALSE)
  invisible(dev.off())
  
  
  ## All data
  
  pdf(paste0(HOMEDIR, '/plots/',args[2],'/volcano_all.pdf'), width = 10.5, height = 7)
  ggplot(data=all.anovas.unfiltered,
         aes(x=effect, y=p.value,
             colour=tissue, size = n.samples)) +
    scale_y_continuous(trans=reverselog_trans(10),
                       label=scientific_10,
                       limits=c(1, min(all.anovas.unfiltered$p.value)/10)) +
    scale_x_continuous(breaks =seq(floor(min(all.anovas.unfiltered$effect)),
                                   ceiling(max(all.anovas.unfiltered$effect)+2), 2),
                       limits=c(floor(min(all.anovas.unfiltered$effect)),
                                ceiling(max(all.anovas.unfiltered$effect)+2))) +
    geom_point(aes(size=(n.samples+10)/10)) +
    geom_vline(xintercept=0, lwd=0.3) +
    geom_vline(xintercept=c(-1,1), lwd=0.3, linetype=2) +
    geom_hline(yintercept=max.pval, lwd=0.3, linetype=2) +
    annotate("text", x=min(all.anovas.unfiltered$effect) + 0.65, y=max.pval-0.0003, size=3.5,
             label=paste('P-value = ', max.pval, sep=""), colour="black") +
    annotate("text", x=min.effect.size + 1.5,
             y=min(all.anovas.unfiltered$p.value)/9, size=3.5,
             label=paste('Effect size of ', 1, sep=""), colour="black") +
    xlab("Drug effect") +
    ylab("P-value") +
    theme_bw() +
    theme(plot.title=element_text(size=14, face="bold", hjust = 0.5)) +
    scale_colour_manual('Tissue', values = tissue.cols,
                        breaks = names(tissue.cols)) +
    scale_size('N. altered samples',
               breaks = c((40+10)/10, (20+10)/10, (3+10)/10),
               labels = c('40', '20', '3'))
  invisible(dev.off())
  
}
