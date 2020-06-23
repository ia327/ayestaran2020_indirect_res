library(tidyverse)
library(ggpubr)
library(gridExtra)
library(ggbeeswarm)
library(ggrepel)
library(scales)
library(readxl)
library(RColorBrewer)
library(cowplot)

HOMEDIR <- '../../'
DATADIR <- paste0(HOMEDIR, 'data/combined_for_pipeline/')
load(paste0(DATADIR, 'gdsc_dr_data.RData'))
load(paste0(DATADIR, 'ctrp_dr_data.RData'))
load(paste0(DATADIR, 'utils.RData'))
load(paste0(DATADIR, 'bems_tidy.RData'))
load(paste0(DATADIR, 'wes_annotation.RData'))
load(paste0(DATADIR, 'cn_annotation.RData'))
load(paste0(DATADIR, 'cancer_gene_lists.RData'))
source(paste0(HOMEDIR, 'manuscript/Figures/figures_utils.R'))


################ GDSC Volcano plot from ANOVA results ################
OUTDIR <- paste0(HOMEDIR, '/results/gdsc/ANOVA_out/')
used.tissues <- all.gdsc %>%
  select(TCGA_DESC, COSMIC_ID) %>%
  unique() %>%
  count(TCGA_DESC) %>%
  filter(n >= 15 & TCGA_DESC %in% all.tissues) %>%
  select(TCGA_DESC) %>% unlist() %>% unname() %>% as.character()

all.anovas.df.to.plot <- map_dfr(paste0(OUTDIR, used.tissues, '_anova_out.csv'),
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
# max.fdr <- 0.20
min.effect.size <- -1 # Min in absolute terms

all.anovas.df.to.plot$tissue[all.anovas.df.to.plot$p.value > max.pval |
                               # all.anovas.df.to.plot$corrected.p.value...q.value. > max.fdr |
                               abs(all.anovas.df.to.plot$effect) < abs(min.effect.size)] <- 'Not significant'
all.anovas.df.to.plot <- all.anovas.df.to.plot[all.anovas.df.to.plot$tissue != 'Not significant',]

all.anovas.df.to.plot <- all.gdsc %>%
  select(DRUG_ID, DRUG_NAME) %>%
  unique() %>%
  mutate(DRUG_ID=as.character(DRUG_ID)) %>%
  right_join(all.anovas.df.to.plot, by = c('DRUG_ID' = 'drug')) %>%
  arrange(p.value) %>%
  mutate(alteration.nice=sapply(alteration, function(x) ifelse(grepl('_mut', x),
                                                               gsub('_', '', x),
                                                               rename.cna(x)))) %>%
  mutate(labels=ifelse(corrected.p.value...q.value. < 0.01,
                       paste0(DRUG_NAME, '\n', alteration.nice),
                       ''),
         labels=ifelse(effect > 0, paste0(labels, '\n', tissue), labels))
# Upper limit of labels
all.anovas.df.to.plot$labels[9:nrow(all.anovas.df.to.plot)] <- ''
all.anovas.df.to.plot$tissue[all.anovas.df.to.plot$effect > 0] <- 'Not significant'



ggvolcano <- ggplot(all.anovas.df.to.plot,
                    aes(x=effect, y=p.value,
                        colour=tissue, size = n.samples)) +
  scale_y_continuous(trans=reverselog_trans(10),
                     label=scientific_10)+
                     # limits=c(0.05, min(all.anovas.df.to.plot$p.value)/10)) +
  scale_x_continuous(breaks =seq(floor(min(all.anovas.df.to.plot$effect)),
                                 ceiling(max(all.anovas.df.to.plot$effect)+2), 2),
                     limits=c(floor(min(all.anovas.df.to.plot$effect)),
                              ceiling(max(all.anovas.df.to.plot$effect)+2))) +
  # expand_limits(x=c(-8, 6)) +
  geom_point(aes(size=(n.samples+10)/10), alpha=0.7) +
  geom_vline(xintercept=0, lwd=0.75) +
  geom_vline(xintercept=c(-1,1), lwd=0.3, linetype=2, alpha=0.6) +
  geom_hline(yintercept=max.pval, lwd=0.3, linetype=2, alpha=0.6) +
  annotate("text", x=min(all.anovas.df.to.plot$effect) + 0.65, y=max.pval, vjust=1.3, size=3.5,
           label=paste('P-value = ', max.pval, sep=""), colour="black") +
  # annotate("text", x=min.effect.size + 1.5,
  #          y=min(all.anovas.df.to.plot$p.value)/9, size=3.5,
  #          label=paste('Effect size of ', 1, sep=""), colour="black") +
  xlab("Drug effect") +
  ylab("P-value") +
  theme_classic2() +
  theme(plot.title=element_text(size=14, face="bold", hjust = 0.5),
        axis.line.y = element_blank(),
        panel.grid.major.y = element_line(linetype = 3, colour = '#00000055')) +
  scale_colour_manual('Tissue', values = tissue.cols,
                      breaks = used.tissues) +
  geom_text_repel(mapping=aes(x=effect,
                              y=p.value,
                              label=labels), colour="gray40", size=3,
                  box.padding = unit(0.2, 'lines'),
                  point.padding = 0.2,
                  segment.alpha = 0.4) +
  scale_size('N. altered samples',
             breaks = c((40+10)/10, (20+10)/10, (3+10)/10),
             labels = c('40', '20', '3')) +
  ggtitle('GDSC sensitivity markers')

ggvolcano <- ggarrange(ggvolcano, labels = 'A')

ggsave(paste0(HOMEDIR, 'manuscript/Figures/figure1v.pdf'),
       ggvolcano, height = 5.5, width = 8)

################  ################  ################

################ CTRP Volcano plot from ANOVA results ################
OUTDIR <- paste0(HOMEDIR, '/results/ctrp/ANOVA_out/')
used.tissues <- all.ctrp %>%
  select(TCGA_DESC, COSMIC_ID) %>%
  unique() %>%
  count(TCGA_DESC) %>%
  filter(n >= 15 & TCGA_DESC %in% all.tissues) %>%
  select(TCGA_DESC) %>% unlist() %>% unname() %>% as.character()

all.anovas.df.to.plot <- map_dfr(paste0(OUTDIR, used.tissues, '_anova_out.csv'),
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
# max.fdr <- 0.20
min.effect.size <- -1 # Min in absolute terms

all.anovas.df.to.plot$tissue[all.anovas.df.to.plot$p.value > max.pval |
                               # all.anovas.df.to.plot$corrected.p.value...q.value. > max.fdr |
                               abs(all.anovas.df.to.plot$effect) < abs(min.effect.size)] <- 'Not significant'
all.anovas.df.to.plot <- all.anovas.df.to.plot[all.anovas.df.to.plot$tissue != 'Not significant',]

all.anovas.df.to.plot <- all.ctrp %>%
  select(DRUG_ID, DRUG_NAME) %>%
  unique() %>%
  mutate(DRUG_ID=as.character(DRUG_ID)) %>%
  right_join(all.anovas.df.to.plot, by = c('DRUG_ID' = 'drug')) %>%
  arrange(p.value) %>%
  mutate(alteration.nice=sapply(alteration, function(x) ifelse(grepl('_mut', x),
                                                               gsub('_', '', x),
                                                               rename.cna(x)))) %>%
  mutate(labels=ifelse(corrected.p.value...q.value. < 0.01,
                       paste0(DRUG_NAME, '\n', alteration.nice),
                       ''),
         labels=ifelse(effect > 0, paste0(labels, '\n', tissue), labels))
# Upper limit of labels
all.anovas.df.to.plot$labels[9:nrow(all.anovas.df.to.plot)] <- ''
all.anovas.df.to.plot$tissue[all.anovas.df.to.plot$effect > 0] <- 'Not significant'



ggvolcano <- ggplot(all.anovas.df.to.plot,
                    aes(x=effect, y=p.value,
                        colour=tissue, size = n.samples)) +
  scale_y_continuous(trans=reverselog_trans(10),
                     label=scientific_10) +
                     # limits=c(0.05, min(all.anovas.df.to.plot$p.value)/10)) +
  scale_x_continuous(breaks =seq(floor(min(all.anovas.df.to.plot$effect)),
                                 ceiling(max(all.anovas.df.to.plot$effect)+2), 2),
                     limits=c(floor(min(all.anovas.df.to.plot$effect)),
                              ceiling(max(all.anovas.df.to.plot$effect)+2))) +
  # expand_limits(x=c(-8, 6)) +
  geom_point(aes(size=(n.samples+10)/10), alpha=0.7) +
  geom_vline(xintercept=0, lwd=0.75) +
  geom_vline(xintercept=c(-1,1), lwd=0.3, linetype=2, alpha=0.6) +
  geom_hline(yintercept=max.pval, lwd=0.3, linetype=2, alpha=0.6) +
  annotate("text", x=min(all.anovas.df.to.plot$effect), y=max.pval, vjust=1.3, size=3.5,
           label=paste('P-value = ', max.pval, sep=""), colour="black") +
  # annotate("text", x=min.effect.size + 1.5,
  #          y=min(all.anovas.df.to.plot$p.value)/9, size=3.5,
  #          label=paste('Effect size of ', 1, sep=""), colour="black") +
  xlab("Drug effect") +
  ylab("P-value") +
  theme_classic2() +
  theme(plot.title=element_text(size=14, face="bold", hjust = 0.5),
        axis.line.y = element_blank(),
        panel.grid.major.y = element_line(linetype = 3, colour = '#00000055')) +
  scale_colour_manual('Tissue', values = tissue.cols,
                      breaks = used.tissues) +
  geom_text_repel(mapping=aes(x=effect,
                              y=p.value,
                              label=labels), colour="gray40", size=3,
                  box.padding = unit(0.2, 'lines'),
                  point.padding = 0.2,
                  segment.alpha = 0.4) +
  scale_size('N. altered samples',
             breaks = c((40+10)/10, (20+10)/10, (3+10)/10),
             labels = c('40', '20', '3')) +
  ggtitle('CTRP sensitivity markers')

ggvolcano <- ggarrange(ggvolcano, labels = 'B')

ggsave(paste0(HOMEDIR, 'manuscript/Figures/figureS2a.pdf'),
       ggvolcano, height = 5.5, width = 8)

################  ################  ################


################ Overview of results of pipeline ################
load(paste0(HOMEDIR, 'results/gdsc/putative_res_markers.RData'))
gdsc.assoc.filt <- associations %>% arrange(outlier.q.value, n.outliers)
gdsc.assoc.filt$clean.name <- unname(apply(gdsc.assoc.filt, 1,
                                           function(yy) {
                                             drug.name <- all.gdsc %>%
                                               filter(DRUG_ID == yy['drug']) %>%
                                               dplyr::select(DRUG_NAME) %>%
                                               unique() %>% unlist() %>% unname()
                                             paste0(yy['tissue'], ' - ',
                                                    drug.name, ' - ', yy['alteration'],
                                                    ' : ', yy['n.outliers'], ' out')
                                           }))
gdsc.assoc.filt$clean.name <- paste0(1:nrow(gdsc.assoc.filt), '. ',
                                     gdsc.assoc.filt$clean.name)

load(paste0(HOMEDIR, 'results/ctrp/putative_res_markers.RData'))
ctrp.assoc.filt <- associations %>% arrange(outlier.q.value, n.outliers)
ctrp.assoc.filt$clean.name <- unname(apply(ctrp.assoc.filt, 1,
                                           function(yy) {
                                             drug.name <- all.ctrp %>%
                                               filter(DRUG_ID == yy['drug']) %>%
                                               dplyr::select(DRUG_NAME) %>%
                                               unique() %>% unlist() %>% unname()
                                             paste0(yy['tissue'], ' - ',
                                                    drug.name, ' - ', yy['alteration'],
                                                    ' : ', yy['n.outliers'], ' out')
                                           }))
ctrp.assoc.filt$clean.name <- paste0(1:nrow(ctrp.assoc.filt), '. ',
                                     ctrp.assoc.filt$clean.name)


plot.putative.volcano <- function(dat, putative.out,
                                  label.qval = 0.05,
                                  label.effect.size = 0.3) {
  putative.out <- putative.out %>%
    mutate(tissue=ifelse(outlier.q.value > 0.15, 'Not significant', tissue),
           tissue=as.factor(tissue)) %>%
    left_join(unique(dplyr::select(dat, DRUG_ID, DRUG_NAME)),
              by = c('drug' = 'DRUG_ID')) %>%
    rownames_to_column(var = 'key.id')

  p <- putative.out %>%
    filter(effect.size >= 0) %>%
    mutate(alteration.nice=sapply(alteration, function(x) ifelse(grepl('_mut', x),
                                                                 gsub('_', '', x),
                                                                 rename.cna(x)))) %>%
    mutate(lbls=paste0(DRUG_NAME, ' resist. in\n', alteration.nice)) %>%
    mutate(lbls=ifelse(outlier.q.value < label.qval &
                         effect.size > label.effect.size,
                       lbls, NA)) %>%
    ggplot(aes(x=effect.size,
               y=outlier.q.value,
               col=tissue,
               size=n.outliers)) +
    geom_line(aes(group=paste0(tissue,association)),
              data=filter(putative.out, tissue != 'Not significant'),
              size=0.5, alpha=0.6, linetype = 2) +
    geom_point(alpha=0.8) +
    theme_classic2() +
    theme(plot.title=element_text(size=14, face="bold", hjust = 0.5),
          axis.line.y = element_blank(),
          panel.grid.major.y = element_line(linetype = 3, colour = '#00000055')) +
    scale_colour_manual('Tissue', values = tissue.cols,
                        breaks = names(tissue.cols)) +
    geom_vline(xintercept=0, lwd=0.75) +
    scale_size_continuous('N of outliers', range=c(2,6)) +
    scale_x_continuous('Normalised SD decrease') +
    scale_y_continuous('Adjusted p-val',
                       trans=reverselog_trans(10),
                       label=scientific_10) +
    geom_text_repel(aes(label=lbls), colour="#444444", # fill = '#FFFFFF44',
                    size=3, lineheight = 0.75,
                    point.padding = 0.1, box.padding = 0.4,
                    min.segment.length = 0,
                    segment.alpha = 0.8,force = 2)
  return(p)
}


p1 <- plot.putative.volcano(dat = all.gdsc,
                            putative.out = gdsc.assoc.filt,
                            label.effect.size = 0.5,
                            label.qval = 0.15) +
  ggtitle('GDSC') +
  theme(legend.position = 'none')

p2 <- plot.putative.volcano(dat = all.ctrp,
                            putative.out = ctrp.assoc.filt,
                            label.effect.size = 0.5,
                            label.qval = 0.08) +
  ggtitle('CTRP') +
  theme(legend.position = 'none')

# joint one just for the legend
pleg <- plot.putative.volcano(dat = all.gdsc,
                             putative.out = bind_rows(gdsc.assoc.filt,
                                                      ctrp.assoc.filt),
                             label.effect.size = 0.4) +
  theme(legend.position = 'right', legend.box = "vertical") +
  guides(color = guide_legend(order=1),
          size = guide_legend(order=2))
both_res <- ggarrange(ggarrange(p1, p2, ncol = 1,
                      labels = c('B', 'C')),
                      get_legend(pleg), nrow = 1, widths = c(9, 2))
ggsave(paste0(HOMEDIR, 'manuscript/Figures/figure2a_vert.pdf'),
       both_res, height = 10, width = 8)

################  ################  ################




################  Plot individual boxplots with outliers highlighted  ########

system(paste0('mkdir -p ', HOMEDIR, 'manuscript/Figures/outlier_indiv_gdsc'))
system(paste0('mkdir -p ', HOMEDIR, 'manuscript/Figures/outlier_indiv_ctrp'))


for (p in seq(nrow(gdsc.assoc.filt))){
  if (gdsc.assoc.filt$outlier.q.value[p] > 0.15){
    next
  }
  pdf(paste0(HOMEDIR, 'manuscript/Figures/outlier_indiv_gdsc/',
             gdsc.assoc.filt$tissue[p], '_',
             gdsc.assoc.filt$drug[p], '_',
             make.names(gdsc.assoc.filt$alteration[p]), '_',
             gdsc.assoc.filt$n.outliers[p],
             '.pdf'), height = 4, width = 2.5)
  g <- plot.comp.boxplot.nice(tissue = gdsc.assoc.filt$tissue[p],
                              drug.id = gdsc.assoc.filt$drug[p],
                              gene = gdsc.assoc.filt$alteration[p],
                              ref.data = all.gdsc,
                              n.outliers = gdsc.assoc.filt$n.outliers[p],
                              highlight = str_split(gdsc.assoc.filt$out.cell.lines[p], ', ')[[1]]) +
    scale_y_continuous(expression(log[10](IC[50]))) +
    annotate('text', x = 1.2, y= Inf,
             vjust = 1,
             label = paste0('q = ',
                            scientific(gdsc.assoc.filt$outlier.q.value[p])))
  print(g)
  dev.off()

}

all.bems.tidy <- all.bems.tidy.ctrp
for (p in seq(nrow(ctrp.assoc.filt))){
  if (ctrp.assoc.filt$outlier.q.value[p] > 0.15){
    next
  }
  pdf(paste0(HOMEDIR, 'manuscript/Figures/outlier_indiv_ctrp/',
             ctrp.assoc.filt$tissue[p], '_',
             ctrp.assoc.filt$drug[p], '_',
             make.names(ctrp.assoc.filt$alteration[p]), '_',
             ctrp.assoc.filt$n.outliers[p],
             '.pdf'), height = 4, width = 2.5)
  g <- plot.comp.boxplot.nice(tissue = ctrp.assoc.filt$tissue[p],
                              drug.id = ctrp.assoc.filt$drug[p],
                              gene = ctrp.assoc.filt$alteration[p],
                              ref.data = all.ctrp,
                              n.outliers = ctrp.assoc.filt$n.outliers[p],
                              highlight = str_split(ctrp.assoc.filt$out.cell.lines[p], ', ')[[1]]) +
    scale_y_continuous(expression(log[10](AUC))) +
    annotate('text', x = 1.2, y= Inf,
             vjust = 1,
             label = paste0('q = ',
                            scientific(ctrp.assoc.filt$outlier.q.value[p])))

  print(g)
  dev.off()
}



system(paste0('mkdir -p ', HOMEDIR, 'manuscript/Figures/luad_cases_label'))
all.bems.tidy <- all.bems.tidy.gdsc
for (p in seq(nrow(gdsc.assoc.filt))){
  if (gdsc.assoc.filt$outlier.q.value[p] > 0.15 |
      gdsc.assoc.filt$tissue[p] != 'LUAD'){
    next
  }
  label.data <- all.gdsc %>%
    filter(COSMIC_ID %in% c('924244','687800')) %>%
    filter(DRUG_ID== gdsc.assoc.filt$drug[p]) %>%
    mutate(mut.status=1,
           log10.dr=log10(exp(LN_IC50)),
           is.highlight=ifelse(COSMIC_ID == '687800', '1', '2'))

  pdf(paste0(HOMEDIR, 'manuscript/Figures/luad_cases_label/',
             gdsc.assoc.filt$tissue[p], '_',
             gdsc.assoc.filt$drug[p], '_',
             make.names(gdsc.assoc.filt$alteration[p]), '_',
             gdsc.assoc.filt$n.outliers[p],
             '.pdf'), height = 4, width = 2.5)
  g <- plot.comp.boxplot.nice(tissue = gdsc.assoc.filt$tissue[p],
                              drug.id = gdsc.assoc.filt$drug[p],
                              gene = gdsc.assoc.filt$alteration[p],
                              ref.data = all.gdsc,
                              n.outliers = gdsc.assoc.filt$n.outliers[p],
                              highlight = str_split(gdsc.assoc.filt$out.cell.lines[p], ', ')[[1]]) +
    scale_y_continuous(expression(log[10](IC[50]))) +
    geom_text_repel(data = label.data, aes(label=GDSC_cell_line),
                    point.padding=0.2,
                    direction = 'y',
                    nudge_x=0.6,
                    segment.size = 0.2,
                    hjust=0,
                    min.segment.length = unit(0, 'lines'),
                    size = 3) +
    scale_colour_manual(breaks = c('0', '1', '2'),
                        values = c(unname(tissue.cols['LUAD']),
                                   '#004d7f', '#FF0000')) +
    annotate('text', x = 1.2, y= Inf,
             vjust = 1,
             label = paste0('q = ',
                            scientific(gdsc.assoc.filt$outlier.q.value[p])))
  print(g)
  dev.off()
}

all.bems.tidy <- all.bems.tidy.ctrp
for (p in seq(nrow(gdsc.assoc.filt))){
  if (ctrp.assoc.filt$outlier.q.value[p] > 0.15 |
      ctrp.assoc.filt$tissue[p] != 'LUAD'){
    next
  }
  label.data <- all.ctrp %>%
    filter(GDSC_COSMIC %in% c('924244','687800')) %>%
    filter(DRUG_ID== ctrp.assoc.filt$drug[p]) %>%
    mutate(mut.status=1,
           log10.dr=log10(exp(LN_IC50)),
           is.highlight=ifelse(COSMIC_ID == '687800', '1', '2'))

  pdf(paste0(HOMEDIR, 'manuscript/Figures/luad_cases_label/',
             ctrp.assoc.filt$tissue[p], '_',
             ctrp.assoc.filt$drug[p], '_',
             make.names(ctrp.assoc.filt$alteration[p]), '_',
             ctrp.assoc.filt$n.outliers[p],
             '.pdf'), height = 4, width = 2.5)
  g <- plot.comp.boxplot.nice(tissue = ctrp.assoc.filt$tissue[p],
                              drug.id = ctrp.assoc.filt$drug[p],
                              gene = ctrp.assoc.filt$alteration[p],
                              ref.data = all.ctrp,
                              n.outliers = ctrp.assoc.filt$n.outliers[p],
                              highlight = str_split(ctrp.assoc.filt$out.cell.lines[p], ', ')[[1]]) +
    scale_y_continuous(expression(log[10](IC[50]))) +
    geom_text_repel(data = label.data, aes(label=GDSC_cell_line),
                    point.padding=0.2,
                    direction = 'y',
                    nudge_x=0.6,
                    segment.size = 0.2,
                    hjust=0,
                    min.segment.length = unit(0, 'lines'),
                    size = 3) +
    scale_colour_manual(breaks = c('0', '1', '2'),
                        values = c(unname(tissue.cols['LUAD']),
                                   '#004d7f', '#FF0000')) +
    annotate('text', x = 1.2, y= Inf,
             vjust = 1,
             label = paste0('q = ',
                            scientific(ctrp.assoc.filt$outlier.q.value[p])))
  print(g)
  dev.off()
}


#### CRISPR stuff
cr.data <- read_excel(paste0(HOMEDIR, 'extra_analyses/supp_gad.291948.116_Supplemental_Table_S1.xlsx'),
                      skip = 2, na = 'NA')

## Check for all mutations and CN altered genes unique to H1650 (vs the other EGFR mut cell lines)
unique.h1650.gene.list <- read_tsv(paste0(HOMEDIR, 'extra_analyses/h1650_genelist.txt'),
                                   col_names = FALSE)

liao.et.al.crispr <- cr.data %>%
  mutate(is.in.list=ifelse(Gene %in% c(unique.h1650.gene.list$X2, 'PTEN'), 'Y', 'N')) %>% # Add PTEN even though it's missing from annotation
  mutate(lbls=ifelse(-log10(p.value.CRISPR) > 4 & log2(FC.CRISPR) > 1,
                     Gene, '')) %>%
  # filter(is.in.list == 'Y') %>%
  mutate(lbls=ifelse(is.in.list=='Y', Gene, lbls)) %>%
  ggplot(aes(x=log2(FC.CRISPR),
             y=-log10(p.value.CRISPR),
             col=is.in.list)) +
  geom_point(alpha=0.6) +
  geom_text_repel(aes(label=lbls), show.legend = FALSE,
                  segment.alpha = 0.3, box.padding = 0.3) +
  ggtitle('CRISPR enrichment data - Gefitinib vs DMSO',
          'Liao et al - PC-9 cell line') +
  geom_vline(xintercept = 0) +
  theme_classic2() +
  scale_color_manual('Mutated in\nNCI-H1650',
                     breaks = c('Y', 'N'),
                     labels = c('Yes', 'No'),
                     values=c('#FF0000', '#000000')) +
  scale_y_continuous(expression(-log[10](Pval)),
                     limits=c(0, max(-log10(cr.data$p.value.CRISPR))+1)) +
  scale_x_continuous(expression(log[2](FC))) +
  theme(axis.line.y = element_blank(),
        legend.position = c(0.9, 0.2),
        legend.background = element_blank(),
        panel.grid.major.y = element_line(linetype = 2,
                                          colour = '#00000033'))

ggsave(paste0(HOMEDIR, 'manuscript/Figures/figure2c.pdf'),
       liao.et.al.crispr, height = 4.25, width = 6)


################# ################# #################
