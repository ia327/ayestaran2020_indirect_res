library(tidyverse)
library(ggpubr)
library(gridExtra)
library(ggbeeswarm)
library(ggrepel)
library(cowplot)
library(stringdist)

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


############### General stats ###############
library(VennDiagram)
gdsc.cosm.dr <- all.gdsc %>%
  filter(!is.na(LN_IC50)) %>%
  select(COSMIC_ID, DRUG_ID, DRUG_NAME)

cell.lines.per.cancer.type <- all.gdsc %>%
  filter(!is.na(LN_IC50)) %>%
  select(TCGA_DESC, COSMIC_ID) %>%
  distinct() %>%
  count(TCGA_DESC)

unique.gdsc.drugs <- gdsc.cosm.dr %>%
  select(DRUG_NAME) %>%
  distinct()
## 495 screened, but 397 unique drugs (there were replicates and so on)


gdsc.stats <- gdsc.cosm.dr %>%
  select(-DRUG_NAME) %>%
  gather(key, value) %>%
  distinct() %>%
  count(key) %>%
  mutate(DATASET='GDSC')

ctrp.cosm.dr <- all.ctrp %>%
  filter(!is.na(LN_IC50)) %>%
  select(COSMIC_ID, DRUG_ID, DRUG_NAME) %>%
  mutate(COSMIC_ID=map_chr(str_split(COSMIC_ID, '_'), 1)) %>%
  distinct()

unique.ctrp.drugs <- ctrp.cosm.dr %>%
  select(DRUG_NAME, DRUG_ID) %>%
  distinct()
## 545 screened, 545 unique drugs (no replicates)


ctrp.stats <- ctrp.cosm.dr %>%
  select(-DRUG_NAME) %>%
  gather(key, value) %>%
  distinct() %>%
  count(key) %>%
  mutate(DATASET='CTRP')

gdsc_drug_info <- read_csv(paste0(HOMEDIR,'data/annotated/screened_compounds_rel_8.0.csv'),
                           col_types = cols())
ctrp_drug_info <- read.table(paste0(HOMEDIR, 'data/annotated/CTRP/v20.meta.per_compound.txt'),
                             sep = '\t', header = TRUE,
                             quote = '') %>%
  select(c(master_cpd_id, cpd_name, top_test_conc_umol,
           gene_symbol_of_protein_target, target_or_activity_of_compound))

### Cross compare drug names in GDSC and CTRP
gdsc.ctrp.drug.lut <- map_dfr(1:nrow(ctrp_drug_info), function(i) {
  row <- ctrp_drug_info[i,]
  x <- ctrp_drug_info$cpd_name[i]

  ## Check DRUG_NAME field
  match <- which(grepl(paste0('^',tolower(x), '$'), tolower(gdsc_drug_info$DRUG_NAME)))

  if (length(match) != 0){
    matchtype <- 'exact'
  } else {
    match <- which(stringdist(tolower(x), tolower(gdsc_drug_info$DRUG_NAME)) <= 1)
    if (length(match) != 0){
      matchtype <- 'close'
    }
  }

  if (length(match) == 0) {
    ## Check SYNONYMS
    match <- which(grepl(paste0('(^|[[:space:]])',tolower(x), '($|,| )'), tolower(gdsc_drug_info$SYNONYMS)))
    if (length(match) != 0){
      matchtype <- 'exact'
    } else {
      all.synonyms <- str_split(tolower(gdsc_drug_info$SYNONYMS), ', ')
      match <- which(map_lgl(all.synonyms, ~ any(stringdist(tolower(x), .x) <= 1)))
      if (length(match) != 0){
        matchtype <- 'close'
      }
    }
  }

  if (length(match) == 0){
    return(NULL)
  }

  for (tt in 1:length(match)){
    if (tt == 1){
      res.to.bind <- row
    } else {
      res.to.bind <- rbind(res.to.bind, row)
    }
  }
  res <- gdsc_drug_info[match,]
  names(res) <- paste0('GDSC_', names(res))
  names(res.to.bind) <- paste0('CTRP_', names(res.to.bind))
  res <- bind_cols(res, res.to.bind, data.frame(matchtype=as.character(rep(matchtype, length(match)))))
  return(res)
})

## Collapse replicates and manually filter inaccurate close matches
unique.gdsc.ctrp.drug.lut <- gdsc.ctrp.drug.lut %>%
  nest(GDSC_DRUG_REPLICATE_LIST = c(GDSC_DRUG_ID, GDSC_SCREENING_SITE))%>%
  filter(!(matchtype=='close' & GDSC_DRUG_NAME %in% c('MIM1', 'CAY10566', 'SB52334',
                                                      'WZ4003', 'ML323'))) %>%
  select(-matchtype)


cross.stats <- data.frame(COSMIC_ID=length(intersect(unique(gdsc.cosm.dr$COSMIC_ID),
                                                     unique(map(str_split(ctrp.cosm.dr$COSMIC_ID, '_'), 1)))))

all.stats <- bind_rows(gdsc.stats, ctrp.stats) %>%
  spread(key, n) %>%
  rename(`Cell lines`=COSMIC_ID,
         Compounds=DRUG_ID) %>%
  arrange(desc(DATASET))

cairo_pdf(paste0(HOMEDIR, '/manuscript/Figures/figure1d.pdf'), height = 4.5, width = 9)
v1 <- draw.pairwise.venn(area1 = gdsc.stats$n[gdsc.stats$key == 'COSMIC_ID'],
                         area2 = ctrp.stats$n[ctrp.stats$key == 'COSMIC_ID'],
                         cross.area = cross.stats$COSMIC_ID,
                         category = c("GDSC", "CTRP"),
                         fill = c("#385bfd", "#FDDA38"),
                         alpha = rep(0.7, 2),
                         cex = 2,
                         # col = rep('#555555', 2),
                         col = rep(NA, 2),
                         cat.dist = c(-0.04, -0.04),
                         cat.cex = rep(1.5, 2),
                         ext.dist = -0.14,
                         ext.pos = 45,
                         ext.length = 0.7,
                         cat.fontfamily = rep('Arial', 2),
                         fontfamily = 'Arial',
                         margin = rep(0.035, 4),
                         ind = FALSE)
v1 <- add.title(v1, 'Cell lines',
                pos = c(0.5, 1),
                fontfamily = 'Arial',
                cex = 2.5)

v2 <- draw.pairwise.venn(area1 = nrow(unique.gdsc.drugs),
                         area2 = nrow(unique.ctrp.drugs),
                         cross.area = nrow(unique.gdsc.ctrp.drug.lut),
                         category = c("GDSC", "CTRP"),
                         fill = c("#385bfd", "#FDDA38"),
                         alpha = rep(0.7, 2),
                         cex = 2,
                         # col = rep('#555555', 2),
                         col = rep(NA, 2),
                         cat.dist = c(-0.06, -0.06),
                         cat.cex = rep(1.5, 2),
                         cat.pos = c(45, -45),
                         ext.dist = -0.14,
                         ext.pos = 45,
                         inverted = T,
                         ext.length = 0.7,
                         cat.fontfamily = rep('Arial', 2),
                         fontfamily = 'Arial',
                         margin = rep(0.035, 4),
                         ind = FALSE)
v2 <- add.title(v2, 'Screened compounds',
                pos = c(0.5, 1),
                fontfamily = 'Arial',
                cex = 2.5)

grid.arrange(grobTree(v1), grobTree(v2), nrow = 1, padding = 2)
dev.off()

###############  ###############  ###############


################ Figure 1.C curves from raw data GDSC1 ###############
load(paste0(HOMEDIR, 'data/drug/curve_fitting_stats.RData'))
all.bems.tidy <- all.bems.tidy.gdsc

stad.cosmic.savolitinib <- all.bems.tidy.gdsc %>%
  filter(Tissue=='STAD' & alteration == 'gain:cnaSTAD45 (MET)') %>%
  select(COSMIC_ID, mut.status) %>%
  distinct() %>%
  mutate(mut.status=as.character(mut.status)) %>%
  mutate(COSMIC_ID=as.numeric(COSMIC_ID))

curve_sens <- plotResponse_custom(model_stats = gdsc_model_stats_2,
                                  cell_line = stad.cosmic.savolitinib,
                                  drug_identifier = "1936_1") +
  ggtitle('STAD - Savolitinib (MET)',
          'Dose response') +
  theme(legend.position = c(0.8, 0.8),
        axis.title.x = element_blank(),
        plot.margin = unit(c(4.8, 9.8, 14.4, 9.8), 'pt')) +
  scale_color_manual('MET status', breaks=c(0,1),
                     labels=c('WT', 'CN gain'),
                     values = c('#00004488', '#FF680088')) +
  geom_segment(aes(x=log10(exp(max(IC50))), xend = log10(exp(max(IC50)))),
               y = 0, yend = -0.1, linetype = 2, alpha = 0.005) +
  geom_segment(aes(x=log10(exp(min(IC50))), xend = log10(exp(min(IC50)))),
               y = 0, yend = -0.1, linetype = 2, alpha = 0.005)


coread.cosmic.nutlin <- all.bems.tidy.gdsc %>%
  filter(Tissue=='COREAD' & alteration == 'TP53_mut') %>%
  select(COSMIC_ID, mut.status) %>%
  distinct() %>%
  mutate(mut.status=as.character(mut.status)) %>%
  mutate(COSMIC_ID=as.numeric(COSMIC_ID))

curve_resis <- plotResponse_custom(model_stats = gdsc_model_stats_1,
                                   cell_line = coread.cosmic.nutlin,
                                   drug_identifier = "1047_8") +
  ggtitle('COREAD - Nutlin-3a',
          'Dose response') +
  theme(legend.position = c(0.8, 0.8),
        axis.title.x = element_blank(),
        plot.margin = unit(c(4.8, 9.8, 14.4, 9.8), 'pt')) +
  scale_color_manual('TP53 status', breaks=c(0,1),
                     labels=c('WT', 'Mut'),
                     values = c('#00004488', '#FF680088')) +
  geom_segment(aes(x=log10(exp(max(IC50))), xend = log10(exp(max(IC50)))),
               y = 0, yend = -0.1, linetype = 2, alpha = 0.005) +
  geom_segment(aes(x=log10(exp(min(IC50))), xend = log10(exp(min(IC50)))),
               y = 0, yend = -0.1, linetype = 2, alpha = 0.005)


luad.cosmic.mut <- all.bems.tidy.gdsc %>%
  filter(Tissue=='LUAD' & alteration == 'EGFR_mut') %>%
  select(COSMIC_ID, mut.status) %>%
  distinct() %>%
  mutate(mut.status=as.character(mut.status)) %>%
  mutate(COSMIC_ID=as.numeric(COSMIC_ID))

curve_outlier <- plotResponse_custom(model_stats = gdsc_model_stats_1,
                                     cell_line = luad.cosmic.mut,
                                     drug_identifier = "1010_0.5") +
  ggtitle('LUAD - Gefitinib',
          'Dose response') +
  theme(legend.position = c(0.8, 0.8),
        axis.title.x = element_blank(),
        plot.margin = unit(c(4.8, 9.8, 14.4, 9.8), 'pt')) +
  scale_color_manual('EGFR status', breaks=c(0,1),
                     labels=c('WT', 'Mut'),
                     values = c('#00004488', '#FF680088')) +
  geom_segment(aes(x=log10(exp(max(IC50))), xend = log10(exp(max(IC50)))),
               y = 0, yend = -0.1, linetype = 2, alpha = 0.005) +
  geom_segment(aes(x=log10(exp(min(IC50))), xend = log10(exp(min(IC50)))),
               y = 0, yend = -0.1, linetype = 2, alpha = 0.005)



bp_dens_sens <- plot_bp_and_dens('STAD', '1936.2', 'gain:cnaSTAD45 (MET)')
bp_dens_res <- plot_bp_and_dens('COREAD', '1047.1', 'TP53_mut')
bp_dens_outlier <- plot_bp_and_dens('LUAD', '1010.1', 'EGFR_mut')


curves_and_bp_and_dens <- ggarrange(
  ggarrange(curve_sens, curve_resis, curve_outlier, ncol=3,
            labels = c('B', 'C', 'D')),
  ggarrange(bp_dens_sens, bp_dens_res, bp_dens_outlier, ncol = 3),
  nrow = 2,
  heights = c(2.5, 3))

ggsave(paste0(HOMEDIR, 'manuscript/Figures/figure1b_d.pdf'),
       curves_and_bp_and_dens, height = 5.5, width = 12)

################  ################  ################


################ Figure 1.X show extrapolation in GDSC vs high dose in CTRP ################
luad.cl <- all.gdsc %>%
  filter(TCGA_DESC == 'LUAD') %>%
  select(GDSC_cell_line, COSMIC_ID) %>%
  distinct() %>%
  mutate(COSMIC_ID=as.numeric(COSMIC_ID))

two.resistant.cl <- tibble(COSMIC_ID=c(724874, 905942),
                               mut.status = c('0', '1')) %>%
  left_join(luad.cl, by = 'COSMIC_ID')

two.resistant.cl.ctrp <- two.resistant.cl %>%
  left_join(all.ctrp, by = 'GDSC_cell_line') %>%
  mutate(COSMIC_ID=COSMIC_ID.x) %>%
  select(CL, mut.status, GDSC_cell_line) %>%
  distinct() %>%
  rename(COSMIC_ID = CL)

p1 <- plotResponse_custom(model_stats = gdsc_model_stats_1,
                          cell_line = two.resistant.cl,
                          drug_identifier = "1010_0.5") +
  scale_color_manual('Cell line', breaks=c(0,1),
                     labels = two.resistant.cl$GDSC_cell_line,
                     values = c('#0000CC77', '#FF000077')) +
  geom_point(aes(x = lx, y=1-y, colour = mut.status)) +
  geom_segment(aes(x=log10(exp(IC50)), y=0.5, xend=2.2, yend= 0.7), size= 0.1) +
  annotate('text', x=2, y=0.7, label = 'Extrapolated IC50',
           vjust = -0.5, hjust=0.4, size = 4) +
  scale_x_continuous(limits = c(-4, 4)) +
  ggtitle('GDSC1 - Gefitinib')

p2 <- plotResponse_custom(model_stats = ctrp_model_stats,
                          cell_line = two.resistant.cl.ctrp,
                          drug_identifier = c('52926_225',
                                              '52926_557')) +
  scale_color_manual('', breaks=c(0,1),
                     values = c('#0000CC77', '#FF000077')) +
  geom_point(aes(x = lx, y=1-y, colour = mut.status)) +
  scale_x_continuous(limits = c(-4, 4)) +
  ggtitle('CTRP - Gefitinib')


grid_extra_curves <- ggarrange(p1, p2, nrow = 2,
                               ncol = 1,
                               labels = c('j', 'k'),
                               common.legend = TRUE,
                               legend = 'right')
ggsave(paste0(HOMEDIR, 'manuscript/Figures/figure1c_extra.pdf'),
       grid_extra_curves,
       height = 5, width = 5)
################  ################  ################



################ Figure 1H Frequency of mutations pancancer ################
rename.cna.mut <- function(name.cna.or.mut) {
  if(grepl('_mut', name.cna.or.mut)){
    sp <- strsplit(name.cna.or.mut, '_')[[1]]
    return(data.frame(type=sp[2],
                      gene=sp[1]))
  } else {
    # Skip cases with no drivers (no genes in the cna name)
    if (!grepl('\\(', name.cna.or.mut)){
      return(data.frame(type=NA, gene=NA))
    }
    sp <- strsplit(name.cna.or.mut, ':')[[1]]
    # cna.genes <- strsplit(sub('\\)', '', sub('cna.* \\(', '', sp[2])), ',')[[1]]
    cna.genes <- sub('\\)', '', sub('cna.* \\(', '', sp[2]))
    return(data.frame(type=sp[1], gene=cna.genes))
  }
}
bar.palette <- c('#EA8223', '#FDDA38','#328582', '#9A8636', '#4F3528')
names(bar.palette) <- c('mut', 'gain', 'loss', 'mut+gain', 'mut+loss')



# n.cells.pancan <- all.bems.tidy %>%
#   select(COSMIC_ID) %>%
#   unique() %>%
#   nrow()
#
# mut.freqs.pancan <- all.bems.tidy %>%
#   filter(!grepl('HypMET', alteration)) %>%
#   mutate(alt.nest=map(alteration, rename.cna.mut)) %>%
#   unnest() %>%
#   na.omit() %>%
#   filter(mut.status == 1) %>%
#   select(COSMIC_ID, gene, type) %>%
#   group_by(COSMIC_ID, gene) %>%
#   summarise(type=paste0(type, collapse = '+')) %>%
#   mutate(type=ifelse(type=='gain+gain', 'gain', type)) %>%
#   ungroup() %>%
#   select(-COSMIC_ID) %>%
#   mutate(total=n.cells.pancan) %>%
#   add_count(gene, type, name='case') %>%
#   unique() %>%
#   group_by(gene, type) %>%
#   summarise(freq=case/total) %>%
#   arrange(-freq) %>%
#   ungroup()
#
# pdf(paste0(HOMEDIR, '/manuscript/Figures/figure1h.pdf'), height = 3.2, width = 5)
# mut.freqs.pancan %>%
#   group_by(gene) %>%
#   mutate(total.freq=sum(freq)) %>%
#   ungroup() %>%
#   arrange(total.freq) %>%
#   mutate(gene=fct_inorder(gene)) %>%
#   filter(total.freq >= 0.05) %>%
#   ggplot(aes(x=gene, y=freq, fill=type)) +
#   geom_bar(stat = 'identity') +
#   coord_flip() +
#   scale_y_continuous('Frequency accross GDSC cell lines',
#                      breaks =  seq(0, 0.8, by = 0.1),
#                      expand = c(0.01, 0)) +
#   theme_classic() +
#   theme(axis.title.y = element_blank(),
#         axis.line.x = element_blank(),
#         panel.grid.major.x = element_line(linetype = 3, colour = '#00000068'),
#         # plot.margin = unit(c(0, 0, 0, 0), 'lines'),
#         legend.position = c(0.8, 0.3)) +
#   scale_fill_manual('Mutation type',
#                     values = bar.palette,
#                     breaks = c('mut', 'gain', 'loss', 'mut+gain', 'mut+loss'),
#                     labels = c('SNV', 'CN gain', 'CN loss',
#                                'CN gain + SNV', 'CN loss + SNV'))
# dev.off()



bar1 <- plot_mut_freq('STAD', min_freq = 0.15)

bar2 <- plot_mut_freq('COREAD', min_freq = 0.15)

bar3 <- plot_mut_freq('LUAD', min_freq = 0.15)

bar_legend <- bar3 +
  theme(legend.direction = 'horizontal',
        legend.position = 'bottom')

bar_legend <-  get_legend(bar_legend)

bar1_print <- bar1 + theme(legend.position = 'none')
bar2_print <- bar2 + theme(legend.position = 'none')
bar3_print <- bar3 + theme(legend.position = 'none')

bar_all <- ggarrange(ggarrange(bar1_print,
                               bar2_print,
                               bar3_print, ncol = 3,
                               labels = c('A', 'B', 'C')),
                     bar_legend, nrow = 2,
                     heights = c(9, 1))


ggsave(paste0(HOMEDIR, '/manuscript/Figures/figure1g_i.pdf'),
       bar_all, height = 4, width = 12)


################  ################  ################



################ Nutlin 3a vs MDM2 exp ################
exp.data <- read_tsv(paste0(HOMEDIR,
                            'data/expression/Cell_line_RMA_proc_basalExp.txt.zip'))
exp.data <- exp.data %>%
  gather(COSMIC_ID, exp.value, -GENE_SYMBOLS, -GENE_title) %>%
  mutate(COSMIC_ID=gsub('DATA.', '', COSMIC_ID)) %>%
  select(-GENE_title)


mdm2_nutlin <- all.gdsc %>%
  filter(DRUG_NAME == 'Nutlin-3a (-)') %>%
  select(COSMIC_ID, LN_IC50, TCGA_DESC) %>%
  inner_join(filter(exp.data, GENE_SYMBOLS == 'MDM2'),
             by = 'COSMIC_ID') %>%
  inner_join(filter(all.bems.tidy, alteration == 'TP53_mut'),
             by = 'COSMIC_ID') %>%
  ggplot(aes(x=exp(LN_IC50),
             y=exp.value,
             col=as.factor(mut.status))) +
  geom_point() +
  # facet_wrap(~Tissue, scales = 'free') +
  scale_color_manual('TP53 status', breaks=c(0,1),
                     labels=c('WT', 'Mut'),
                     values = c('#00004488', '#FF680088')) +
  scale_x_log10(expression('Nutlin-3a' ~ IC[50] ~ (mu * M))) +
  scale_y_continuous('MDM2 expression') +
  theme_pubr() +
  theme(legend.position = c(0.85, 0.85))

ggsave(paste0(HOMEDIR, 'manuscript/Figures/mdm2_nutlin_pancan.pdf'),
       mdm2_nutlin, height = 5, width = 5)

################  ################  ################
