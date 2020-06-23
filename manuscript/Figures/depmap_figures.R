library(tidyverse)
library(ggrepel)
library(ggpubr)
library(RColorBrewer)
library(ComplexHeatmap)
library(GGally)
library(circlize)
library(gridExtra)

HOMEDIR <- '../../'

sanger.crispr.dat <- read_csv(paste0(HOMEDIR, 'extra_analyses/sanger_crispr/gene_effect.csv'))
model.lut <- read_csv(paste0(HOMEDIR, 'extra_analyses/sanger_crispr/model_list_latest.csv.gz'))
DATADIR <- paste0(HOMEDIR, 'data/combined_for_pipeline/')
load(paste0(DATADIR, 'gdsc_dr_data.RData'))
load(paste0(DATADIR, 'ctrp_dr_data.RData'))
load(paste0(DATADIR, 'utils.RData'))
load(paste0(DATADIR, 'bems_tidy.RData'))
load(paste0(DATADIR, 'wes_annotation.RData'))
load(paste0(DATADIR, 'cn_annotation.RData'))
load(paste0(DATADIR, 'cancer_gene_lists.RData'))

########## Functions ############


get_delta_ess <- function(outlier, tissue, alt, exclude = NA, check_min_diff = 2) {

  all.cosmic.id <- all.bems.tidy.gdsc %>%
    filter(Tissue == tissue) %>%
    filter(alteration == alt)

  all.cosmic.id.sp <- split(all.cosmic.id$COSMIC_ID, all.cosmic.id$mut.status)

  wt <- model.lut %>%
    filter(COSMIC_ID %in% all.cosmic.id.sp$`0`) %>%
    select(BROAD_ID) %>%
    na.omit() %>%
    unlist() %>% unname()

  mut <- model.lut %>%
    filter(COSMIC_ID %in% all.cosmic.id.sp$`1`) %>%
    select(BROAD_ID) %>%
    na.omit() %>%
    unlist() %>% unname()

  if (!is.na(exclude)) {
    wt <- wt[wt != exclude]
    mut <- mut[mut != exclude]
  }

  crispr.out <- sanger.crispr.dat %>%
    filter(X1 %in% outlier) %>%
    as.data.frame() %>%
    column_to_rownames('X1') %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column('Gene')
  crispr.mut <- sanger.crispr.dat %>%
    filter(X1 %in% mut) %>%
    filter(X1 != outlier) %>%
    as.data.frame() %>%
    column_to_rownames('X1') %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column('Gene')
  crispr.wt <- sanger.crispr.dat %>%
    filter(X1 %in% wt) %>%
    as.data.frame() %>%
    column_to_rownames('X1') %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column('Gene')

  mut.diff <- crispr.out[,outlier] - crispr.mut[,-1]
  mut.diff.z <- scale(mut.diff)
  mut.min.diff <- rowSums(abs(mut.diff.z) > check_min_diff)/ncol(mut.diff.z)
  mean.mut.diff.z <- rowMeans(mut.diff.z)
  delta_ess <- data.frame(Gene=crispr.mut[,1], delta_ess_mut=mean.mut.diff.z,
                          mut_comps_over_min=mut.min.diff) %>%
    bind_cols(as.data.frame(mut.diff.z)) %>%
    mutate(Gene=map_chr(str_split(Gene, ' \\('), 1))


  wt.diff <- crispr.out[,outlier] - crispr.wt[,-1]
  wt.diff.z <- scale(wt.diff)
  wt.min.diff <-  rowSums(abs(wt.diff.z) > check_min_diff)/ncol(wt.diff.z)
  mean.wt.diff.z <- rowMeans(wt.diff.z)
  delta_ess$delta_ess_wt <- mean.wt.diff.z
  delta_ess$wt_comps_over_min <- wt.min.diff

  delta_ess <- delta_ess %>%
    bind_cols(as.data.frame(wt.diff.z))

  return(delta_ess)
}


plot_diff_essentiality_heatmap <- function(diffs, resis, sens) {
  diffs <- arrange(diffs, delta_ess_mut)
  mycol <- colorRamp2(c(-2, 0, 2), hcl.colors(3, 'Blue-Red'))
  cl.labels <- model.lut$model_name[sapply(c(resis, sens),
                                           function(x) which(model.lut$BROAD_ID == x))]

  all.mean.diffs <- diffs$delta_ess_mut

  sidecols <- colorRamp2(c(min(all.mean.diffs), 0, max(all.mean.diffs)),
                         c('#E57E00','#FFFFFF', '#7E00E5'))

  topcols <- c(rep('#888888', length(resis)),
               rep('#000000', length(sens)))

  topcols <- c(rep('Resistant', length(resis)),
               rep('Sensitive', length(sens)))

  directioncol <- diffs$delta_ess_mut > 0

  to.hm <- diffs %>%
    mutate(Gene=map_chr(str_split(Gene, ' \\('), 1)) %>%
    select(Gene, resis, sens) %>%
    as.data.frame() %>%
    column_to_rownames('Gene') %>%
    as.matrix() %>%
    `colnames<-`(cl.labels)

  hh <- Heatmap(to.hm, name = "Essentiality",
                col = mycol,
                column_names_rot = 30,
                row_split = directioncol,
                row_gap = unit(3, 'mm'),
                column_split = topcols,
                column_gap = unit(0, 'mm'),
                column_title = NULL,
                row_title = c('More essential', 'Less essential'),
                top_annotation = HeatmapAnnotation(Response = topcols,
                                                   col = list(Response = c('Resistant' = '#000000',
                                                                           'Sensitive' = '#999999')),
                                                   show_annotation_name = FALSE),
                left_annotation = rowAnnotation(`∆ess` = diffs$delta_ess_mut,
                                                `wt-like` = diffs$is.wt.like,
                                                col = list(`∆ess`=sidecols,
                                                           `wt-like`=c('FALSE'='#006100', 'TRUE'='#b5b5b5')),
                                                show_annotation_name = FALSE),
                show_parent_dend_line = FALSE,
                # cluster_rows = FALSE,
                border = TRUE)
  return(draw(hh, merge_legend = TRUE, newpage = FALSE))
}

plot_essentialities <- function(cl1, cl2,
                                label_thres = 5,
                                label_size = 3,
                                dat) {
  p <- dat %>%
    mutate(Gene.nice=map_chr(str_split(Gene, ' \\('), 1)) %>%
    mutate(diff.value=get(cl1) - get(cl2)) %>%
    mutate(diff.value.z=(diff.value-mean(diff.value))/sd(diff.value)) %>%
    mutate(lbls=ifelse(abs(diff.value.z) > label_thres, Gene.nice, '')) %>%
    mutate(is.labelled = lbls != '') %>%
    ggplot(aes(x=get(cl1),
               y=get(cl2),
               col=is.labelled)) +
    geom_point() +
    geom_text_repel(aes(label = lbls),
                    size = label_size,
                    segment.alpha = 0.3) +
    theme_classic2() +
    theme(legend.position = 'none',
          panel.grid.major = element_line(colour = '#88888844', size = 0.3),
          axis.line = element_blank()) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    scale_color_manual(values=c('#00000033', '#000000')) +
    geom_abline(slope = 1, intercept = 0, linetype=2, colour='#00000055')
  return(p)
}

plot_essentiality_diffs <- function(diff_data, diff_val1, diff_val2,
                                    diff_thres = 3) {

  g <-  diff_data %>%
    mutate(lbls=ifelse(abs(get(diff_val1)) > diff_thres &
                         abs(get(diff_val2)) > diff_thres,
                       Gene, '')) %>%
    ggplot(aes(x=get(diff_val1), y=get(diff_val2), col=(lbls==''))) +
    scale_colour_manual('',
                        values = c('#000000', '#00000022'),
                        guide = FALSE) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    geom_abline(slope = 1, intercept = 0,
                linetype = 2, colour = '#00000077') +
    geom_point() +
    geom_text_repel(aes(label=lbls),
                    col='#000000',
                    size=3,
                    point.padding = 0.3,
                    box.padding = 0,
                    min.segment.length = 0,
                    segment.alpha = 0.4,
                    alpha = 0.8,
                    show.legend = FALSE) +
    theme_classic2() +
    theme(axis.line = element_blank())
}



########################## Figure 3i ##########################
# DepMap identifiers
h1975 <- 'ACH-000587'
h1650 <- 'ACH-000035'
pc14 <- 'ACH-000030'

h1650_delta_ess <- get_delta_ess(h1650, 'LUAD', 'EGFR_mut', check_min_diff = 4)


luad_crispr_raw <- sanger.crispr.dat %>%
  filter(X1 %in% c(h1650, names(h1650_delta_ess))) %>%
  as.data.frame() %>%
  column_to_rownames('X1') %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column('Gene') %>%
  mutate(Gene=map_chr(str_split(Gene, ' \\('), 1))

p1 <- plot_essentialities(h1650, h1975,
                          dat = luad_crispr_raw,
                          label_thres = 6,
                          label_size = 5) +
  scale_y_continuous('Gene essentiality NCI-H1975 (T790M)') +
  scale_x_continuous('Gene essentiality NCI-H1650')

p2 <- plot_essentialities(h1650, pc14,
                          dat = luad_crispr_raw,
                          label_thres = 6,
                          label_size = 5) +
  scale_y_continuous('Gene essentiality PC-14 (sensitive)') +
  scale_x_continuous('Gene essentiality NCI-H1650')


genes_to_highlight <- data.frame(lbls = c('EGFR', 'GRB2',
                                          'SHC1', 'CDC37'),
                                 group = 'EGFR signaling') %>%
  bind_rows(data.frame(lbls = 'GPX4', group = 'GPX4')) %>%
  bind_rows(data.frame(lbls = c('USP18', 'ISG15', 'ADAR', 'EIF41A'),
                       group = 'Interferon α/β signaling'))


to.plot <- h1650_delta_ess %>%
  mutate(is.wt.like=abs(delta_ess_wt) < 3 & wt_comps_over_min < 1/3) %>%
  mutate(lbls=ifelse(abs(delta_ess_mut) > 4 &
                       mut_comps_over_min == 1,
                     Gene, '')) %>%
  left_join(genes_to_highlight, by = 'lbls') %>%
  mutate(group=ifelse(is.na(group), 'notsignificant', group)) %>%
  mutate(group=ifelse(lbls != '' & group == 'notsignificant', 'Other', group))

lbl.cols.df <- data.frame(cols=c(brewer.pal(length(unique(genes_to_highlight$group)),
                                            'Set1'), '#000000'),
                          group=c(unique(genes_to_highlight$group), 'Other')) %>%
  full_join(genes_to_highlight, by = "group") %>%
  full_join(to.plot, by = c('lbls', 'group')) %>%
  mutate(cols=ifelse(is.na(cols), '#000000', as.character(cols))) %>%
  select(cols, lbls) %>%
  filter(lbls != '')

lbl.cols <- lbl.cols.df$cols
names(lbl.cols) <- lbl.cols.df$lbls

diffs_without_wt <- to.plot %>%
  ggplot(aes(x=get(pc14), y=get(h1975),
             col=group)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_abline(slope = 1, intercept = 0,
              linetype = 2, colour = '#00000077') +
  geom_point() +
  geom_text_repel(aes(label=lbls),
                  col=lbl.cols[to.plot$lbls],
                  size=5,
                  point.padding = 0.3,
                  box.padding = 0,
                  min.segment.length = 0,
                  segment.alpha = 0.4,
                  alpha = 0.8) +
  theme_classic2() +
  theme(axis.line = element_blank()) +
  scale_color_manual('',
                     breaks = c(unique(genes_to_highlight$group), 'Other'),
                     values = c(brewer.pal(length(unique(genes_to_highlight$group)),
                                           'Set1'), '#00000022', '#000000')) +
  scale_x_continuous('Δess NCI-H1650 vs PC-14 (Z-score)') +
  scale_y_continuous('Δess NCI-H1650 vs NCI-H1975 (Z-score)') +
  theme(legend.position = c(0.82, 0.15),
        legend.margin = margin(0, 0, 0, 0),
        legend.spacing = unit(0, 'cm'),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.key.height = unit(12, 'pt'))


empty <- ggplot() + theme_void()

diffs_and_comps <- ggarrange(p2, empty, p1, empty, diffs_without_wt,
                             ncol = 5, nrow = 1, widths = c(8, 1, 8, 1, 8))

ggsave(paste0(HOMEDIR, 'manuscript/Figures/figure3i.pdf'),
       diffs_and_comps, height = 5, width = 16.25, device = cairo_pdf)



########################### Suppl Fig 3 ##########################
wt_like_vs_not <- get_delta_ess(h1650, 'LUAD', 'EGFR_mut', check_min_diff = 3) %>%
  mutate(is.wt.like=abs(delta_ess_wt) < 4 & wt_comps_over_min < 1/3) %>%
  mutate(is.wt.like=ifelse(is.wt.like, 'wt-like', 'Not wt-like')) %>%
  mutate(lbls=ifelse(abs(delta_ess_mut) > 4 &
                       mut_comps_over_min == 1,
                     Gene, '')) %>%
  ggplot(aes(x=get(pc14), y=get(h1975),
             col=is.wt.like)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_abline(slope = 1, intercept = 0,
              linetype = 2, colour = '#00000077') +
  geom_point() +
  geom_text_repel(aes(label=lbls),
                  size=5,
                  point.padding = 0.3,
                  box.padding = 0,
                  min.segment.length = 0,
                  segment.alpha = 0.4,
                  alpha = 0.8,
                  show.legend=FALSE) +
  theme_bw() +
  theme(axis.line = element_blank()) +
  scale_color_manual('',
                     labels = c('Not wt-like', 'wt-like'),
                     values = c('#006100', '#b5b5b5'), guide =FALSE) +
  facet_wrap(~is.wt.like, ncol=1) +
  scale_x_continuous('Δess NCI-H1650 vs PC-14 (Z-score)') +
  scale_y_continuous('Δess NCI-H1650 vs NCI-H1975 (Z-score)')

ggsave(paste0(HOMEDIR, 'manuscript/Figures/figureS5a.pdf'),
       wt_like_vs_not, height = 12, width = 6, device = cairo_pdf)


########################## Suppl Fig 4 ##########################

tov21g <- 'ACH-000885'
oaw42 <- 'ACH-000704'


########## TOV-21G
tov21g_delta_ess <- get_delta_ess(tov21g, 'OV', 'PIK3CA_mut', check_min_diff = 3, exclude = oaw42) %>%
  mutate(is.wt.like=abs(delta_ess_wt) < 3 & wt_comps_over_min < 1/3)

tov21g_crispr_raw <- sanger.crispr.dat %>%
  filter(X1 %in% c(tov21g, names(tov21g_delta_ess))) %>%
  as.data.frame() %>%
  column_to_rownames('X1') %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column('Gene') %>%
  mutate(Gene=map_chr(str_split(Gene, ' \\('), 1))

tov_direct_comp <- lapply(names(tov21g_delta_ess)[4:7], function(cl){
  cl.label2 <- model.lut$model_name[which(model.lut$BROAD_ID == cl)]
  plot_essentialities(tov21g, cl,
                      dat = tov21g_crispr_raw,
                      label_thres = 5) +
    scale_y_continuous(paste0('Gene essentiality ', cl.label2)) +
    scale_x_continuous('Gene essentiality TOV-21G')
})

comparisons <- data.frame(i1=rep(1:4, each = 4),
                          i2=rep(1:4, times= 4)) %>%
  filter(i1 < i2)

diff_g_list <- apply(comparisons, 1, function(cc){
  cl1 <- names(tov21g_delta_ess[4:7])[cc[1]]
  cl2<- names(tov21g_delta_ess[4:7])[cc[2]]
  cl.label1 <- model.lut$model_name[which(model.lut$BROAD_ID == cl1)]
  cl.label2 <- model.lut$model_name[which(model.lut$BROAD_ID == cl2)]

  plot_essentiality_diffs(tov21g_delta_ess,
                          diff_val1 = cl1,
                          diff_val2 = cl2,
                          diff_thres = 4) +
    labs(x=paste0('Δess TOV-21G vs ', cl.label1, ' (Z-score)'),
         y=paste0('Δess TOV-21G vs ', cl.label2, ' (Z-score)'))
})

tov.hm <- sanger.crispr.dat %>%
  filter(X1 %in% c(tov21g, names(tov21g_delta_ess))) %>%
  as.data.frame() %>%
  column_to_rownames('X1') %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column('Gene') %>%
  mutate(Gene=map_chr(str_split(Gene, ' \\('), 1)) %>%
  full_join(tov21g_delta_ess, by = 'Gene', suffix = c('', 'diff')) %>%
  filter(abs(delta_ess_mut) > 4 &
           mut_comps_over_min >= 1/2) %>%
  plot_diff_essentiality_heatmap(resis = tov21g, sens = names(tov21g_delta_ess)[4:7])

# lay <- rbind(c(1,1,11,11,11,11),
#              c(2,2,11,11,11,11),
#              c(3,3,11,11,11,11),
#              c(4,4,11,11,11,11),
#              c(5,5,7,7,9,9),
#              c(6,6,8,8,10,10))

lay <- rbind(c(1,11,11,5,8),
             c(1,11,11,5,8),
             c(1,11,11,5,8),
             c(2,11,11,5,8),
             c(2,11,11,6,9),
             c(2,11,11,6,9),
             c(3,11,11,6,9),
             c(3,11,11,6,9),
             c(3,11,11,7,10),
             c(4,11,11,7,10),
             c(4,11,11,7,10),
             c(4,11,11,7,10))

gt <- grid.arrange(grobs = lapply(c(tov_direct_comp, diff_g_list), '+',
                                  theme(plot.margin = unit(c(0.75,0.75,0.75,0.75), "cm"))),
                   layout_matrix = lay)
ggsave(paste0(HOMEDIR, 'manuscript/Figures/figureS3_new.pdf'),
       gt, width = 20, height = 13.33, device = cairo_pdf)
cairo_pdf(paste0(HOMEDIR, 'manuscript/Figures/figureS3_new_hm.pdf'),width = 6, height = 10)
tov.hm
dev.off()



####### OAW-42
oaw42_delta_ess <- get_delta_ess(oaw42, 'OV', 'PIK3CA_mut', check_min_diff = 3, exclude = tov21g) %>%
  mutate(is.wt.like=abs(delta_ess_wt) < 3 & wt_comps_over_min < 1/3)

oaw42_crispr_raw <- sanger.crispr.dat %>%
  filter(X1 %in% c(oaw42, names(oaw42_delta_ess))) %>%
  as.data.frame() %>%
  column_to_rownames('X1') %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column('Gene') %>%
  mutate(Gene=map_chr(str_split(Gene, ' \\('), 1))

oaw42_direct_comp <- lapply(names(oaw42_delta_ess)[4:7], function(cl){
  cl.label2 <- model.lut$model_name[which(model.lut$BROAD_ID == cl)]
  plot_essentialities(oaw42, cl,
                      dat = oaw42_crispr_raw,
                      label_thres = 5) +
    scale_y_continuous(paste0('Gene essentiality ', cl.label2)) +
    scale_x_continuous('Gene essentiality OAW-42')
})

comparisons <- data.frame(i1=rep(1:4, each = 4),
                          i2=rep(1:4, times= 4)) %>%
  filter(i1 < i2)

diff_g_list <- apply(comparisons, 1, function(cc){
  cl1 <- names(oaw42_delta_ess[4:7])[cc[1]]
  cl2<- names(oaw42_delta_ess[4:7])[cc[2]]
  cl.label1 <- model.lut$model_name[which(model.lut$BROAD_ID == cl1)]
  cl.label2 <- model.lut$model_name[which(model.lut$BROAD_ID == cl2)]

  plot_essentiality_diffs(oaw42_delta_ess,
                          diff_val1 = cl1,
                          diff_val2 = cl2,
                          diff_thres = 4) +
    labs(x=paste0('Δess OAW-42 vs ', cl.label1, ' (Z-score)'),
         y=paste0('Δess OAW-42 vs ', cl.label2, ' (Z-score)'))
})

oaw42.hm <- sanger.crispr.dat %>%
  filter(X1 %in% c(oaw42, names(oaw42_delta_ess))) %>%
  as.data.frame() %>%
  column_to_rownames('X1') %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column('Gene') %>%
  mutate(Gene=map_chr(str_split(Gene, ' \\('), 1)) %>%
  full_join(oaw42_delta_ess, by = 'Gene', suffix = c('', 'diff')) %>%
  filter(abs(delta_ess_mut) > 4 &
           mut_comps_over_min >= 1/2) %>%
  plot_diff_essentiality_heatmap(resis = oaw42, sens = names(oaw42_delta_ess)[4:7])

# lay <- rbind(c(1,1,11,11,11,11),
#             c(2,2,11,11,11,11),
#             c(3,3,11,11,11,11),
#             c(4,4,11,11,11,11),
#             c(5,5,7,7,9,9),
#             c(6,6,8,8,10,10))

lay <- rbind(c(1,11,11,5,8),
             c(1,11,11,5,8),
             c(1,11,11,5,8),
             c(2,11,11,5,8),
             c(2,11,11,6,9),
             c(2,11,11,6,9),
             c(3,11,11,6,9),
             c(3,11,11,6,9),
             c(3,11,11,7,10),
             c(4,11,11,7,10),
             c(4,11,11,7,10),
             c(4,11,11,7,10))

go <- grid.arrange(grobs = lapply(c(oaw42_direct_comp, diff_g_list), '+',
                                  theme(plot.margin = unit(c(0.75,0.75,0.75,0.75), "cm"))),
                   layout_matrix = lay)
ggsave(paste0(HOMEDIR, 'manuscript/Figures/figureS4_new.pdf'),
       go,  width = 20, height = 13.33, device = cairo_pdf)

cairo_pdf(paste0(HOMEDIR, 'manuscript/Figures/figureS4_new_hm.pdf'), width = 6, height = 10)
oaw42.hm
dev.off()



########## H1650 heatmap
h1650_delta_ess <- get_delta_ess(h1650, 'LUAD', 'EGFR_mut', check_min_diff = 3) %>%
  mutate(is.wt.like=abs(delta_ess_wt) < 3 & wt_comps_over_min < 1/3)

h1650.hm <- sanger.crispr.dat %>%
  filter(X1 %in% c(h1650, names(h1650_delta_ess))) %>%
  as.data.frame() %>%
  column_to_rownames('X1') %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column('Gene') %>%
  mutate(Gene=map_chr(str_split(Gene, ' \\('), 1)) %>%
  full_join(h1650_delta_ess, by = 'Gene', suffix = c('', 'diff')) %>%
  filter(abs(delta_ess_mut) > 4 &
           mut_comps_over_min == 1) %>%
  plot_diff_essentiality_heatmap(resis = h1650, sens = names(h1650_delta_ess)[4:5])


cairo_pdf(paste0(HOMEDIR, 'manuscript/Figures/figureS5_hm.pdf'), width = 6, height = 10)
h1650.hm
dev.off()


delta_ess_res <- list(h1650_delta_ess, tov21g_delta_ess, oaw42_delta_ess)
save(list = 'delta_ess_res',
     file = paste0(HOMEDIR, 'manuscript/Figures/delta_ess_res.RData'))
