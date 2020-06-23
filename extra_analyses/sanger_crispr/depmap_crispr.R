library(tidyverse)
library(ggrepel)
library(ggpubr)
library(RColorBrewer)
library(ComplexHeatmap)
library(GGally)
library(circlize)

HOMEDIR <- '~/Google_Drive/indirect_resistance_stat_framework/'

sanger.crispr.dat <- read_csv(paste0(HOMEDIR, 'extra_analyses/sanger_crispr/gene_effect.csv'))
model.lut <- read_csv(paste0(HOMEDIR, 'extra_analyses/sanger_crispr/model_list_latest.csv.gz'))


plot_essentialities <- function(cl1, cl2, 
                                label_thres = 5,
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
                    size = 5,
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
                                    diff_thres = 3,
                                    highlight_data=NA,
                                    group.color=!all(is.na(highlight_data)),
                                    diff.z.colour=FALSE,
                                    highlight_wt_like=FALSE) {
  
  if (group.color) {
    if (highlight_wt_like){
      g <-  diff_data %>%
        mutate(lbls=ifelse(abs(get(diff_val1)) > diff_thres & 
                             abs(get(diff_val2)) > diff_thres,
                           Gene.nice, '')) %>%
        left_join(highlight_data, by = 'lbls') %>%
        mutate(group=ifelse(is.na(group), 'notsignificant', group)) %>%
        mutate(group=ifelse(lbls != '' & group == 'notsignificant', 'Other', group)) %>%
        ggplot(aes(x=get(diff_val1), y=get(diff_val2),
                   col=is.wt.like, fill=group)) +
        scale_color_manual('',
                           values = c('#ffbd2366', '#00000044'),
                           guide = FALSE) +
        scale_fill_manual('', 
                          breaks = c(unique(highlight_data$group), 'Other'),
                          values = c(brewer.pal(length(unique(highlight_data$group)),
                                                'Set1'), '#00000022', '#000000'))
      
    } else {
      g <-  diff_data %>%
        mutate(lbls=ifelse(abs(get(diff_val1)) > diff_thres & 
                             abs(get(diff_val2)) > diff_thres,
                           Gene.nice, '')) %>%
        left_join(highlight_data, by = 'lbls') %>%
        mutate(group=ifelse(is.na(group), 'notsignificant', group)) %>%
        mutate(group=ifelse(lbls != '' & group == 'notsignificant', 'Other', group)) %>%
        ggplot(aes(x=get(diff_val1), y=get(diff_val2), col=group)) +
        scale_colour_manual('', 
                            breaks = c(unique(highlight_data$group), 'Other'),
                            values = c(brewer.pal(length(unique(highlight_data$group)),
                                                  'Set1'), '#00000022', '#000000'))
    }
  } else if (diff.z.colour){ 
    g <-  diff_data %>%
      mutate(lbls=ifelse(abs(get(diff_val1)) > diff_thres & 
                           abs(get(diff_val2)) > diff_thres,
                         Gene.nice, '')) %>%
      ggplot(aes(x=get(diff_val1), y=get(diff_val2), 
                 col = mean.diff.zs,alpha=(lbls==''))) +
      scale_color_gradient2(low='#0000CC', mid='#999999', high = '#CC0000') +
      scale_alpha_manual('', 
                         values = c(1, 0.25),
                         guide = FALSE)
  } else if (highlight_wt_like) {
    g <-  diff_data %>%
      mutate(lbls=ifelse(abs(get(diff_val1)) > diff_thres & 
                           abs(get(diff_val2)) > diff_thres,
                         Gene.nice, '')) %>%
      ggplot(aes(x=get(diff_val1), y=get(diff_val2), col=is.wt.like,
                 fill=(lbls==''))) +
      scale_color_manual('',
                         values = c('#ffbd2366', '#9600f866'),
                         guide = FALSE) +
      scale_fill_manual('', 
                        values = c('#000000', '#00000022'),
                        guide = FALSE)
  } else {
    g <-  diff_data %>%
      mutate(lbls=ifelse(abs(get(diff_val1)) > diff_thres & 
                           abs(get(diff_val2)) > diff_thres,
                         Gene.nice, '')) %>%
      ggplot(aes(x=get(diff_val1), y=get(diff_val2), col=(lbls==''))) +
      scale_colour_manual('', 
                          values = c('#000000', '#00000022'),
                          guide = FALSE)
    
  }
  g <- g +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    geom_abline(slope = 1, intercept = 0,
                linetype = 2, colour = '#00000077') +
    geom_point(shape=21) +
    geom_text_repel(aes(label=lbls),
                    col='#000000',
                    size=5,
                    point.padding = 0.3,
                    box.padding = 0,
                    min.segment.length = 0,
                    segment.alpha = 0.4,
                    alpha = 0.8,
                    show.legend = FALSE) +
    
    theme_classic2() +
    theme(axis.line = element_blank()) 
}


# DepMap identifiers
h1975 <- 'ACH-000587' 
h1650 <- 'ACH-000035'
pc14 <- 'ACH-000030'

pc3 <- 'ACH-002184'
h3255 <- 'ACH-000109'
hcc827 <- 'ACH-000012'

all.luad.egfr <- c(h1975, h1650, pc14,
                   pc3, h3255, hcc827)


sanger.crispr.tidy.luad <- sanger.crispr.dat %>%
  filter(X1 %in% all.luad.egfr) %>%
  as.data.frame() %>%
  column_to_rownames('X1') %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column('Gene')

p1 <- plot_essentialities(h1650, h1975,
                          dat = sanger.crispr.tidy.luad,
                          label_thres = 6) +
  scale_y_continuous('Gene essentiality NCI-H1975 (T790M)') +
  scale_x_continuous('Gene essentiality NCI-H1650')

p2 <- plot_essentialities(h1650, pc14,
                          dat = sanger.crispr.tidy.luad,
                          label_thres = 6) +
  scale_y_continuous('Gene essentiality PC-14 (sensitive)') +
  scale_x_continuous('Gene essentiality NCI-H1650')

# both_indiv_comp <- ggarrange(p1, p2, ncol = 2)
# ggsave(paste0(HOMEDIR, 'extra_analyses/sanger_crispr/h1650_vs_h1975_and_pc4_crispr.pdf'),
#        both_indiv_comp,   height = 6, width = 16)
h1650.diff.score <- sanger.crispr.tidy.luad %>%
  mutate(diff1=get(h1650) - get(pc14),
         diff2=get(h1650) - get(h1975),
         diff1.z=(diff1-mean(diff1))/sd(diff1),
         diff2.z=(diff2-mean(diff2))/sd(diff2),
         avg.diff=(diff1.z+diff2.z)/2) %>%
  arrange(avg.diff) %>%
  mutate(Gene.nice=map_chr(str_split(Gene, ' \\('), 1)) %>%
  mutate(mean.diff.zs=(diff1.z+diff2.z)/2)



# h1650.diff.score %>%
#   arrange(-mean.diff.zs) %>%
#   select(Gene.nice, mean.diff.zs) %>%
#   write_tsv(path=paste0(HOMEDIR, 'extra_analyses/sanger_crispr/h1650_diff_ess_ranked.rnk'),
#             col_names = FALSE)

genes_to_highlight <- data.frame(lbls = c('EGFR', 'GRB2',
                                          'SHC1', 'CDC37'),
                                 group = 'EGFR signaling') %>%
  # bind_rows(data.frame(lbls = c('SKA1', 'SKA2', 'CENPM', 'CENPW',
  #                               'PPP1R12A'),
  #                      group = 'Kinetochore')) %>%
  bind_rows(data.frame(lbls = c('USP18', 'ISG15', 'ADAR', 'EIF41A'),
                       group = 'Interferon signaling')) %>%
  # bind_rows(data.frame(lbls = c('ALDOA', 'ENO1'),
  #                      group = 'Gluconegenesis')) %>%
  bind_rows(data.frame(lbls = 'GPX4', group = 'GPX4'))

depmap_diffs <- h1650.diff.score %>%
  plot_essentiality_diffs(diff_val1 = 'diff1.z', diff_val2 = 'diff2.z',
                          highlight_data = genes_to_highlight, diff_thres = 4) +
  scale_x_continuous('Δess NCI-H1650 vs PC-14 (Z-score)') +
  scale_y_continuous('Δess NCI-H1650 vs NCI-H1975 (Z-score)') +
  theme(legend.position = c(0.85, 0.15))



# depmap_diffs <- h1650.diff.score %>%
#   plot_essentiality_diffs(diff_val1 = 'diff1.z', diff_val2 = 'diff2.z',
#                           diff.z.colour = TRUE, diff_thres = 4) +
#   scale_x_continuous('Δ essentiality NCI-H1650 vs PC-14 (Z-score)') +
#   scale_y_continuous('Δ essentiality NCI-H1650 vs NCI-H1975 (Z-score)')

empty <- ggplot() + theme_void()

diffs_and_comps <- ggarrange(p2, empty, p1, empty, depmap_diffs,
                             ncol = 5, nrow = 1, widths = c(8, 1, 8, 1, 8))

ggsave(paste0(HOMEDIR, 'manuscript/Figures/figure3i.pdf'),
       diffs_and_comps, height = 5, width = 16.25, device = cairo_pdf)

# library(plotly)
# plot_ly(h1650.diff.score,
#         x=~diff1.z, y=~diff2.z,
#         color=~abs(avg.diff),
#         text=~Gene.nice) %>%
#   add_markers(hoverinfo='text')
# 
# 
# h1650.diff.score %>%
#   filter(diff1.z > 4 & diff2.z > 4) %>%
#   select(Gene.nice) %>%
#   unlist %>% unname %>% cat(sep='\n')
# 
# h1650.diff.score %>%
#   filter(diff1.z < -4 & diff2.z < -4) %>%
#   select(Gene.nice) %>%
#   unlist %>% unname %>% cat(sep='\n')


### Same for OV case

# DepMap identifiers
tov21g <- 'ACH-000885'
oaw42 <- 'ACH-000704'

oc314 <- 'ACH-000962'
ovmiu <- 'ACH-002183'
ovise <- 'ACH-000527'


all.ov.pik3ca <- c(tov21g, oaw42, oc314,
                   ovmiu, ovise)

sanger.crispr.tidy.ov <- sanger.crispr.dat %>%
  filter(X1 %in% all.ov.pik3ca) %>%
  as.data.frame() %>%
  column_to_rownames('X1') %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column('Gene')


### TOV-21G

tov21g_oc314 <- plot_essentialities(tov21g, oc314,
                                    dat = sanger.crispr.tidy.ov) +
  scale_y_continuous('Gene essentiality OC-314') +
  scale_x_continuous('Gene essentiality TOV-21G')

tov21g_ovise <- plot_essentialities(tov21g, ovise,
                                    dat = sanger.crispr.tidy.ov) +
  scale_y_continuous('Gene essentiality OVISE') +
  scale_x_continuous('Gene essentiality TOV-21G')

tov21g_ovmiu <- plot_essentialities(tov21g, ovmiu,
                                    dat = sanger.crispr.tidy.ov) +
  scale_y_continuous('Gene essentiality OVMIU') +
  scale_x_continuous('Gene essentiality TOV-21G')


tov21g.diff.score <- sanger.crispr.tidy.ov %>%
  mutate(vsoc314=get(tov21g) - get(oc314),
         vsoc314.z=(vsoc314-mean(vsoc314))/sd(vsoc314),
         vsovise=get(tov21g) - get(ovise),
         vsovise.z=(vsovise-mean(vsovise))/sd(vsovise),
         vsovmiu=get(tov21g) - get(ovmiu),
         vsovmiu.z=(vsovmiu-mean(vsovmiu))/sd(vsovmiu)) %>%
  mutate(Gene.nice=map_chr(str_split(Gene, ' \\('), 1)) %>%
  mutate(mean.diff.zs=(vsoc314.z+vsovise.z+vsovmiu.z)/3)

comparisons <- rbind(c('vsoc314.z', 'vsovise.z'),
                     c('vsoc314.z', 'vsovmiu.z'),
                     c('vsovise.z', 'vsovmiu.z'))
diff_g_list <- apply(comparisons, 1, function(cc){
  plot_essentiality_diffs(tov21g.diff.score,
                          diff_val1 = cc[1],
                          diff_val2 = cc[2],
                          diff_thres = 4,
                          diff.z.colour = FALSE,
                          group.color = FALSE)
})

d1 <- diff_g_list[[1]] + labs(x='Δess TOV-21G vs OC-314 (Z-score)',
                              y='Δess TOV-21G vs OVISE (Z-score)')
d2 <- diff_g_list[[2]] + labs(x='Δess TOV-21G vs OC-314 (Z-score)',
                              y='Δess TOV-21G vs OVMIU (Z-score)')
d3 <- diff_g_list[[3]] + labs(x='Δess TOV-21G vs OVISE (Z-score)',
                              y='Δess TOV-21G vs OVMIU (Z-score)')

empty <- ggplot() + theme_void()

gt <- ggarrange(tov21g_oc314, empty, empty,
                d1, tov21g_ovise, empty,
                d2, d3, tov21g_ovmiu,
                labels = c('A', '', '',
                           'B', 'C', '',
                           'D', 'E', 'F'),
                ncol=3, nrow=3)

ggsave(paste0(HOMEDIR, 'manuscript/Figures/figureS3.pdf'),
       gt, width = 14, height = 14, device = cairo_pdf)



## OAW-42

oaw42_oc314 <- plot_essentialities(oaw42, oc314,
                                   dat = sanger.crispr.tidy.ov) +
  scale_y_continuous('Gene essentiality OC-314') +
  scale_x_continuous('Gene essentiality OAW-42')

oaw42_ovise <- plot_essentialities(oaw42, ovise,
                                   dat = sanger.crispr.tidy.ov) +
  scale_y_continuous('Gene essentiality OVISE') +
  scale_x_continuous('Gene essentiality OAW-42')

oaw42_ovmiu <- plot_essentialities(oaw42, ovmiu,
                                   dat = sanger.crispr.tidy.ov) +
  scale_y_continuous('Gene essentiality OVMIU') +
  scale_x_continuous('Gene essentiality OAW-42')


oaw42.diff.score <- sanger.crispr.tidy.ov %>%
  mutate(vsoc314=get(oaw42) - get(oc314),
         vsoc314.z=(vsoc314-mean(vsoc314))/sd(vsoc314),
         vsovise=get(oaw42) - get(ovise),
         vsovise.z=(vsovise-mean(vsovise))/sd(vsovise),
         vsovmiu=get(oaw42) - get(ovmiu),
         vsovmiu.z=(vsovmiu-mean(vsovmiu))/sd(vsovmiu)) %>%
  mutate(Gene.nice=map_chr(str_split(Gene, ' \\('), 1)) %>%
  mutate(mean.diff.zs=(vsoc314.z+vsovise.z+vsovmiu.z)/3)

comparisons <- rbind(c('vsoc314.z', 'vsovise.z'),
                     c('vsoc314.z', 'vsovmiu.z'),
                     c('vsovise.z', 'vsovmiu.z'))
diff_g_list <- apply(comparisons, 1, function(cc){
  plot_essentiality_diffs(oaw42.diff.score,
                          diff_val1 = cc[1],
                          diff_val2 = cc[2],
                          diff_thres = 4,
                          group.color = FALSE)
})

d1 <- diff_g_list[[1]] + labs(x='Δess OAW-42 vs OC-314 (Z-score)',
                              y='Δess OAW-42 vs OVISE (Z-score)') +
  lims(x=c(-7, 7), y=c(-7, 7))
d2 <- diff_g_list[[2]] + labs(x='Δess OAW-42 vs OC-314 (Z-score)',
                              y='Δess OAW-42 vs OVMIU (Z-score)')
d3 <- diff_g_list[[3]] + labs(x='Δess OAW-42 vs OVISE (Z-score)',
                              y='Δess OAW-42 vs OVMIU (Z-score)')

empty <- ggplot() + theme_void()

go <- ggarrange(oaw42_oc314, empty, empty,
                d1, oaw42_ovise, empty,
                d2, d3, oaw42_ovmiu,
                labels = c('A', '', '',
                           'B', 'C', '',
                           'D', 'E', 'F'),
                ncol=3, nrow=3)

ggsave(paste0(HOMEDIR, 'manuscript/Figures/figureS4.pdf'),
       go, width = 14, height = 14, device = cairo_pdf)

# oaw42.diff.score %>%
#   filter(vsoc314.z > 3 & 
#            vsovise.z > 3 &
#            vsovmiu.z > 3) %>%
#   select(Gene.nice) %>%
#   unlist %>% unname %>% cat(sep='\n')
# 
# oaw42.diff.score %>%
#   filter(vsoc314.z < -3 & 
#            vsovise.z < -3 &
#            vsovmiu.z < -3) %>%
#   select(Gene.nice) %>%
#   unlist %>% unname %>% cat(sep='\n')


#######################################################################

plot_diff_essentiality_heatmap <- function(diffs, resis, sens) {
  
  mycol <- colorRamp2(c(-2, 0, 2), hcl.colors(3, 'Blue-Red'))
  cl.labels <- model.lut$model_name[sapply(c(resis, sens), 
                                           function(x) which(model.lut$BROAD_ID == x))]
  
  all.mean.diffs <- diffs$mean.diff.zs
  
  sidecols <- colorRamp2(c(min(all.mean.diffs), 0, max(all.mean.diffs)), 
                         c('#E57E00','#FFFFFF', '#7E00E5'))
  
  topcols <- c(rep('#888888', length(resis)),
               rep('#000000', length(sens)))
  
  topcols <- c(rep('Resistant', length(resis)),
               rep('Sensitive', length(sens)))
  
  directioncol <- diffs$mean.diff.zs > 0
  
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
                left_annotation = rowAnnotation(Difference = diffs$mean.diff.zs,
                                                col = list(Difference=sidecols),
                                                show_annotation_name = FALSE),
                show_parent_dend_line = FALSE,
                border = TRUE)
  draw(hh, merge_legend = TRUE)
}


# pca.ov <- sanger.crispr.tidy.ov %>%
#   column_to_rownames('Gene') %>%
#   as.matrix() %>%
#   prcomp()
# 
# pca.ov$rotation %>%
#   as.data.frame() %>%
#   rownames_to_column('sample') %>%
#   mutate(cl = sapply(sample, function(x) model.lut$model_name[which(model.lut$BROAD_ID == x)])) %>%
#   ggplot(aes(x=PC1, y=PC2)) +
#   geom_text(aes(label = cl))





pdf(paste0(HOMEDIR, 'extra_analyses/sanger_crispr/heatmap_oaw42.pdf'),
    height = 8, width = 5)
oaw42.diff.score %>%
  arrange(-abs(mean.diff.zs)) %>%
  filter(abs(mean.diff.zs) > 4) %>%
  plot_diff_essentiality_heatmap(oaw42, all.ov.pik3ca[3:5])
dev.off()

pdf(paste0(HOMEDIR, 'extra_analyses/sanger_crispr/heatmap_tov21g.pdf'),
    height = 8, width = 5)
tov21g.diff.score %>%
  arrange(-abs(mean.diff.zs)) %>%
  filter(abs(mean.diff.zs) > 4) %>%
  plot_diff_essentiality_heatmap(tov21g, all.ov.pik3ca[3:5])
dev.off()


pdf(paste0(HOMEDIR, 'extra_analyses/sanger_crispr/heatmap_h1650.pdf'),
    height = 8, width = 5)
h1650.diff.score %>%
  arrange(-abs(mean.diff.zs)) %>%
  filter(abs(mean.diff.zs) > 4) %>%
  plot_diff_essentiality_heatmap(h1650, c(h1975, pc14))
dev.off()


### Trying something different
load(paste0(HOMEDIR, 'data/combined_data.RData'))

get_delta_ess <- function(outlier, tissue, alt, exclude = NA, check_min_diff = 2) {
  
  all.cosmic.id <- all.bems.tidy %>%
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


h1650 <- 'ACH-000035'
library(plotly)
vline <- function(x = 0, color = "black") {
  list(
    type = "line", 
    y0 = 0, 
    y1 = 1, 
    yref = "paper",
    x0 = x, 
    x1 = x, 
    line = list(color = color, dash = 'dot')
  )
}

hline <- function(y = 0, color = "black") {
  list(
    type = "line", 
    x0 = 0, 
    x1 = 1, 
    xref = "paper",
    y0 = y, 
    y1 = y, 
    line = list(color = color, dash = 'dot')
  )
}



h1650_delta_ess <- get_delta_ess(h1650, 'LUAD', 'EGFR_mut', check_min_diff = 3)

h1650_delta_ess %>%
  # filter(abs(delta_ess_wt) < 3 & abs(delta_ess_mut) > 3) %>%
  plot_ly(x=~delta_ess_wt, y=~delta_ess_mut, text=~Gene) %>%
  add_markers(hoverinfo ='text') %>%
  layout(shapes = list(hline(3), hline(-3), vline(3), vline(-3)),
         title = 'H1650')

h1650.wt.like <- h1650_delta_ess %>% 
  mutate(is.wt.like=abs(delta_ess_wt) < 3) %>%
  rename(Gene.nice = Gene)

## TODO
h1650.diff.score %>%
  full_join(h1650.wt.like, by = 'Gene.nice') %>%
  plot_essentiality_diffs(diff_val1 = 'diff1.z', diff_val2 = 'diff2.z',
                          highlight_data = genes_to_highlight,
                          highlight_wt_like = TRUE, diff_thres = 4) +
  # geom_point(aes(col=is.wt.like))
  scale_x_continuous('Δess NCI-H1650 vs PC-14 (Z-score)') +
  scale_y_continuous('Δess NCI-H1650 vs NCI-H1975 (Z-score)') +
  theme(legend.position = c(0.85, 0.15))



get_delta_ess(h1650, 'LUAD', 'EGFR_mut') %>%
  filter(abs(delta_ess_wt) < 3 & (delta_ess_mut) <  -3) %>%
  select(Gene) %>%
  unlist() %>% as.character() %>% cat(sep='\n')

get_delta_ess(h1650, 'LUAD', 'EGFR_mut') %>%
  filter(abs(delta_ess_wt) < 3 & (delta_ess_mut) > 3) %>%
  select(Gene) %>%
  unlist() %>% as.character() %>% cat(sep='\n')


tov21g <- 'ACH-000885'
oaw42 <- 'ACH-000704'

get_delta_ess(tov21g, 'OV', 'PIK3CA_mut', exclude = oaw42) %>%
  # filter(abs(delta_ess_wt) < 3 & abs(delta_ess_mut) > 3) %>%
  plot_ly(x=~delta_ess_wt, y=~delta_ess_mut, text=~Gene) %>%
  add_markers(hoverinfo ='text') %>%
  layout(shapes = list(hline(3), hline(-3), vline(3), vline(-3)),
         title = 'TOV-21G')


get_delta_ess(tov21g, 'OV', 'PIK3CA_mut', exclude = oaw42) %>%
  filter(abs(delta_ess_wt) < 3 & (delta_ess_mut) > 3) %>%
  select(Gene) %>%
  unlist() %>% as.character() %>% cat(sep='\n')

get_delta_ess(tov21g, 'OV', 'PIK3CA_mut', exclude = oaw42) %>%
  filter(abs(delta_ess_wt) < 3 & (delta_ess_mut) < -3) %>%
  select(Gene) %>%
  unlist() %>% as.character() %>% cat(sep='\n')


get_delta_ess(oaw42, 'OV', 'PIK3CA_mut', exclude = tov21g) %>%
  # filter(abs(delta_ess_wt) < 3 & abs(delta_ess_mut) > 3) %>%
  plot_ly(x=~delta_ess_wt, y=~delta_ess_mut, text=~Gene) %>%
  add_markers(hoverinfo ='text') %>%
  layout(shapes = list(hline(3), hline(-3), vline(3), vline(-3)),
         title= 'OAW-42')



get_delta_ess(oaw42, 'OV', 'PIK3CA_mut', exclude = tov21g) %>%
  filter(abs(delta_ess_wt) < 3 & (delta_ess_mut) > 3) %>%
  select(Gene) %>%
  unlist() %>% as.character() %>% cat(sep='\n')


get_delta_ess(oaw42, 'OV', 'PIK3CA_mut', exclude = tov21g) %>%
  filter(abs(delta_ess_wt) < 3 & (delta_ess_mut) < -3) %>%
  select(Gene) %>%
  unlist() %>% as.character() %>% cat(sep='\n')





# Redo Figure 3I
genes_to_highlight <- data.frame(lbls = c('EGFR', 'GRB2',
                                          'SHC1', 'CDC37'),
                                 group = 'EGFR signaling') %>%
  # bind_rows(data.frame(lbls = c('SKA1', 'SKA2', 'CENPM', 'CENPW',
  #                               'PPP1R12A'),
  #                      group = 'Kinetochore')) %>%
  # bind_rows(data.frame(lbls = c('ALDOA', 'ENO1'),
  #                      group = 'Gluconegenesis')) %>%
  bind_rows(data.frame(lbls = 'GPX4', group = 'GPX4')) %>%
  bind_rows(data.frame(lbls = c('USP18', 'ISG15', 'ADAR', 'EIF41A'),
                       group = 'Interferon α/β signaling'))
diff_thres <- 4

h1650_delta_ess <- get_delta_ess(h1650, 'LUAD', 'EGFR_mut', check_min_diff = 4)

to.plot <- h1650_delta_ess %>%
  mutate(is.wt.like=abs(delta_ess_wt) < 3 & wt_comps_over_min < 1/3) %>%
  mutate(lbls=ifelse(abs(delta_ess_mut) > diff_thres &
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

diffs_with_wt <- to.plot %>%
  ggplot(aes(x=get(pc14), y=get(h1975),
             col=is.wt.like, fill=group)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_abline(slope = 1, intercept = 0,
              linetype = 2, colour = '#00000077') +
  geom_point(shape=21) +
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
                     labels = c('Not wt-like', 'wt-like'),
                     values = c('#FFA500', '#00000044')) +
  scale_fill_manual('', 
                    breaks = c(unique(genes_to_highlight$group), 'Other'),
                    values = c(brewer.pal(length(unique(genes_to_highlight$group)),
                                          'Set1'), '#00000011', '#000000')) +
  scale_x_continuous('Δess NCI-H1650 vs PC-14 (Z-score)') +
  scale_y_continuous('Δess NCI-H1650 vs NCI-H1975 (Z-score)') +
  theme(legend.position = c(0.82, 0.15),
        legend.margin = margin(0, 0, 0, 0), 
        legend.spacing = unit(0, 'cm'), 
        legend.background = element_blank(),
        legend.title = element_blank(), 
        legend.key.height = unit(12, 'pt'))

ggsave(paste0(HOMEDIR, 'manuscript/Figures/figure3i_with_wt.pdf'),
       diffs_with_wt, height = 5, width = 5, device = cairo_pdf)

