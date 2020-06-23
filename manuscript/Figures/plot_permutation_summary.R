library(tidyverse)
library(patchwork)

HOMEDIR <- '~/Google_Drive/indirect_resistance_stat_framework/'
perm_dir <- paste0(HOMEDIR, 'permutation_test/')


gdsc_perm_anova_sign <- map(dir(paste0(perm_dir, 'gdsc/ANOVA_filtered'), 
                                pattern = 'sign_anova_out.csv',
                                full.names = TRUE),
                            read_csv)

gdsc_perm_anova_sign_numbers <- map_dfr(gdsc_perm_anova_sign, ~data.frame(n=nrow(.x)), .id = 'ID') %>%
  mutate(ID=as.numeric(ID)) %>%
  mutate(perm=dir(paste0(perm_dir, 'gdsc/ANOVA_filtered'), 
                  pattern = 'sign_anova_out.csv')[ID]) %>%
  select(-ID) %>%
  mutate(perm=map_chr(str_split(perm, '_'), 1)) %>%
  full_join(tibble(perm=paste0('P', 1:100)),  by = "perm") %>% # Missing cases mean 0 detected cases 
  replace_na(list(n=0))

gdsc_real_anova_sign <- read_csv('~/Google_Drive/indirect_resistance_stat_framework/results/gdsc/sign_anova_out.csv')

ag <- gdsc_perm_anova_sign_numbers %>%
  ggplot(aes(x=n)) +
  geom_bar(width=0.8, fill = '#1c2b63') +
  theme_bw() +
  scale_x_continuous('Number of significant sensitivity biomarkers after filtering',
                     limits = c(0, 60)) +
  scale_y_continuous(limits = c(0, 17)) +
  geom_vline(xintercept = nrow(gdsc_real_anova_sign),
             color = '#AA0000') +
  annotate('text', x = nrow(gdsc_real_anova_sign), y = Inf, col = '#AA0000',
           hjust = 1.1, vjust = 2, label = 'Number observed\nin real data') + 
  ggtitle(paste0('Sensitivity associations in ',
                 nrow(gdsc_perm_anova_sign_numbers), 
                 ' GDSC permutations'))


gdsc_perm_unres <- map2(dir(paste0(perm_dir, 'gdsc/UNRES_detection'), 
                           pattern = 'high_var_results.RData',
                           full.names = TRUE),
                        dir(paste0(perm_dir, 'gdsc/UNRES_detection'), 
                            pattern = 'fdr.tree.value.RData',
                            full.names = TRUE),
                       function(x, y) {
                         load(x)
                         load(y)
                         high.var.associations %>%
                           filter(outlier.q.value <= 0.15)
                       })


gdsc_perm_unres_numbers <- map_dfr(gdsc_perm_unres, ~.x %>%
                                     select(tissue, association) %>%
                                     distinct() %>% summarise(n=n()), .id = 'ID') %>%
  mutate(ID=as.numeric(ID)) %>%
  mutate(perm=dir(paste0(perm_dir, 'gdsc/UNRES_detection'), 
                  pattern = 'high_var_results.RData')[ID]) %>%
  select(-ID) %>%
  mutate(perm=map_chr(str_split(perm, '_'), 1)) %>%
  full_join(tibble(perm=paste0('P', 1:100)),  by = "perm") %>% # Missing cases mean 0 detected cases 
  replace_na(list(n=0))


load('~/Google_Drive/indirect_resistance_stat_framework/results/gdsc/high_var_results.RData')
load('~/Google_Drive/indirect_resistance_stat_framework/results/gdsc/fdr.tree.value.RData')
gdsc_real_unres <- high.var.associations %>% 
  filter(outlier.q.value <= 0.15)




ug <- gdsc_perm_unres_numbers %>%
  ggplot(aes(x=n)) +
  geom_bar(width=0.8, fill = '#b4772b') +
  theme_bw() +
  scale_x_continuous('Unique associations with UNRES cases identified',
                     limits = c(-0.6, 25)) +
  scale_y_continuous(limits = c(0, 60)) +
  geom_vline(xintercept = gdsc_real_unres %>%
               select(tissue, association) %>%
               distinct() %>% nrow(),
             color = '#AA0000') +
  annotate('text', x = gdsc_real_unres %>%
             select(tissue, association) %>%
             distinct() %>% nrow(), y = Inf, col = '#AA0000',
           hjust = 1.1, vjust = 2, label = 'Number observed\nin real data') + 
  ggtitle(paste0('UNRES discovery in ',
                 nrow(gdsc_perm_unres_numbers), 
                 ' GDSC permutations'))


maxh_gdsc_perm <- gdsc_perm_unres %>%
  map_dfr(~.x %>% group_by(tissue, association) %>% summarise(nmax=max(n.outliers))) %>%
  ggplot(aes(x=nmax)) +
  geom_bar(width=0.8) +
  scale_y_continuous(limits = c(0, 45)) +
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  ggtitle('Max number of UNRES cell lines per association', '100 permutated GDSC datasets')


maxh_gdsc_real <- gdsc_real_unres %>%
  group_by(tissue, association) %>% 
  summarise(nmax=max(n.outliers)) %>%
  ggplot(aes(x=nmax)) +
  geom_bar(width=0.8, fill = '#AA0000') +
  scale_y_continuous(limits = c(0, 45)) +
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  ggtitle('Max number of UNRES cell lines per association', 'Observed GDSC data')



ctrp_perm_anova_sign <- map(dir(paste0(perm_dir, 'ctrp/ANOVA_filtered'), 
                                pattern = 'sign_anova_out.csv',
                                full.names = TRUE),
                            read_csv)

ctrp_perm_anova_sign_numbers <- map_dfr(ctrp_perm_anova_sign, ~data.frame(n=nrow(.x)), .id = 'ID') %>%
  mutate(ID=as.numeric(ID)) %>%
  mutate(perm=dir(paste0(perm_dir, 'ctrp/ANOVA_filtered'), 
                  pattern = 'sign_anova_out.csv')[ID]) %>%
  select(-ID) %>%
  mutate(perm=map_chr(str_split(perm, '_'), 1)) %>%
  full_join(tibble(perm=paste0('P', 1:100)),  by = "perm") %>% # Missing cases mean 0 detected cases 
  replace_na(list(n=0))

ctrp_real_anova_sign <- read_csv('~/Google_Drive/indirect_resistance_stat_framework/results/ctrp/sign_anova_out.csv')

ac <- ctrp_perm_anova_sign_numbers %>%
  ggplot(aes(x=n)) +
  geom_bar(width=0.8, fill = '#1c2b63') +
  theme_bw() +
  scale_x_continuous('Number of significant sensitivity biomarkers after filtering',
                     limits = c(0, 60)) +
  scale_y_continuous(limits = c(0, 17)) +
  geom_vline(xintercept = nrow(ctrp_real_anova_sign),
             color = '#AA0000') +
  annotate('text', x = nrow(ctrp_real_anova_sign), y = Inf, col = '#AA0000',
           hjust = 1.1, vjust = 2, label = 'Number observed\nin real data') + 
  ggtitle(paste0('Sensitivity associations in ',
                 nrow(ctrp_perm_anova_sign_numbers), 
                 ' CTRP permutations'))


ctrp_perm_unres <- map2(dir(paste0(perm_dir, 'ctrp/UNRES_detection'), 
                            pattern = 'high_var_results.RData',
                            full.names = TRUE),
                        dir(paste0(perm_dir, 'ctrp/UNRES_detection'), 
                            pattern = 'fdr.tree.value.RData',
                            full.names = TRUE),
                        function(x, y) {
                          load(x)
                          load(y)
                          high.var.associations %>%
                            filter(outlier.q.value <= 0.15)
                        })

ctrp_perm_unres_numbers <- map_dfr(ctrp_perm_unres, ~.x %>%
                                     select(tissue, association) %>%
                                     distinct() %>% summarise(n=n()), .id = 'ID') %>%
  mutate(ID=as.numeric(ID)) %>%
  mutate(perm=dir(paste0(perm_dir, 'ctrp/UNRES_detection'), 
                  pattern = 'high_var_results.RData')[ID]) %>%
  select(-ID) %>%
  mutate(perm=map_chr(str_split(perm, '_'), 1)) %>%
  full_join(tibble(perm=paste0('P', 1:100)),  by = "perm") %>% # Missing cases mean 0 detected cases 
  replace_na(list(n=0))


load('~/Google_Drive/indirect_resistance_stat_framework/results/ctrp/high_var_results.RData')
load('~/Google_Drive/indirect_resistance_stat_framework/results/ctrp/fdr.tree.value.RData')
ctrp_real_unres <- high.var.associations %>% 
  filter(outlier.q.value <= 0.15)




uc <- ctrp_perm_unres_numbers %>%
  ggplot(aes(x=n)) +
  geom_bar(width=0.8, fill = '#b4772b') +
  theme_bw() +
  scale_x_continuous('Unique associations with UNRES cases identified',
                     limits = c(-0.6, 25)) +
  scale_y_continuous(limits = c(0, 60)) +
  geom_vline(xintercept = ctrp_real_unres %>%
               select(tissue, association) %>%
               distinct() %>% nrow(),
             color = '#AA0000') +
  annotate('text', x = ctrp_real_unres %>%
             select(tissue, association) %>%
             distinct() %>% nrow(), y = Inf, col = '#AA0000',
           hjust = 1.1, vjust = 2, label = 'Number observed\nin real data') + 
  ggtitle(paste0('UNRES discovery in ',
                 nrow(ctrp_perm_unres_numbers), 
                 ' CTRP permutations'))


maxh_ctrp_perm <- ctrp_perm_unres %>%
  map_dfr(~.x %>% group_by(tissue, association) %>% summarise(nmax=max(n.outliers))) %>%
  ggplot(aes(x=nmax)) +
  geom_bar(width=0.8) +
  scale_y_continuous(limits = c(0, 45)) +
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  ggtitle('Max number of UNRES cell lines per association', '100 permutated CTRP datasets')


maxh_ctrp_real <- ctrp_real_unres %>%
  group_by(tissue, association) %>% 
  summarise(nmax=max(n.outliers)) %>%
  ggplot(aes(x=nmax)) +
  geom_bar(width=0.8, fill = '#AA0000') +
  scale_y_continuous(limits = c(0, 45)) +
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  ggtitle('Max number of UNRES cell lines per association', 'Observed CTRP data')


# perm_summ <- (ag | ac)/(ug | uc)
# ggsave(perm_summ, filename = paste0(HOMEDIR, 'manuscript/Figures/perm_summary_fig.pdf'),
#        width = 11, height = 4.5)
# 
# 
# max_cl_per_assoc <- (maxh_gdsc_perm | maxh_ctrp_perm) / (maxh_gdsc_real | maxh_ctrp_real)
# ggsave(max_cl_per_assoc, filename = paste0(HOMEDIR, 'manuscript/Figures/perm_summary_max_cl.pdf'),
#        width = 11, height = 4.5)

allp <- ggarrange(ggarrange(ag, ac, ug, uc, nrow = 2, ncol = 2, 
                    labels = c('A', 'B', 'C', 'D')),
          ggarrange(maxh_gdsc_perm, maxh_ctrp_perm,
                    maxh_gdsc_real, maxh_ctrp_real, nrow = 2, ncol = 2,
                    labels = c('E', 'F', 'G', 'H')),
          nrow = 2)

ggsave(allp, filename = paste0(HOMEDIR, 'manuscript/Figures/figureS4_perm_summary.pdf'),
       width = 11, height = 9)
