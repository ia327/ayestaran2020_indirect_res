library(tidyverse)
library(scales)
library(openxlsx)

HOMEDIR <- '../../'

#### Suppl Table 1: GDSC ANOVA results
anova_out_files <- dir(paste0(HOMEDIR, 'results/gdsc/ANOVA_out'),
                       pattern = '_anova_out.csv',
                       full.names = TRUE)

gdsc_anova_out <- map(anova_out_files, read_csv)
names(gdsc_anova_out) <- map_chr(str_split(map_chr(str_split(anova_out_files, '/'), tail, n=1), '_'), 1)

load(paste0(HOMEDIR, 'results/gdsc/sign_anova_out.RData'))
after_filter <- selected.sign.anovas.df %>%
  select(association, tissue) %>%
  mutate(filters_pass='Y')

gdsc_anova_out <- map2(gdsc_anova_out, names(gdsc_anova_out), ~mutate(.x, tissue = .y)) %>%
  bind_rows() %>%
  full_join(after_filter, by = c("association", "tissue")) %>%
  replace_na(list(filters_pass = 'N')) %>%
  split(.$tissue)

write.xlsx(gdsc_anova_out, file =  paste0(HOMEDIR, 'manuscript/Figures/tableS1_gdsc.xlsx'))


#### Suppl Table 2: CTRP ANOVA results
anova_out_files <- dir(paste0(HOMEDIR, 'results/ctrp/ANOVA_out'),
                       pattern = '_anova_out.csv',
                       full.names = TRUE)

ctrp_anova_out <- map(anova_out_files, read_csv)
names(ctrp_anova_out) <- map_chr(str_split(map_chr(str_split(anova_out_files, '/'), tail, n=1), '_'), 1)

load(paste0(HOMEDIR, 'results/ctrp/sign_anova_out.RData'))
after_filter <- selected.sign.anovas.df %>%
  select(association, tissue) %>%
  mutate(filters_pass='Y')

ctrp_anova_out <- map2(ctrp_anova_out, names(ctrp_anova_out), ~mutate(.x, tissue = .y)) %>%
  bind_rows() %>%
  full_join(after_filter, by = c("association", "tissue")) %>%
  replace_na(list(filters_pass = 'N')) %>%
  split(.$tissue)

write.xlsx(ctrp_anova_out, file =  paste0(HOMEDIR, 'manuscript/Figures/tableS2_ctrp.xlsx'))



#### Suppl Table 3: GDSC UNRES
gdsc_put_res <- read_tsv(paste0(HOMEDIR, 'results/gdsc/putative_res_markers.tsv'))

gdsc_put_res %>%
  rename(`UNRES FDR` = `Outlier FDR`,
         `UNRES p-value` = `Outlier p-value`,
         `UNRES cell lines` = `Outlier cell lines`,
         `Normalised SD decrease` = `Effect size`) %>%
  mutate(`UNRES p-value`= scientific(`UNRES p-value`),
         `UNRES FDR` = scientific(`UNRES FDR`),
         `Normalised SD decrease` = round(`Normalised SD decrease`, 3)) %>%
  write_csv(path = paste0(HOMEDIR, 'manuscript/Figures/tableS3_gdsc.csv'))

#### Suppl Table 4: CTRP UNRES
ctrp_put_res <- read_tsv(paste0(HOMEDIR, 'results/ctrp/putative_res_markers.tsv'))

ctrp_put_res %>%
  rename(`UNRES FDR` = `Outlier FDR`,
         `UNRES p-value` = `Outlier p-value`,
         `UNRES cell lines` = `Outlier cell lines`,
         `Normalised SD decrease` = `Effect size`) %>%
  mutate(`UNRES p-value`= scientific(`UNRES p-value`),
         `UNRES FDR` = scientific(`UNRES FDR`),
         `Normalised SD decrease` = round(`Normalised SD decrease`, 3)) %>%
  write_csv(path = paste0(HOMEDIR, 'manuscript/Figures/tableS4_ctrp.csv'))


### Suppl Table 5: Differential essentiality results
load(paste0(HOMEDIR, 'manuscript/Figures/delta_ess_res.RData'))

delta_res <- map(delta_ess_res, function(l){
  names_l <- names(l)
  names_l[grepl('ACH-', names_l)] <- paste0(names_l[grepl('ACH-', names_l)], '_diff')
  names(l) <- names_l
  l$mut_comps_over_min <- NULL
  l$wt_comps_over_min <- NULL
  l$is.wt.like <- NULL
  return(l)
})


names(delta_res) <- c('NCI-H1650_LUAD', 'TOV-21G_OV', 'OAW-42_OV')
write.xlsx(delta_res, file =  paste0(HOMEDIR, 'manuscript/Figures/tableS5_delta_ess.xlsx'))



### Suppl Table 6: CTRP refitted data
load(paste0(HOMEDIR, 'data/combined_for_pipeline/ctrp_dr_data.RData'))

all.ctrp %>%
  select(CL, COSMIC_ID, experiment_id, GDSC_cell_line, TCGA_DESC,
         DRUG_ID, DRUG_NAME, PUTATIVE_TARGET, MAX_CONC, LN_IC50, auc, RMSE) %>%
  rename(CCLE_cell_ID = CL) %>%
  write.xlsx(file =  paste0(HOMEDIR, 'manuscript/Figures/tableS6_ctrp_refit.xlsx'))
