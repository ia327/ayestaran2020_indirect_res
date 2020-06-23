library(gdscIC50)
library(tidyverse)

HOMEDIR <- '../../'
DATADIR <- paste0(HOMEDIR,'/data/')

load(paste0(HOMEDIR, 'data/drug/ctrp_all_fitted.RData'))

ctrp_model_stats %>%
  ggplot(aes(x=IC50)) +
  geom_histogram(bins=200)

ctrp_model_stats %>%
  ggplot(aes(x=auc)) +
  geom_histogram(bins=200)

plotResponse(ctrp_model_stats,
             cell_line = 272,
             drug_identifier = '362338_37')







FITTED_DR_1 <- read_csv(paste0(DATADIR, "drug/GDSC1_fitted_dose_response_17Jul19.csv"), col_types = cols())
FITTED_DR_2 <- read_csv(paste0(DATADIR, "drug/GDSC2_fitted_dose_response_17Jul19.csv"), col_types = cols())
FITTED_DR <- bind_rows(FITTED_DR_1, FITTED_DR_2)

FITTED_DR %>%
  ggplot(aes(x=IC50)) +
  geom_histogram(bins=200)


cell.line.lut <- read.csv(paste0(DATADIR, 'annotated/cell_line_master.csv'), row.names = 1)
cell.line.lut$GDSC_COSMIC <- as.numeric(cell.line.lut$GDSC_COSMIC)


FITTED_DR_CTRP <- read.table(paste0(DATADIR, 'annotated/CTRP/v20.data.curves_post_qc.txt'),
                             sep = '\t', header = TRUE) %>%
  select(c(experiment_id, apparent_ec50_umol, area_under_curve, master_cpd_id))
ctrp.exp.meta <- read.table(paste0(DATADIR, 'annotated/CTRP/v20.meta.per_experiment.txt'),
                            sep = '\t', header = TRUE) %>%
  select(c(experiment_id, culture_media, growth_mode, master_ccl_id))

ctrp.cl.meta <- read.table(paste0(DATADIR, 'annotated/CTRP/v20.meta.per_cell_line.txt'),
                           sep = '\t', header = TRUE) %>%
  select(-ccl_availability)

ctrp.compound.meta <- read.table(paste0(DATADIR, 'annotated/CTRP/v20.meta.per_compound.txt'),
                                 sep = '\t', header = TRUE,
                                 quote = '') %>%
  select(c(master_cpd_id, cpd_name, top_test_conc_umol,
           gene_symbol_of_protein_target, target_or_activity_of_compound))
cl.details <- read_csv(paste0(DATADIR,'annotated/Cell_Lines_Details.csv'), col_types = cols())
colnames(cl.details) <- make.names(colnames(cl.details))
cl.details <- cl.details %>%
  mutate(COSMIC.identifier=as.numeric(COSMIC.identifier)) %>%
  mutate_at(vars(Whole.Exome.Sequencing..WES.:Growth.Properties), as.factor) %>%
  # Remove rows with no TCGA label
  filter(!is.na(Cancer.Type..matching.TCGA.label.)) %>%
  mutate(Cancer.Type..matching.TCGA.label.=recode(Cancer.Type..matching.TCGA.label., `COAD/READ`='COREAD'))


ctrp_their_fit <- inner_join(FITTED_DR_CTRP, ctrp.exp.meta,
                             by = 'experiment_id') %>%
  inner_join(ctrp.compound.meta, by = 'master_cpd_id') %>%
  inner_join(ctrp.cl.meta, by = 'master_ccl_id') %>%
  select(experiment_id, apparent_ec50_umol, area_under_curve,
         master_cpd_id, master_ccl_id, cpd_name, ccl_name,
         top_test_conc_umol) %>%
  distinct()



ctrp_our_fit <- ctrp_model_stats %>%
  select(drug, DRUG_ID, cpd_name, maxc,
         experiment_id, CL, culture_media,
         growth_mode, IC50, auc) %>%
  distinct()


comparison <- ctrp_their_fit %>%
  as_tibble() %>%
  full_join(ctrp_our_fit,
            by = c('master_cpd_id'= 'DRUG_ID',
                   'master_ccl_id' = 'CL',
                   'cpd_name', 'experiment_id',
                   'top_test_conc_umol' = 'maxc'))

cor_ctrp <- cor(comparison$area_under_curve,
                comparison$auc, method = 'pearson')

comparison %>%
  ggplot(aes(x=area_under_curve, y=auc)) +
  geom_point(alpha=0.4) +
  scale_x_continuous('CTRP curve fitting') +
  scale_y_continuous('CTRP data with GDSC curve fitting') +
  annotate('text', x = 25, y = 0.5,
           label = round(cor_ctrp, 3)) +
  theme_bw()
