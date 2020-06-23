# devtools::install_github("cancerrxgene/gdscIC50", build_vignettes=FALSE)
library(gdscIC50)
library(tidyverse)

HOMEDIR <- '../../'
# load(paste0(HOMEDIR, 'data/combined_data.RData'))

## GDSC1
gdsc1_raw <- read_csv(paste0(HOMEDIR, 'data/drug/GDSC1_public_raw_data_17Jul19.csv.zip'),
                      col_types=cols(
                        RESEARCH_PROJECT = col_character(),
                        BARCODE = col_character(),
                        SCAN_ID = col_double(),
                        DATE_CREATED = col_datetime(format = ""),
                        SCAN_DATE = col_logical(),
                        CELL_ID = col_double(),
                        MASTER_CELL_ID = col_double(),
                        COSMIC_ID = col_double(),
                        CELL_LINE_NAME = col_character(),
                        SEEDING_DENSITY = col_logical(),
                        DRUGSET_ID = col_character(),
                        ASSAY = col_character(),
                        DURATION = col_double(),
                        POSITION = col_double(),
                        TAG = col_character(),
                        DRUG_ID = col_double(),
                        CONC = col_double(),
                        INTENSITY = col_double()
                      ))

gdsc1_raw <- removeFailedDrugs(gdsc1_raw)
gdsc1_raw <- removeMissingDrugs(gdsc1_raw)

normalized_gdsc1 <- normalizeData(gdsc1_raw,
                                  trim = T,
                                  neg_control = "NC-0",
                                  pos_control = "B")

scaled_gdsc1 <- setConcsForNlme(normalized_gdsc1, group_conc_ranges = F)
nlme_data_1 <- prepNlmeData(scaled_gdsc1,
                            cl_id = "COSMIC_ID",
                            drug_specifiers = c("DRUG_ID_lib", "maxc"))
nlme_model_1 <- fitModelNlmeData(nlme_data_1, isLargeData = F)
gdsc_model_stats_1 <- calcNlmeStats(nlme_model_1, nlme_data_1)



## GDSC2
gdsc2_raw <- read_csv(paste0(HOMEDIR, 'data/drug/GDSC2_public_raw_data_17Jul19.csv.zip'),
                      col_types=cols(
                        RESEARCH_PROJECT = col_character(),
                        BARCODE = col_character(),
                        SCAN_ID = col_double(),
                        DATE_CREATED = col_datetime(format = ""),
                        SCAN_DATE = col_character(),
                        CELL_ID = col_double(),
                        MASTER_CELL_ID = col_double(),
                        COSMIC_ID = col_double(),
                        CELL_LINE_NAME = col_character(),
                        SEEDING_DENSITY = col_double(),
                        DRUGSET_ID = col_character(),
                        ASSAY = col_character(),
                        DURATION = col_double(),
                        POSITION = col_double(),
                        TAG = col_character(),
                        DRUG_ID = col_double(),
                        CONC = col_double(),
                        INTENSITY = col_double()
                      ))

gdsc2_raw <- removeFailedDrugs(gdsc2_raw)
gdsc2_raw <- removeMissingDrugs(gdsc2_raw)

normalized_gdsc2 <- normalizeData(gdsc2_raw,
                                  trim = T,
                                  neg_control = "NC-0",
                                  pos_control = "B")

scaled_gdsc2 <- setConcsForNlme(normalized_gdsc2, group_conc_ranges = F)
nlme_data_2 <- prepNlmeData(scaled_gdsc2,
                          cl_id = "COSMIC_ID",
                          drug_specifiers = c("DRUG_ID_lib", "maxc"))
nlme_model_2 <- fitModelNlmeData(nlme_data_2, isLargeData = F)
gdsc_model_stats_2 <- calcNlmeStats(nlme_model_2, nlme_data_2)


## CTRP
ctrp_raw <- read_tsv(paste0(HOMEDIR, 'data/annotated/CTRP/v20.data.per_cpd_post_qc.txt.zip'))

ctrp_compound_meta <- read_tsv(paste0(HOMEDIR, 'data/annotated/CTRP/v20.meta.per_compound.txt')) %>%
  select(c(master_cpd_id, cpd_name, top_test_conc_umol,
           gene_symbol_of_protein_target, target_or_activity_of_compound))

# ctrp_raw <- ctrp_compound_meta %>%
#   filter(cpd_name=='gefitinib') %>%
#   left_join(ctrp_raw, by = "master_cpd_id") %>%
#   distinct()


ctrp_exp_meta <- read_tsv(paste0(HOMEDIR, 'data/annotated/CTRP/v20.meta.per_experiment.txt')) %>%
  select(c(experiment_id, culture_media, growth_mode, master_ccl_id)) %>%
  distinct()


ctrp_model_data <- ctrp_raw %>%
  left_join(ctrp_exp_meta, by = 'experiment_id') %>%
  left_join(ctrp_compound_meta, by = "master_cpd_id") %>%
  rename(CONC=cpd_conc_umol,
         maxc=top_test_conc_umol,
         CL = master_ccl_id,
         DRUG_ID=master_cpd_id) %>%
  mutate(x = (log(CONC/maxc)/log(2)) +
           9,
         y = 1-cpd_avg_pv) %>%
  mutate(drug=paste0(DRUG_ID, '_', experiment_id),
         drug_spec="DRUG_ID+experiment_id",
         DRUG_ID = as.numeric(DRUG_ID),
         DRUG_ID_lib = DRUG_ID)

ctrp_model <- fitModelNlmeData(ctrp_model_data, isLargeData = F)
ctrp_model_stats <- calcNlmeStats(ctrp_model, ctrp_model_data)


ctrp_dr <- ctrp_model_stats %>%
  select(DRUG_ID, cpd_name, maxc,
         experiment_id, CL, culture_media,
         growth_mode, IC50, auc, RMSE,
         gene_symbol_of_protein_target,
         target_or_activity_of_compound) %>%
  distinct()

save(list = c('gdsc_model_stats_1', 'gdsc_model_stats_2', 'ctrp_model_stats'),
     file = paste0(HOMEDIR, 'data/drug/curve_fitting_stats.RData'))

save(ctrp_model_stats,
     file = paste0(HOMEDIR, 'data/drug/ctrp_all_fitted.RData'))

write_csv(ctrp_dr,path = paste0(HOMEDIR, 'data/drug/CTRP_refitted_dose_response.csv'))
