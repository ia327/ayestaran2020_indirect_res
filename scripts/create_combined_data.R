#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)

### This script creates an .RData object (combined_data.Rdata)
### that contains data and functions used by many of the different scripts.

library(readxl)
suppressMessages(library(tidyverse))
HOMEDIR <- args[1]
DATADIR <- paste0(HOMEDIR,'/data/')
system(paste0('mkdir -p ', DATADIR, 'combined_for_pipeline'))


# Information about screened drugs
# drug_info_old <- read_excel(paste0(DATADIR,'annotated/Screened_Compounds.xlsx'))
drug_info <- read_csv(paste0(DATADIR,'annotated/screened_compounds_rel_8.0.csv'), 
                      col_types = cols())

# Lookup for comparing GDSC and CTRP data
cell.line.lut <- read.csv(paste0(DATADIR, 'annotated/cell_line_master.csv'), row.names = 1)
cell.line.lut$GDSC_COSMIC <- as.character(cell.line.lut$GDSC_COSMIC)


##### GDSC molecular characterisation
# Sequencing data (whole-exome)
WES_variants <- read_csv(paste0(DATADIR, "sequencing/WES_variants.csv"),
                         col_types = cols(.default = "c"))
colnames(WES_variants) <- make.names(colnames(WES_variants))

# Might create a warning regarding values at column 9, but these are not
# relevant to our analysis

# Copy number data
RACS_in_cell_lines <- read_csv(paste0(DATADIR, "copy_number/RACS_in_cell_lines.csv"),
                               col_types = cols())
colnames(RACS_in_cell_lines) <- make.names(colnames(RACS_in_cell_lines))


Gene_level_CN <- suppressWarnings(read_csv(paste0(DATADIR,'copy_number/Gene_level_CN.csv.zip'),
                                           skip = 1, col_types = cols(.default = "c")))
names(Gene_level_CN)[1:4] <- c('gene', 'chr', 'start', 'stop')

# Gene_level_CN_tidy <- Gene_level_CN %>%
#   gather(COSMIC_ID, CN, -gene, -chr, -start, -stop) %>%
#   mutate(CN=map(CN, str_split, ','))
# Gene_level_CN <- read_csv(paste0(DATADIR,'copy_number/gene_level_cnv_latest.csv.gz'))


# Lookup table for RACS and the genes they contain
RACS_to_gene <- suppressWarnings(read_csv(paste0(DATADIR, 'copy_number/RACS_to_gene.csv'),
                                          col_types = cols()))

# BEMs (binary event matrix) and which tissues
all.bems <-  list.files(paste0(DATADIR,'BEMs/'), pattern = 'MOBEM') %>% 
  grep('PANCAN', ., value = TRUE, invert = TRUE)
all.tissues <- unname(sapply(all.bems,
                             function(x) unlist(strsplit(x, split = '_'))[1]))

tissue.cols <- c('#E6194B', '#3cb44b', '#ffe119', '#0082c8', '#f58231',
                 '#911eb4', '#46f0f0', '#f032e6', '#d2f53c', '#fabebe',
                 '#008080', '#C972FF', '#aa6e28', '#CCC263', '#800000',
                 '#aaffc3', '#808000', '#ffd8b1', '#000080', '#000000',
                 '#C8C8C8')
names(tissue.cols) <- c(all.tissues, 'Not significant')

all.bems.tidy <- lapply(seq(length(all.bems)), function(x) {
  bem <- read.table(paste0(DATADIR,'BEMs/',all.bems[x]),
                    sep = '\t', row.names = 1, header = TRUE,
                    check.names = FALSE)
  bem %>% 
    rownames_to_column('alteration') %>%
    gather(COSMIC_ID, mut.status, -alteration) %>%
    mutate(Tissue=all.tissues[x],
           COSMIC_ID=as.character(COSMIC_ID)) %>%
    return()
})
all.bems.tidy.gdsc <- do.call('rbind', all.bems.tidy)

##### GDSC
# Drug response data
FITTED_DR_1 <- read_csv(paste0(DATADIR, "drug/GDSC1_fitted_dose_response_17Jul19.csv"), col_types = cols())
FITTED_DR_2 <- read_csv(paste0(DATADIR, "drug/GDSC2_fitted_dose_response_17Jul19.csv"), col_types = cols())
FITTED_DR <- bind_rows(FITTED_DR_1, FITTED_DR_2)
FITTED_DR$COSMIC_ID <- as.character(FITTED_DR$COSMIC_ID)

# Information about each cell line
cl.details <- read_csv(paste0(DATADIR,'annotated/Cell_Lines_Details.csv'), col_types = cols())
colnames(cl.details) <- make.names(colnames(cl.details))
cl.details <- cl.details %>%
  mutate(COSMIC.identifier=as.character(COSMIC.identifier)) %>%
  mutate_at(vars(Whole.Exome.Sequencing..WES.:Growth.Properties), as.factor) %>%
  # Remove rows with no TCGA label
  filter(!is.na(Cancer.Type..matching.TCGA.label.)) %>%
  mutate(Cancer.Type..matching.TCGA.label.=recode(Cancer.Type..matching.TCGA.label., `COAD/READ`='COREAD'))


# Combine into one data frame
all.gdsc <- inner_join(FITTED_DR, cl.details,
                       by = c('COSMIC_ID' = 'COSMIC.identifier')) %>%
  inner_join(drug_info, by = c('DRUG_ID', 'DRUG_NAME', 'PUTATIVE_TARGET', 'PATHWAY_NAME')) %>%
  left_join(cell.line.lut, by = c('COSMIC_ID' = 'GDSC_COSMIC',
                                  'GDSC.Tissue.descriptor.1' = 'GDSC_tissue')) %>%
  select(-Sample.Name, -CELL_LINE_NAME, -Cancer.Type..matching.TCGA.label.) %>%
  mutate(DRUG_ID=paste0(DRUG_ID, '.', str_sub(DATASET, -1)))


# Cancer gene census (COSMIC)
cancer.gene.census <- read_csv(file = paste0(DATADIR, 'annotated/cancer_gene_census.csv'), 
                               col_types = cols())

# Cancer genes (from Iorio et al, 2016)
cancer.genes.gdsc <- read_csv(paste0(DATADIR, 'annotated/cancer_genes_gdsc.csv'),
                              col_types=cols())
cancer.genes.gdsc$`Cancer Type`[cancer.genes.gdsc$`Cancer Type` == 'COAD/READ'] <- 'COREAD'


##### CTRP
# ctrp.exp.meta <- read.table(paste0(DATADIR, 'annotated/CTRP/v20.meta.per_experiment.txt'),
#                             sep = '\t', header = TRUE) %>%
#   select(c(experiment_id, culture_media, growth_mode, master_ccl_id))
#
# ctrp.cl.meta <- read.table(paste0(DATADIR, 'annotated/CTRP/v20.meta.per_cell_line.txt'),
#                            sep = '\t', header = TRUE) %>%
#   select(-ccl_availability)
#
# ctrp.compound.meta <- read.table(paste0(DATADIR, 'annotated/CTRP/v20.meta.per_compound.txt'),
#                                  sep = '\t', header = TRUE,
#                                  quote = '') %>%
#   select(c(master_cpd_id, cpd_name, top_test_conc_umol,
#            gene_symbol_of_protein_target, target_or_activity_of_compound))
#
# FITTED_DR_CTRP <- read.table(paste0(DATADIR, 'annotated/CTRP/v20.data.curves_post_qc.txt'),
#                              sep = '\t', header = TRUE) %>%
#   select(c(experiment_id, apparent_ec50_umol, area_under_curve, master_cpd_id))
#
# all.ctrp <- inner_join(FITTED_DR_CTRP, ctrp.exp.meta,
#                        by = 'experiment_id') %>%
#   inner_join(ctrp.compound.meta, by = 'master_cpd_id') %>%
#   inner_join(ctrp.cl.meta, by = 'master_ccl_id') %>%
#   inner_join(cell.line.lut, by = c('ccl_name' = 'CTD2_cell_line')) %>%
#   filter(!is.na(GDSC_cell_line)) %>%
#   inner_join(cl.details, by = c('GDSC_COSMIC' = 'COSMIC.identifier')) %>%
#   mutate_at('apparent_ec50_umol', log) %>%
#   mutate(master_cpd_id=as.character(master_cpd_id),
#          gene_symbol_of_protein_target=as.character(gene_symbol_of_protein_target),
#          target_or_activity_of_compound=as.character(target_or_activity_of_compound),
#          gene_symbol_of_protein_target=ifelse(is.na(gene_symbol_of_protein_target),
#                                               target_or_activity_of_compound,
#                                               gene_symbol_of_protein_target)) %>%
#   rename(PUTATIVE_TARGET = gene_symbol_of_protein_target,
#          DRUG_NAME = cpd_name,
#          DRUG_ID = master_cpd_id,
#          COSMIC_ID = GDSC_COSMIC,
#          TCGA_DESC = Cancer.Type..matching.TCGA.label.,
#          # LN_IC50 = apparent_ec50_umol,
#          # AUC instead IC_50 for CTRP
#          LN_IC50 = area_under_curve,
#          MAX_CONC = top_test_conc_umol)


## Use refitted data
FITTED_DR_CTRP <- read_csv(paste0(DATADIR, 'drug/CTRP_refitted_dose_response.csv'),
                           col_types = cols())

ctrp.cl.meta <- read.table(paste0(DATADIR, 'annotated/CTRP/v20.meta.per_cell_line.txt'),
                           sep = '\t', header = TRUE) %>%
  select(-ccl_availability)

all.ctrp <- FITTED_DR_CTRP  %>%
  inner_join(ctrp.cl.meta, by = c('CL' = 'master_ccl_id')) %>%
  inner_join(cell.line.lut, by = c('ccl_name' = 'CTD2_cell_line')) %>%
  filter(!is.na(GDSC_cell_line)) %>%
  inner_join(cl.details, by = c('GDSC_COSMIC' = 'COSMIC.identifier')) %>%
  select(-Screen.Medium, -Growth.Properties) %>% # Use CTRP media and growth info
  mutate(DRUG_ID=as.character(DRUG_ID),
         gene_symbol_of_protein_target=as.character(gene_symbol_of_protein_target),
         target_or_activity_of_compound=as.character(target_or_activity_of_compound),
         gene_symbol_of_protein_target=ifelse(is.na(gene_symbol_of_protein_target), 
                                              target_or_activity_of_compound,
                                              gene_symbol_of_protein_target),
         COSMIC_ID=paste0(GDSC_COSMIC,'_', experiment_id)) %>%
  rename(PUTATIVE_TARGET = gene_symbol_of_protein_target,
         DRUG_NAME = cpd_name,
         TCGA_DESC = Cancer.Type..matching.TCGA.label.,
         LN_IC50 = IC50,
         MAX_CONC = maxc,
         Screen.Medium = culture_media,
         Growth.Properties = growth_mode)

### CTRP COSMIC IDs have appended the corresponding experiment ID, as there are
## repeated experiments on the same cell line, sometimes with different screen media
all.bems.tidy.ctrp <- all.ctrp %>%
  select(GDSC_COSMIC, COSMIC_ID) %>%
  distinct() %>%
  inner_join(all.bems.tidy.gdsc, by = c('GDSC_COSMIC' = 'COSMIC_ID')) %>%
  select(-GDSC_COSMIC)

######################### FUNCTIONS #########################

# Function to obtain IC50 values to compare between
# wt and mut for a determined gene-drug association
wt.vs.mut <- function(tissue, drug.id, gene, ref.data,
                      plot.beeswarm = TRUE, colour = NA){
  # Get BEM
  simple.bem <- all.bems.tidy %>%
    filter(Tissue == tissue) %>%
    select(-Tissue)
  
  # Get Drug response and mutation status data
  selected.dr <- ref.data %>%
    filter(DRUG_ID == drug.id) %>%
    right_join(simple.bem, by = 'COSMIC_ID') %>%
    filter(alteration == gene & !is.na(LN_IC50)) %>%
    arrange(COSMIC_ID) %>%
    unique() %>%
    mutate(log10.dr = log10(exp(LN_IC50)))
  
  # Get max conc value (if there is a single different value, treat it as and error)
  max.conc <- log10(as.numeric(names(sort(table(selected.dr$MAX_CONC),
                                          decreasing = TRUE)))[1])
  if (plot.beeswarm){
    gene.label <- ifelse(grepl('mut', gene),
                         strsplit(gene, '_mut')[[1]][1], gene)
    if (is.na(colour)){
      colour <- tissue.cols[tissue]
    }
    
    
    p <- selected.dr %>%
      ggplot(aes(x=mut.status, y=log10.dr, group=mut.status)) +
      geom_boxplot(col=colour, fill=colour, alpha=0.6) +
      geom_beeswarm(col=colour, fill=colour) +
      geom_hline(yintercept = max.conc, linetype = 2) +
      theme_minimal() +
      scale_x_continuous(gene.label, breaks = c(0,1), labels = c('WT', 'Mut')) +
      scale_y_continuous('log10(IC50)') +
      ggtitle(paste0(tissue, ' - ', unique(selected.dr$DRUG_NAME),
                     ' (', unique(selected.dr$PUTATIVE_TARGET), ')'))
    print(p)
  }
  
  to.return <- selected.dr %>%
    mutate(tissue=tissue, log10.max.conc=max.conc) %>%
    rename(gene=alteration) %>%
    select(tissue, COSMIC_ID, DRUG_ID, DRUG_ID, PUTATIVE_TARGET,
           gene, mut.status, log10.dr, log10.max.conc) %>%
    mutate_if(is.character, as.factor)
  
  return(to.return)
}


save('all.gdsc', file = paste0(DATADIR, 'combined_for_pipeline/gdsc_dr_data.RData'))
save('all.ctrp', file = paste0(DATADIR, 'combined_for_pipeline/ctrp_dr_data.RData'))
save(list= c('all.tissues', 'all.bems', 'all.bems.tidy.gdsc', 'all.bems.tidy.ctrp'),
     file = paste0(DATADIR, 'combined_for_pipeline/bems_tidy.RData'))
save('WES_variants', file = paste0(DATADIR, 'combined_for_pipeline/wes_annotation.RData'))
save(list = c('RACS_to_gene', 'Gene_level_CN', 'RACS_in_cell_lines'),
     file = paste0(DATADIR, 'combined_for_pipeline/cn_annotation.RData'))
save(list = c('cancer.gene.census', 'cancer.genes.gdsc'),
     file = paste0(DATADIR, 'combined_for_pipeline/cancer_gene_lists.RData'))
save(list = c('wt.vs.mut', 'tissue.cols'),
     file = paste0(DATADIR, 'combined_for_pipeline/utils.RData'))


