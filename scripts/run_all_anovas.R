#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)

## Arguments:
# 1. Home directory
# 2. Dataset
# 3. Output prefix
# 4. (optional) Do permutation?


### Data used:
## TCGA_DESC, Drug.Response,
## Microsatellite..instability.Status..MSI.,
## Screen.Medium, Growth.Properties, COSMIC_ID
## COSMIC_ID, DRUG_ID, LN_IC50




# Run ANOVA for all cancer types.
library(methods)
suppressMessages(library(tidyverse))
suppressMessages(library(DRANOVA))
library(parallel)

HOMEDIR <- args[1]
DATADIR <- paste0(HOMEDIR, '/data/combined_for_pipeline/')
OUTDIR <- args[3]
# Optional argument to permute IC50s within a tissue
do_permutation <- !is.na(args[4]) & args[4] == 'permute'


load(paste0(DATADIR, args[2], '_dr_data.RData'))
load(paste0(DATADIR, 'bems_tidy.RData'))

dataset <- switch(args[2],
                  gdsc = all.gdsc,
                  ctrp = all.ctrp)
all.bems.tidy <- switch(args[2],
                        gdsc = all.bems.tidy.gdsc,
                        ctrp = all.bems.tidy.ctrp)


# Function for running the tissue specific ANOVA
run.tissue.ANOVA <- function(tissue, ref.data){
  cat(tissue, '\n')

  tissue.details <- ref.data %>%
    filter(TCGA_DESC == tissue) %>%
    select(c(TCGA_DESC, Drug.Response,
             Microsatellite..instability.Status..MSI.,
             Screen.Medium, Growth.Properties, COSMIC_ID)) %>%
    unique() %>%
    mutate_if(is.factor, as.character)

  # What columns correspond to covariates?
  tissue.covariates <- tissue.details %>%
    mutate(COSMIC_ID=as.character(COSMIC_ID)) %>%
    arrange(COSMIC_ID) %>%
    column_to_rownames('COSMIC_ID') %>%
    select(Microsatellite..instability.Status..MSI.,
           Screen.Medium, Growth.Properties)

  names(tissue.covariates) <- c('msi', 'medium', 'growth')
  tissue.covariates[is.na(tissue.covariates)] <- 'NA'

  # Select drug response data corresponding to the tissue
  res.data <- ref.data %>%
    select(c(COSMIC_ID, DRUG_ID, LN_IC50)) %>%
    mutate(COSMIC_ID=as.character(COSMIC_ID)) %>%
    filter(COSMIC_ID %in% rownames(tissue.covariates)) %>%
    unique() %>%
    spread(DRUG_ID, LN_IC50) %>%
    column_to_rownames('COSMIC_ID') %>%
    as.data.frame() %>%
    as.matrix()

  # Permute the data randomly if stated
  if (do_permutation){
    perm.res.data <- apply(res.data, 2,sample)
    rownames(perm.res.data) <- rownames(res.data)
    write.csv(perm.res.data, file = paste0(OUTDIR, tissue,'_dr.csv'))
    res.data <- perm.res.data
  }

  # Select BEM data:
  # Remove Hyper methylation data, and get overlap in case there is
  # some extra cell line in the BEM (as in COREAD, where a cell line
  # has no drug response data)
  bem.data <- all.bems.tidy %>%
    filter(Tissue == tissue) %>%
    spread(alteration, mut.status) %>%
    dplyr::select(-Tissue, -ends_with('_HypMET')) %>%
    as.data.frame() %>%
    mutate(COSMIC_ID=as.character(COSMIC_ID)) %>%
    filter(COSMIC_ID %in%rownames(res.data)) %>%
    arrange(COSMIC_ID) %>%
    column_to_rownames('COSMIC_ID') %>%
    as.matrix()


  # Create DRANOVA object and run analysis with the covariates
  dranova <- new('DRANOVA', res = res.data,
                 bem = bem.data,
                 covariates = tissue.covariates)

  dranova.out <- runANOVA(dranova, cov.names = colnames(tissue.covariates))

  # Save it as a data frame
  saveDRANOVAout(dranova.out, file = paste0(OUTDIR, tissue, '_anova_out.csv'),
                 overwrite = TRUE)
  return('OK')
}

tissues.to.anova <- dataset %>%
  filter(Drug.Response == 'Y') %>%
  select(TCGA_DESC, COSMIC_ID) %>%
  distinct() %>%
  count(TCGA_DESC) %>%
  filter(TCGA_DESC %in% unique(all.bems.tidy$Tissue)) %>%
  filter(n >= 8)

foo <- mclapply(tissues.to.anova$TCGA_DESC, function(tiss) {
  run.tissue.ANOVA(tiss, ref.data = dataset)
  }, mc.cores = 10)
