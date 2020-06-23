#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)

# args <- c(HOMEDIR, 'gdsc', paste0(HOMEDIR, 'results/gdsc/'), paste0(HOMEDIR, 'results/gdsc/'))
# args <- c(HOMEDIR, 'ctrp', paste0(HOMEDIR, 'results/ctrp/'), paste0(HOMEDIR, 'results/ctrp/'))


# args should be a vector of length 4-5:
# 1. Home directory of the repository
# 2. "gdsc" or "ctrp"
# 3. UNRES detection result directory
# 4. Output directory
# 5. (optional) Do permutation?
# 6. (Optional) Where is the permuted drug response data stored?


### This script takes identified outlier groups and looks for unique/enriched
### biomarkers in the resistant subcohort.

HOMEDIR <- args[1]
DATADIR <- paste0(HOMEDIR, '/data/combined_for_pipeline/')
INDIR <- args[3]
OUTDIR <- args[4]
do_permutation <- !is.na(args[5]) & args[5] == 'permute'
ANOVAOUTDIR <-  args[6]

suppressMessages(library(tidyverse))

stopifnot(file.exists(paste0(INDIR, 'high_var_results.RData')))


### Data used:
## DRUG_ID, COSMIC_ID, LN_IC50, DRUG_NAME, PUTATIVE_TARGET,
## target_or_activity_of_compound (only CTRP)

# Load useful data and functions
load(paste0(DATADIR, args[2], '_dr_data.RData'))
load(paste0(DATADIR, 'utils.RData'))
load(paste0(DATADIR, 'bems_tidy.RData'))
load(paste0(DATADIR, 'wes_annotation.RData'))
load(paste0(DATADIR, 'cn_annotation.RData'))
load(paste0(DATADIR, 'cancer_gene_lists.RData'))


dataset <- switch(args[2],
                  gdsc = all.gdsc,
                  ctrp = all.ctrp)

all.bems.tidy <- switch(args[2],
                        gdsc = all.bems.tidy.gdsc,
                        ctrp = all.bems.tidy.ctrp)
## Function to classify the mutant cell lines as outliers or not,
## based on the result of the variance change analysis
separate.cell.lines <- function(x){
  simple.bem <- all.bems.tidy %>%
    filter(Tissue == as.character(x['tissue'])) %>%
    select(-Tissue)
  
  # Get Drug response and mutation status data
  selected.dr <- dataset %>%
    filter(DRUG_ID == as.character(x['drug'])) %>%
    right_join(simple.bem, by = 'COSMIC_ID') %>%
    filter(alteration == as.character(x['alteration'])) %>%
    arrange(COSMIC_ID) %>%
    unique() %>%
    mutate(log10.dr = log10(exp(LN_IC50)))
  
  # Get max conc value (if there is a single different value, treat it as and error)
  max.conc <- log10(as.numeric(names(sort(table(selected.dr$MAX_CONC),
                                          decreasing = TRUE)))[1])
  
  queried.data <- selected.dr %>%
    mutate(tissue=as.character(x['tissue']), 
           log10.max.conc=max.conc,
           gene=as.character(x['alteration'])) %>%
    select(tissue, COSMIC_ID, DRUG_ID, PUTATIVE_TARGET,
           gene, mut.status, log10.dr, log10.max.conc) %>%
    mutate_if(is.character, as.factor) %>%
    mutate(COSMIC_ID = as.character(COSMIC_ID))
  
  if (do_permutation){
    perm.dr <- read.csv(paste0(ANOVAOUTDIR, x['tissue'], '_dr.csv'),
                        row.names = 1, check.names = FALSE)
    drug_target <- unique(na.omit(queried.data$PUTATIVE_TARGET))
    drug_target <- ifelse(length(drug_target) == 0, NA,
                          drug_target)
    queried.data <- tibble(LN_IC50_perm=perm.dr[, as.character(x['drug'])],
                           COSMIC_ID=rownames(perm.dr)) %>%
      full_join(queried.data, by = 'COSMIC_ID') %>%
      mutate(log10.dr=log10(exp(LN_IC50_perm)),
             DRUG_ID = as.character(x['drug']),
             PUTATIVE_TARGET = drug_target) %>%
      select(-LN_IC50_perm)
  }
  
  remove.n <- as.numeric(x['n.outliers'])
  
  res <- queried.data %>%
    filter(!is.na(log10.dr)) %>%
    filter(mut.status == 1) %>%
    arrange(-log10.dr) %>%
    rename(cell_line=COSMIC_ID) %>%
    mutate(index=1:n(),
           subgroup=ifelse(index <= remove.n, 'A', 'B')) %>%
    select(cell_line, subgroup)
  return(res)
}

## Function to get CN data at gene level
get.gene.cn.data <- function(assoc, cell.lines){
  region <- strsplit(as.character(assoc['alteration']), ':')[[1]][2]
  region <- strsplit(region, ' ')[[1]][1]
  genes <- strsplit(RACS_to_gene[RACS_to_gene$Identifier == region,]$`Contained genes`,
                    ',')[[1]]
  cn.list <- lapply(genes, function(g) {
    Gene_level_CN[Gene_level_CN$gene == g,
                  map_chr(str_split(as.character(cell.lines$cell_line), '_'), 1)]
  })
  names(cn.list) <- genes
  cn.df <- do.call('rbind', cn.list)
  colnames(cn.df) <- cell.lines$cell_line
  ## Some genes are missing from 'RACS_to_genes'
  return(t(cn.df))
}


## Function to identify CFEs that differ between subgroups from BEMs
get.cfes.from.bem <- function(assoc, cell.lines, differents = TRUE){
  bem <- all.bems.tidy %>%
    filter(Tissue == as.character(assoc['tissue']))
  relevant.bem <- bem %>%
    filter(COSMIC_ID %in% cell.lines$cell_line) %>%
    select(-Tissue) %>%
    pivot_wider(names_from = 'COSMIC_ID', values_from = 'mut.status') %>%
    as.data.frame() %>%
    column_to_rownames('alteration')
  relevant.bem <- relevant.bem[,as.character(cell.lines$cell_line)]
  if (differents){
    which.differ <- apply(relevant.bem, 1,
                          function(x) all(unique(as.numeric(x)[cell.lines$subgroup == 'A']) !=
                                            unique(as.numeric(x)[cell.lines$subgroup == 'B'])) &
                            length(unique(as.numeric(x)[cell.lines$subgroup == 'A'])) == 1 &
                            length(unique(as.numeric(x)[cell.lines$subgroup == 'B'])) == 1
    )
    return(relevant.bem[which.differ,])
  } else {
    return(relevant.bem)
  }
}


## Function that takes an association and returns the specifics about the mutation
## or CNAs present in the corresponding cell lines.
get.specific.info <- function(assoc){
  cell.lines <- separated.cl.list[[rownames(assoc)]]
  alt.type <- ifelse(grepl('_mut', assoc$alteration), 'mut', 'cna')
  
  ### Mutation
  if (alt.type == 'mut'){
    # seq.data <- WES_variants[WES_variants$COSMIC_ID %in% cell.lines$cell_line,]
    seq.data <- map_dfr(cell.lines$cell_line, ~WES_variants %>% 
                          filter(COSMIC_ID == strsplit(as.character(.x), '_')[[1]][1]) %>%
                          mutate(COSMIC_ID = .x))
    
    gene <- strsplit(assoc$alteration, '_mut')[[1]]
    gene.mut <- lapply(cell.lines$cell_line,
                       function(x) seq.data[seq.data$Gene == gene &
                                              seq.data$COSMIC_ID == x, c('SAMPLE', 'COSMIC_ID',
                                                                         'Cancer.Type', 'Gene',
                                                                         'Transcript', 'cDNA',
                                                                         'AA', 'Classification')])
    names(gene.mut) <- cell.lines$cell_line
    return(gene.mut)
  }
  
  ## Copy number
  if (alt.type == 'cna'){
    cna.data <- get.gene.cn.data(assoc, cell.lines)
    if (sum(cell.lines$subgroup == 'A') > 1 & nrow(unique(cna.data)) > 1){
      cn.comparison <- apply(cna.data, 2, function(x) {
        gene <- do.call('rbind', strsplit(x, ','))
        max.cn <- as.numeric(gene[,1])
        if (length(unique(max.cn[cell.lines$subgroup == 'A'])) == 1 &
            length(unique(max.cn[cell.lines$subgroup == 'B'])) == 1){
          return(1)
        }
        if (grepl('gain', assoc$alteration)){
          max.pval <- t.test(max.cn[cell.lines$subgroup == 'A'],
                             max.cn[cell.lines$subgroup == 'B'],
                             alternative = 'greater')$p.value
        } else {
          
          max.pval <- t.test(max.cn[cell.lines$subgroup == 'A'],
                             max.cn[cell.lines$subgroup == 'B'],
                             alternative = 'less')$p.value
        }
        return(max.pval)
      })
      cn.adj <- p.adjust(cn.comparison, method = 'BH')
      return(cn.adj[cn.adj < 0.05])
    } else {
      return(cna.data)
    }
  }
}

rename.cna <- function(x) {
  if (grepl('mutation', x)){
    return(x)
  }
  region <- strsplit(x, ':')[[1]][2]
  region2 <- strsplit(region, ' ')[[1]][1]
  return(unique(paste0(RACS_to_gene$locus[RACS_to_gene$Identifier == region2], ' ',
                       RACS_in_cell_lines$Alteration.Type[RACS_in_cell_lines$Region.identifier == region],
                       ' ', strsplit(region, ' ')[[1]][2])))
}

save.assoc.table <- function(data, file.name, markers = FALSE){
  ### Save table with resulting cases
  rm.cols <- c('association','p.value','corrected.p.value...q.value.','effect',
               'p.value.msi', 'q.value.msi', 'p.value.medium', 'q.value.medium',
               'p.value.growth', 'q.value.growth')
  save.df <- data[order(data$outlier.q.value),
                  !names(data) %in% rm.cols]
  
  save.df$DRUG_ID <- save.df$drug
  
  if (args[2] == 'ctrp'){
    dataset$PUTATIVE_TARGET <- as.character(dataset$PUTATIVE_TARGET)
    dataset$target_or_activity_of_compound <- as.character(dataset$target_or_activity_of_compound)
    dataset$PUTATIVE_TARGET[is.na(dataset$PUTATIVE_TARGET)] <- dataset$target_or_activity_of_compound[is.na(dataset$PUTATIVE_TARGET)]
  }
  
  save.df <- dataset %>%
    filter(DRUG_ID %in% save.df$DRUG_ID) %>%
    select(c(DRUG_ID, DRUG_NAME, PUTATIVE_TARGET)) %>%
    unique() %>%
    right_join(save.df, by = "DRUG_ID") %>%
    mutate_at('alteration', function(x) gsub('_mut', ' mutation', x)) %>%
    mutate_at('alteration', function(x) sapply(x, rename.cna))
  
  
  if (markers){
    save.df <- save.df[,c('tissue', 'DRUG_NAME', 'PUTATIVE_TARGET', 'alteration',
                          'n.samples', 'n.outliers', 'outlier.p.value',
                          'outlier.q.value', 'effect.size','out.cell.lines',
                          'putative.markers')]
    colnames(save.df) <- c('Tissue', 'Drug', 'Drug Target', 'Sensitivity marker',
                           'N. of sensitive samples', 'N. of resistant outliers',
                           'Outlier p-value', 'Outlier FDR', 'Effect size', 
                           'Outlier cell lines', 'Putative resistance markers')
  } else {
    save.df <- save.df[,c('tissue', 'DRUG_NAME', 'PUTATIVE_TARGET', 'alteration',
                          'n.samples', 'n.outliers', 'outlier.p.value',
                          'outlier.q.value', 'effect.size')]
    colnames(save.df) <- c('Tissue', 'Drug', 'Drug Target', 'Sensitivity marker',
                           'N. of sensitive samples', 'N. of resistant outliers',
                           'Outlier p-value', 'Outlier FDR', 'Effect size')
  }
  save.df <- apply(save.df,2,as.character)
  write.table(save.df, file = file.name, row.names = FALSE, sep = '\t',
              quote = FALSE)
}


### Load results data created by 'anova_out_analysis.R'
## Objects: high.var.associations and separated.cl.list
load(paste0(INDIR, 'high_var_results.RData'))
load(paste0(INDIR, 'fdr.tree.value.RData'))
alpha <- 0.15

associations <- high.var.associations
associations <- associations[order(associations$drug),]
associations <- associations[order(associations$alteration),]
associations <- associations[order(associations$tissue),]

save.assoc.table(associations, paste0(OUTDIR, 'all_found_outliers.tsv'))

########################## TABLE 2 (unique associations) ##########################
associations.sign <- associations %>%  filter(outlier.q.value < alpha)
associations.sign$association <- paste0(associations.sign$tissue, ' - ', associations.sign$association)
alpha.n.outliers <- sapply(unique(associations.sign$association), function(x) paste(associations.sign$n.outliers[associations.sign$association == x],
                                                                               collapse = ', '))

### Save table with resulting cases
rm.cols <- c('p.value','corrected.p.value...q.value.','effect',
             'p.value.msi', 'q.value.msi', 'p.value.medium', 'q.value.medium',
             'p.value.growth', 'q.value.growth')
save.df <- associations.sign %>%
  select(-all_of(rm.cols))

save.df$DRUG_ID <- save.df$drug

save.df <- dataset %>%
  filter(DRUG_ID %in% save.df$DRUG_ID) %>%
  select(c(DRUG_ID, DRUG_NAME, PUTATIVE_TARGET)) %>%
  unique() %>%
  right_join(save.df, by = "DRUG_ID") %>%
  mutate_at('alteration', function(x) gsub('_mut', ' mutation', x)) %>%
  mutate_at('alteration', function(x) sapply(x, rename.cna))

save.df <- save.df[,c('association', 'tissue', 'DRUG_NAME', 'PUTATIVE_TARGET', 'alteration',
                      'n.samples')]

save.df <- save.df[!duplicated(save.df$association),]

colnames(save.df) <- c('Assoc', 'Tissue', 'Drug', 'Drug Target', 'Sensitivity marker',
                       'N. of sensitive samples')

save.df$Assoc <- NULL
save.df$`Number of outliers` <- alpha.n.outliers
save.df <- apply(save.df, 2, as.character)
write.table(save.df, file = paste0(OUTDIR, 'unique_outliers_detected.tsv'),
            row.names = FALSE, sep = '\t', quote = FALSE)
##################################

separated.cl.list <- apply(associations, 1, separate.cell.lines)

## Look at the specifics of mutation/CNA
specific.cfe.list <- map(rownames(associations), ~get.specific.info(associations[.x,]))
names(specific.cfe.list) <- rownames(associations)

## Look at co-ocurring CFEs
other.cfe.list <- lapply(rownames(associations), function(x) {
  res <- get.cfes.from.bem(associations[x,],
                           separated.cl.list[[x]])
  if (nrow(res) == 0){
    return(NULL)
  } else {
    return(res)
  }
}
)
names(other.cfe.list) <- rownames(associations)



## Test for overrepresentation and underrepresentation

## Set a filter for 2-3 mutant cell lines minimum
## Perform a power study to find ideally how many mutations we would need
## and the possible best p-values

test.enrichment <- function(ii, adjust = TRUE, method = 'fisher',
                            test = 'enrichment', min.mut.n = 1) {
  cfes <- get.cfes.from.bem(associations[ii,],
                            separated.cl.list[[ii]], differents = FALSE)
  # Only study mutations with at least 2 affected cell lines
  # min.mut.n <- 3
  cfes <- cfes[rowSums(cfes) >= min.mut.n,]
  res <- apply(cfes, 1, function(x) {
    gen <- as.numeric(x)
    both <- data.frame(gen, separated.cl.list[[ii]]$subgroup)
    if (length(unique(both$gen)) == 2){
      tab <- matrix(table(both), nrow = 2, ncol = 2)
    } else {
      tab <- matrix(nrow = 2, ncol = 2)
      rownames(tab) <- c('0','1')
      tab[rownames(table(both)),] <- as.matrix(table(both))
    }
    tab[is.na(tab)] <- 0
    if (method == 'hyper' & test == 'enrichment'){
      ret <- phyper(tab[1,1], sum(tab[1,]), sum(tab[2,]), sum(tab[,1]),
                    lower.tail = TRUE)
    }
    if (method == 'fisher' & test == 'enrichment'){
      ret <- fisher.test(tab, alternative = 'less')$p.value
    }
    return(unname(ret))
  }
  )
  
  if (adjust){
    q.res <- p.adjust(res, method = 'BH')
    return(res[q.res < 0.05])
  } else {
    return(res)
  }
}


enrichment.fisher <- lapply(1:nrow(associations), test.enrichment,
                            adjust = TRUE, method = 'fisher',
                            test = 'enrichment')
names(enrichment.fisher) <- rownames(associations)





## Get rid of cases with no result whatsoever
has.info <- t(sapply(1:nrow(associations), function(x) {
  alt.type <- ifelse(grepl('_mut', associations$alteration[x]), 'mut', 'cna')
  if (alt.type == 'cna'){
    spec.res <- length(specific.cfe.list[[x]]) == 0
  }
  if (alt.type == 'mut'){
    # Look for mutations unique in outliers, or secondary mutations
    mut.table <- specific.cfe.list[[x]]
    n.outliers <- associations$n.outliers[x]
    out.muts <- matrix(unlist(lapply(1:n.outliers,
                                     function(y) sapply(as.character(mut.table[[y]]$AA),
                                                        function(z) c(y, z)))), nrow = 2)
    rest.muts <- matrix(unlist(lapply((n.outliers+1):length(mut.table),
                                      function(y) sapply(as.character(mut.table[[y]]$AA),
                                                         function(z) c(y, z)))), nrow = 2)
    duplicate <- any(unique(out.muts[2,]) %in% unique(rest.muts[2,]))
    no.secondary <- sapply(1:n.outliers, function(r) sum(out.muts[1,] == r) == 1)
    spec.res <- all(duplicate, no.secondary)
  }
  
  rev <- c(spec.res,
           is.null(other.cfe.list[[x]]),
           length(enrichment.fisher[[x]]) == 0)
  return(!rev)
}
))

colnames(has.info) <- c('specific', 'other', 'enrichment')
rownames(has.info) <- rownames(associations)
has.info <- as.data.frame(has.info)
# associations.filtered <- associations[rowSums(has.info) != 0,]


# putative.biomarkers <- vector(mode = 'list', length = nrow(associations.filtered))
# for (ii in 1:nrow(associations.filtered)){
#   # alt.type <- ifelse(grepl('_mut', associations.filtered$alteration[ii]), 'mut', 'cna')
#   # if (alt.type == 'mut'){
#    list(specific.cfe.list[[rownames(associations.filtered)[ii]]],
#       other.cfe.list[[rownames(associations.filtered)[ii]]],
#       enrichment.fisher[[rownames(associations.filtered)[ii]]])[as.logical(which.info)]
#   # }
# }


## Function to print info about each association
show.results <- function(assoc){
  alt.type <- ifelse(grepl('_mut', assoc$alteration), 'mut', 'cna')
  names <- rownames(assoc)
  which.info <- has.info[names,]
  out.cell.lines <- separated.cl.list[[names]] %>%
    filter(subgroup == 'A') %>%
    mutate(cell_line=as.character(cell_line)) %>%
    left_join(dataset, by = c('cell_line' = 'COSMIC_ID')) %>%
    select(GDSC_cell_line) %>% unique() %>% unlist() %>%
    unname() %>% as.character()
  
  # # Header
  # cat('Printing information about:\n')
  # cat('\tTissue: ', assoc$tissue, '\n')
  # cat('\tAlteration: ', assoc$alteration, '\n')
  # cat('\tDrug: ', as.character(unique(dataset[dataset$DRUG_ID == assoc$drug, 'DRUG_NAME'])),
  #     ' (target: ', as.character(unique(dataset[dataset$DRUG_ID == assoc$drug, 'PUTATIVE_TARGET'])),')\n',
  #     sep = '')
  # cat('\tNumber of outliers: ', assoc$n.outliers, '\n')
  # cat('\tOutlier cell lines: ', paste(out.cell.lines, collapse = ', '), '\n\n\n')
  
  spec <-  specific.cfe.list[[names]]
  other <-    other.cfe.list[[names]]
  enrich <- enrichment.fisher[[names]]
  
  putative.biomark <- character()
  
  ## Specific info
  if (alt.type == 'mut'){
    # cat('Mutation types: \n\n')
    # cat('COSMIC_ID','SAMPLE','Gene',
    #     'cDNA','AA', 'Classification',
    #     sep = '\t')
    # print.spec <- lapply(names(spec), function(x) {
    #   print(unname(as.data.frame(spec[[x]][c('COSMIC_ID','SAMPLE','Gene',
    #                                          'cDNA','AA', 'Classification')])
    #   ), row.names = FALSE, sep = '\t')
    # }
    # )
    if (which.info$specific){
      n.outliers <- assoc$n.outliers
      out.muts <- matrix(unlist(lapply(1:n.outliers,
                                       function(y) sapply(as.character(spec[[y]]$AA),
                                                          function(z) c(y, z)))), nrow = 2)
      rest.muts <- matrix(unlist(lapply((n.outliers+1):length(spec),
                                        function(y) sapply(as.character(spec[[y]]$AA),
                                                           function(z) c(y, z)))), nrow = 2)
      putative.biomark <- paste(assoc$alteration, out.muts[2, !out.muts[2,] %in% rest.muts[2,]],
                                sep = '-')
      # cat('\nPutative resistance marker:\t', putative.biomark)
      
    }
  }
  
  
  if (alt.type == 'cna'){
    if (which.info$specific){
      # cat('Copy number details: \n')
      if (assoc$n.outliers == 1){
        # cat('\tVisual analysis of copy numbers is required, as there is a single outlier\n')
      } else {
        # cat('\tGenes with a significant change in copy number: \n')
        # cat('\t\t', paste(names(spec), collapse = ', '))
        # Append driver genes with significant copy number change to the table
        driver.cn.genes <- names(spec)[names(spec) %in% cancer.genes.gdsc$Gene]
        if (length(driver.cn.genes > 0)){
          putative.biomark <- paste0('CNA: ', driver.cn.genes,
                                     collapse = ', ')
          
          # cat('\nPutative resistance marker:\t', putative.biomark)
        }
      }
      # } else {
      #   cat('No significant copy number difference between outliers and the rest')
      #   
    }
  }
  
  ## Other BEM info
  # cat('\n\n')
  if (is.null(other)){
    # cat('No unique markers (BEM) corresponding to outliers were found')
  } else {
    # cat('Unique markers for outliers: \n')
    for (i in 1:nrow(other)){
      present <- ifelse(other[i,1] == 1, 'Present in outliers', 'Absent in outliers')
      # cat('\t', rownames(other)[i], '\t(', present, ')\n', sep = '')
      # if (other[i,1] == 1){
        cfes.case <- rownames(other)[i]
        cfes.is.absent <-  other[i,1] == 0
        # Fix HypMet text
        if (grepl('HypMET',cfes.case)){
          hypmet.case <- strsplit(cfes.case[grepl('HypMET',cfes.case)], '[()]')[[1]][2]
          if (hypmet.case %in% cancer.genes.gdsc$Gene){
            putative.biomark <- c(putative.biomark, paste0(ifelse(cfes.is.absent,
                                                                  'Lack of ', ''),
                                                           hypmet.case, '_HypMET'))
          }
        }
        
        
        # Fix CNA text
        if (grepl('cna', cfes.case)){
          cna.case <- sapply(cfes.case,
                             function(x) {
                               region <- strsplit(x, ':')[[1]][2]
                               region2 <- strsplit(region, ' ')[[1]][1]
                               res.str <- unique(paste0(RACS_to_gene$locus[RACS_to_gene$Identifier == region2], ' ',
                                                        RACS_in_cell_lines$Alteration.Type[RACS_in_cell_lines$Region.identifier == region],
                                                        ' ', strsplit(region, ' ')[[1]][2]))
                               # Erase the result if there are no driver genes
                               if (endsWith(res.str, ' NA')){
                                 res.str <- ''
                               }
                               return(res.str)
                             }
          )
          if (nchar(cna.case) > 0){
            putative.biomark <- c(putative.biomark, paste0(ifelse(cfes.is.absent,
                                                                  'Lack of ', ''), cna.case))
          }
        }
        
        # Fix mutations text
        if (grepl('mut', cfes.case)){
          mut.case <- strsplit(cfes.case, '_mut')[[1]][1]
          if (mut.case %in% cancer.genes.gdsc$Gene){
            putative.biomark <- c(putative.biomark, paste0(ifelse(cfes.is.absent,
                                                                  'Lack of ', ''),
                                                           mut.case, ' Mut'))
          }
        }
      # }
    }
  }
  
  
  ## Enrichment
  # cat('\n\n')
  if (length(enrich) == 0){
    # cat('No enriched markers in outliers. \n\n')
  } else {
    # cat('Enriched markers in outliers: \n')
    # cat('\t', paste(names(enrich), collapse = ', '))
    # cat('\n\n')
    # 
    for (ii in 1:length(enrich)){
      cfes.case <- names(enrich)[ii]
      # Fix HypMet text
      if (grepl('HypMET',cfes.case)){
        hypmet.case <- strsplit(cfes.case[grepl('HypMET',cfes.case)], '[()]')[[1]][2]
        if (hypmet.case %in% cancer.genes.gdsc$Gene){
          putative.biomark <- c(putative.biomark, paste0(hypmet.case, '_HypMET'))
        }
      }
      
      
      # Fix CNA text
      if (grepl('cna', cfes.case)){
        cna.case <- sapply(cfes.case,
                           function(x) {
                             region <- strsplit(x, ':')[[1]][2]
                             region2 <- strsplit(region, ' ')[[1]][1]
                             res.str <- unique(paste0(RACS_to_gene$locus[RACS_to_gene$Identifier == region2], ' ',
                                                      RACS_in_cell_lines$Alteration.Type[RACS_in_cell_lines$Region.identifier == region],
                                                      ' ', strsplit(region, ' ')[[1]][2]))
                             # Erase the result if there are no driver genes
                             if (endsWith(res.str, ' NA')){
                               res.str <- ''
                             }
                             return(res.str)
                           }
        )
        if (nchar(cna.case) > 0){
          putative.biomark <- c(putative.biomark, cna.case)
        }
      }
      
      # Fix mutations text
      if (grepl('mut', cfes.case)){
        mut.case <- strsplit(cfes.case, '_mut')[[1]][1]
        if (mut.case %in% cancer.genes.gdsc$Gene){
          putative.biomark <- c(putative.biomark, paste0(mut.case, ' Mut'))
        }
      }
    }
    
  }
  return(putative.biomark)
}


all.putative.biomarkers <- vector(mode = 'list', length = nrow(associations))
names(all.putative.biomarkers) <- rownames(associations)
for (ii in 1:nrow(associations)){
  # cat('\nCase #', ii, ':\n')
  all.putative.biomarkers[[ii]] <- show.results(associations[ii,])
  # cat(rep('-', 80),'\n', sep = '')
}

# Remove empty results
associations$putative.markers <- sapply(all.putative.biomarkers, function(x) paste(x, collapse = ', '))
# associations <- associations[sapply(associations$putative.markers, nchar) != 0,]
associations$out.cell.lines <- sapply(rownames(associations),
                                      function(x){
                                        separated.cl.list[[x]] %>%
                                          filter(subgroup == 'A') %>%
                                          mutate(cell_line=as.character(cell_line)) %>%
                                          left_join(dataset, by = c('cell_line' = 'COSMIC_ID')) %>%
                                          select(GDSC_cell_line) %>% unique() %>% unlist() %>%
                                          unname() %>% as.character() %>%
                                          paste(collapse = ', ')
                                      })



save(associations, file = paste0(OUTDIR, 'putative_res_markers.RData'))

associations.filtered <- associations[with(associations, outlier.q.value <= alpha),]

save.assoc.table(associations.filtered,
                 paste0(OUTDIR, 'putative_res_markers.tsv'),
                 markers = TRUE)

