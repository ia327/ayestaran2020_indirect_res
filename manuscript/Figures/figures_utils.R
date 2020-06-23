library(gdscIC50)


#### Volcano plots
# function for transforming y-axis
reverselog_trans <- function (base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^-x
  trans_new(paste0("log-", format(base)), trans, inv, log_breaks(base = base),
            domain = c(1e-100, Inf))
}

# function for relabeling the y-axis
scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scientific_format()(x)))
}

rename.cna <- function(x) {
  region <- map_chr(str_split(x, ':'), 2)
  region2 <- map_chr(str_split(region, ' '), 1)
  if (grepl(' ', region)){
    genes <-  map_chr(str_split(region, ' '), 2)
  } else {
    genes <- NA
  }
  chr.reg <- RACS_to_gene$locus[RACS_to_gene$Identifier == region2]
  cn.type <- unique(RACS_in_cell_lines$Alteration.Type[RACS_in_cell_lines$Region.identifier == region])
  cn.type <- str_sub(cn.type, 1, 3)
  if (is.na(genes)){
    return(paste0(chr.reg, cn.type))
  } else {
    return(paste0(chr.reg, genes, cn.type))
  }
}


## Function for curve fitting in different log base
l3_model2_base10 <- function (lx, maxc, xmid, scal) 
{
  x <- getXfromConc(10^lx, maxc)
  yhat <- 1 - logist3(x, xmid, scal)
  return(yhat)
}

## Function to plot curve fittings
plotResponse_custom <- function (model_stats, cell_line, drug_identifier) {
  all.plot_data <- model_stats %>% 
    filter(CL %in% cell_line$COSMIC_ID) %>%
    filter(drug %in% drug_identifier) %>% 
    mutate(lx = log10(getConcFromX(x, maxc)),
           lxmid = log10(getConcFromX(xmid, maxc))) %>%
    left_join(cell_line, by = c('CL' = 'COSMIC_ID')) %>%
    distinct()
  plot_low_x_all <- numeric(length = nrow(cell_line))
  plot_high_x_all <- numeric(length = nrow(cell_line))
  counter <- 0
  for (ii in seq(nrow(cell_line))){
    plot_data <- all.plot_data %>% 
      filter(CL == cell_line$COSMIC_ID[ii] &
             drug %in% drug_identifier)
    if(nrow(plot_data) == 0){
      next
    }
    counter <- counter + 1
    IC50 <- unique(plot_data$IC50)
    stopifnot(length(IC50) == 1)
    auc <- unique(plot_data$auc)
    stopifnot(length(auc) == 1)
    rmse <- unique(plot_data$RMSE)
    stopifnot(length(rmse) == 1)
    drug_id <- unique(plot_data$DRUG_ID_lib)
    stopifnot(length(drug_id) == 1)
    cell_line_name <- unique(plot_data$CELL_LINE_NAME)
    cell_line_name <- ifelse(is.null(cell_line_name), "", cell_line_name)
    stopifnot(length(cell_line_name) == 1)
    max_conc <- unique(plot_data$maxc)
    stopifnot(length(max_conc) == 1)
    
    plot_xmid <- plot_data %>% select(xmid) %>% distinct()
    plot_scal <- plot_data %>% select(scal) %>% distinct()
    plot_maxc <- plot_data %>% select(maxc) %>% distinct()
    plot_low_x <- 1 - plot_scal$scal * log10((1 - 0.001)/0.001) +  plot_xmid$xmid
    plot_low_x <- log10(getConcFromX(plot_low_x, + plot_maxc$maxc))
    plot_low_x <- min(c(all.plot_data$lx, plot_low_x))
    plot_low_x_all[ii] <- plot_low_x
    plot_high_x <- 1 - plot_scal$scal * log10(0.001/(1 - 0.001)) + plot_xmid$xmid
    plot_high_x <- log10(getConcFromX(plot_high_x, plot_maxc$maxc))
    plot_high_x <- max(c(all.plot_data$lx, plot_high_x))
    plot_high_x_all[ii] <- plot_high_x
    
    
    if (counter == 1){
      p <- ggplot(all.plot_data) + aes(x = lx, y = 1 - yhat)
      p <- p + theme_classic2() + theme(legend.background = element_rect(fill = '#FFFFFF00'))
      p <- p + annotate("rect", xmin = min(all.plot_data$lx), 
                        xmax = max(all.plot_data$lx), 
                        ymin = 0, ymax = max(c(1, 1 - all.plot_data$yhat)), alpha = 0.2)
      p <- p + ylab("Cell viability") + theme(axis.title.y = element_text(size = 12))
      p <- p + xlab(expression(log[10](Dose) ~ mu * M)) + theme(axis.title.x = element_text(size = 12))
    }
    p <- p + stat_function(data=plot_data,
                           aes(x = lx, colour = mut.status), 
                           fun = l3_model2_base10, 
                           size=0.8,
                           args = list(maxc = plot_maxc$maxc, 
                                       xmid = plot_xmid$xmid, scal = plot_scal$scal))
    
  }
  p <- p + scale_x_continuous(limits = c(min(plot_low_x_all), 
                                         max(plot_high_x_all)))
  p <- p + geom_hline(yintercept = 0.5, linetype = 3, alpha=0.5)
  return(p)
}



#### Funciton to plot nice boxplot highlighting outliers
plot.comp.boxplot.nice <- function(tissue, drug.id, gene, ref.data,
                                   highlight = NA, n.outliers = 0, 
                                   colour.by = 'tissue',
                                   add.title = TRUE) {
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
    mutate(log10.dr = log10(exp(LN_IC50))) %>%
    mutate(is.highlight=as.character(as.numeric(GDSC_cell_line %in% highlight))) %>%
    mutate(lbls=ifelse(is.highlight == '1', as.character(GDSC_cell_line), '')) %>%
    arrange(-mut.status, -log10.dr) %>%
    select(mut.status, log10.dr, is.highlight, MAX_CONC, lbls)
  
  selected.dr$lbls[(n.outliers+1):nrow(selected.dr)] <- ''
  selected.dr$is.highlight[(n.outliers+1):nrow(selected.dr)] <- '0'
  
  
  if (nrow(selected.dr) == 0){
    return(NULL)
  }
  summary.stats <- selected.dr %>%
    group_by(mut.status) %>%
    summarise(med=median(log10.dr),
              q25=quantile(log10.dr, probs = 0.25),
              q75=quantile(log10.dr, probs = 0.75))
  
  gene.label <- ifelse(grepl("mut", gene), 
                       strsplit(gene, "_mut")[[1]][1], gene)
  
  # Get max conc value (if there is a single different value, treat it as and error)
  max.conc <- log10(unique(selected.dr$MAX_CONC))
  stopifnot(length(max.conc) == 1) 
  
  target <- ifelse(unique(ref.data$PUTATIVE_TARGET[ref.data$DRUG_ID == drug.id]) != '',
                   unique(as.character(ref.data$PUTATIVE_TARGET[ref.data$DRUG_ID == drug.id])),
                   unique(as.character(ref.data$target_or_activity_of_compound[ref.data$DRUG_ID == drug.id])))
  
  if (colour.by == 'tissue') {
    g <- ggplot(selected.dr, aes(x=mut.status, y=log10.dr, 
                                 group=mut.status, colour=is.highlight)) +  
      scale_colour_manual(breaks = c('0', '1'),
                          values = c(unname(tissue.cols[tissue]), '#000000'))
  }
  
  if (colour.by == 'mut.status'){
    g <- ggplot(selected.dr, aes(x=mut.status, y=log10.dr, 
                                 group=mut.status, colour=as.factor(mut.status)))
  }
  g <- g + geom_hline(yintercept = max.conc, 
                      linetype = 2, colour = '#00000055') +
    geom_segment(data = summary.stats,
                 aes(x=mut.status-0.2, xend=mut.status+0.2,
                     y=med, yend=med),
                 color='#000000',
                 size = 1) +
    geom_segment(data = summary.stats, inherit.aes = FALSE,
                 aes(x=mut.status, xend=mut.status,
                     y=q25, yend=q75), size = 0.3,
                 color='#333333',
                 arrow = arrow(90, length = unit(0.15, "inches"), ends = 'both')) +
    geom_beeswarm(cex=4, size=2, alpha = 0.8) 
  
  if (any(!is.na(highlight))) {
    g <- g + geom_text_repel(aes(label=lbls), 
                             point.padding=0.2,
                             direction = 'y',
                             nudge_x=0.6,
                             segment.size = 0.2,
                             hjust=0,
                             min.segment.length = unit(0, 'lines'),
                             size = 3) +
      expand_limits(x=c(0, 2))
  }
  
  g <- g + theme_classic2() +
    theme(axis.line.x = element_blank(),
          legend.position = 'none') +
    scale_x_continuous(gene.label, breaks = c(0,1), labels = c('WT', 'Mut')) +
    scale_y_continuous(expression(log[10](IC[50]) ~ mu * M)) 
  if (add.title){
    g <- g + ggtitle(paste0(tissue, ' - ', unique(as.character(ref.data$DRUG_NAME[ref.data$DRUG_ID == drug.id]))),
                     paste0('[', drug.id, ']'))
  }
  return(g)
}


## Function to plot barplots with frequency of mutations
plot_mut_freq <- function(t, min_freq = 0.1) {
  tissue.cls <- all.bems.tidy %>%
    filter(Tissue==t) %>%
    select(COSMIC_ID) %>%
    distinct()
  
  msi.data <- all.gdsc %>%
    select(COSMIC_ID, Microsatellite..instability.Status..MSI.) %>%
    rename(MSI.st=Microsatellite..instability.Status..MSI.) %>%
    distinct() %>% 
    inner_join(tissue.cls, by = 'COSMIC_ID') %>%
    mutate(MSI.st=ifelse(is.na(MSI.st),
                         'Unknown', as.character(MSI.st))) %>%
    count(MSI.st) %>%
    mutate(lbls=paste0(MSI.st, '\n', n, '/', sum(n)))
  
  n.cells.tissue <- sum(msi.data$n)
  
  if (n.cells.tissue < 15){
    next
  }
  
  mut.freqs.tissue <- all.bems.tidy %>%
    filter(Tissue == t) %>%
    filter(!grepl('HypMET', alteration)) %>%
    mutate(alt.nest=map(alteration, rename.cna.mut)) %>%
    unnest() %>%
    na.omit() %>%
    filter(mut.status == 1) %>%
    select(COSMIC_ID, gene, type) %>%
    group_by(COSMIC_ID, gene) %>%
    summarise(type=paste0(type, collapse = '+')) %>% 
    ungroup() %>%
    select(-COSMIC_ID) %>%
    mutate(total=n.cells.tissue) %>%
    add_count(gene, type, name='case') %>%
    unique() %>%
    group_by(gene, type) %>%
    summarise(freq=case/total) %>%
    arrange(-freq) %>%
    ungroup()
  
  to.plot.p <- mut.freqs.tissue %>%
    group_by(gene) %>%
    mutate(total.freq=sum(freq)) %>%
    ungroup() %>%
    arrange(total.freq) %>%
    mutate(gene=fct_inorder(gene)) %>%
    filter(total.freq >= min_freq)
  
  msi.palette <- c('#6cb4b9', '#B9716C', '#777777')
  names(msi.palette) <- c('MSS/MSI-L', 'MSI-H', 'Unknown')
  pie <- ggpie(msi.data, x = 'n', 
               label = 'lbls',
               fill = 'MSI.st',
               color = NULL, lab.font = c(2, "bold", "black"),
               lab.pos = 'in', palette = msi.palette)
  pie <- ggpar(pie, legend = 'none', ggtheme=theme_void(),
               xlab = FALSE, ylab = FALSE,
               tickslab = FALSE,
               main = 'Microsatellite Status', 
               font.main = c(10, 'plain', 'black'))
  
  p <- ggplot(to.plot.p, aes(x=gene, y=freq, fill=type)) +
    geom_bar(stat = 'identity') +
    coord_flip() +
    scale_y_continuous('Mutation frequency', 
                       breaks =  seq(0, 1, by = 0.1),
                       limits = c(0, 1),
                       expand = c(0.01, 0)) +
    ggtitle(paste0(t,' (n = ', n.cells.tissue,')')) +
    theme_classic() +
    theme(axis.title.y = element_blank(),
          axis.line.x = element_blank(), 
          axis.text.y = element_text(size=10),
          panel.grid.major.x = element_line(linetype = 3, colour = '#00000068'),
          legend.position = c(0.8, 0.2)) +
    scale_fill_manual('Mutation type', 
                      values = bar.palette,
                      breaks = c('mut', 'gain', 'loss', 'mut+gain', 'mut+loss'),
                      labels = c('SNV', 'CN gain', 'CN loss', 
                                 'CN gain + SNV', 'CN loss + SNV')) +
    annotation_custom(ggplotGrob(pie), xmin  = 1, 
                      xmax=length(unique(to.plot.p$gene))/1.5,
                      ymin = 0.45, ymax = 0.95)
  return(p)
}



## Function to plot WT vs Mut boxplot with marginal density
plot_bp_and_dens <- function(tissue, drug.id, mutation){
  is.cna <- grepl('cna', mutation)
  x.title <- ifelse(is.cna, 
                    gsub('\\)', '', strsplit(mutation, '\\(')[[1]][2]),
                    gsub('_mut', '', mutation))
  x.labels <- c('WT', ifelse(is.cna,
                             paste0('CN\n', strsplit(mutation, ':')[[1]][1]),
                             'Mut'))
  max.conc <- wt.vs.mut(tissue, drug.id, mutation,
                        all.gdsc, plot.beeswarm = FALSE) %>%
    select(log10.max.conc) %>% distinct() %>% unlist() %>% unname()
  bp <- plot.comp.boxplot.nice(tissue, drug.id, mutation,
                               all.gdsc, colour.by = 'mut.status',
                               add.title = FALSE) +
    scale_color_manual(values = c('#00004488', '#FF680088')) +
    annotate('rect', xmin = -Inf, xmax = Inf, 
             ymin = -Inf, ymax = max.conc, alpha = 0.2) +
    coord_flip() +
    scale_x_reverse(x.title, 
                    breaks = c(0,1),
                    labels = x.labels)
  
  dens <- axis_canvas(bp, axis = "y",coord_flip = T) +
    geom_density(data=wt.vs.mut(tissue, drug.id, mutation,
                                all.gdsc, plot.beeswarm = FALSE),
                 aes(x=log10.dr, fill=as.factor(mut.status),
                     group=mut.status, y=..count..), 
                 alpha=0.8, adjust=2, col = NA) +
    scale_fill_manual(values = c('#00004488', '#FF680088'))
  
  bp_dens <- insert_xaxis_grob(bp, dens, position = "bottom")
  return(bp_dens)
}

