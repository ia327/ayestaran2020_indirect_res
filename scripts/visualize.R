library(shiny)
library(beeswarm)
library(scales)
suppressMessages(library(tidyverse))
suppressMessages(library(plotly))

## FIXME: Add result of ANOVA, with reactive boxplot comparison
## Add reactive boxplot comparison to ANOVA filtered volcano
## Idea, just one volcano plot, with option to toggle filtered vs non filtered?

DATADIR <- '../data/combined_for_pipeline/'
RESDIR <- '../results/'

load(paste0(DATADIR, 'gdsc_dr_data.RData'))
load(paste0(DATADIR, 'ctrp_dr_data.RData'))
load(paste0(DATADIR, 'utils.RData'))
load(paste0(DATADIR, 'bems_tidy.RData'))
load(paste0(DATADIR, 'wes_annotation.RData'))
load(paste0(DATADIR, 'cn_annotation.RData'))
load(paste0(DATADIR, 'cancer_gene_lists.RData'))


load(paste0(RESDIR, 'ctrp/sign_anova_out.RData'))
ctrp.res.list <- list(sign.anova=selected.sign.anovas.df,
                      IC50s=selected.all.IC50s.df)
load(paste0(RESDIR, 'ctrp/putative_res_markers.RData'))
ctrp.assoc.filt <- associations %>% arrange(outlier.q.value, n.outliers)
ctrp.assoc.filt$clean.name <- unname(apply(ctrp.assoc.filt, 1,
                                           function(yy) {
                                             drug.name <- all.ctrp %>%
                                               filter(DRUG_ID == yy['drug']) %>%
                                               dplyr::select(DRUG_NAME) %>%
                                               unique() %>% unlist() %>% unname()
                                             paste0(yy['tissue'], ' - ',
                                                    drug.name, ' - ', yy['alteration'],
                                                    ' : ', yy['n.outliers'], ' out')
                                           }))
ctrp.assoc.filt$clean.name <- paste0(1:nrow(ctrp.assoc.filt), '. ',
                                     ctrp.assoc.filt$clean.name)
load(paste0(RESDIR, 'ctrp/fdr.tree.value.RData'))
fdr.tree.ctrp <- fdr.tree

load(paste0(RESDIR, 'gdsc/sign_anova_out.RData'))
gdsc.res.list <- list(sign.anova=selected.sign.anovas.df,
                      IC50s=selected.all.IC50s.df)
load(paste0(RESDIR, 'gdsc/putative_res_markers.RData'))
gdsc.assoc.filt <- associations %>% arrange(outlier.q.value, n.outliers)
gdsc.assoc.filt$clean.name <- unname(apply(gdsc.assoc.filt, 1,
                                           function(yy) {
                                             drug.name <- all.gdsc %>%
                                               filter(DRUG_ID == yy['drug']) %>%
                                               dplyr::select(DRUG_NAME) %>%
                                               unique() %>% unlist() %>% unname()
                                             paste0(yy['tissue'], ' - ',
                                                    drug.name, ' - ', yy['alteration'],
                                                    ' : ', yy['n.outliers'], ' out')
                                           }))
gdsc.assoc.filt$clean.name <- paste0(1:nrow(gdsc.assoc.filt), '. ',
                                     gdsc.assoc.filt$clean.name)
load(paste0(RESDIR, 'gdsc/fdr.tree.value.RData'))
fdr.tree.gdsc <- fdr.tree



plot.comp.boxplot <- function(tissue, drug.id, gene, ref.data, ref.bems, outliers = 0){
  # Get BEM
  simple.bem <- ref.bems %>%
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

  gene.label <- ifelse(grepl("mut", gene), strsplit(gene,
                                                    "_mut")[[1]][1], gene)

  selected.dr <- beeswarm(log10.dr ~ mut.status, data=selected.dr,
                          do.plot = FALSE, spacing=1) %>%
    mutate(x.orig=as.numeric(x.orig)) %>%
    left_join(selected.dr, by = c('x.orig'= 'mut.status', 'y.orig' = 'log10.dr')) %>%
    mutate(x.orig=x.orig+1)

  target <- ifelse(unique(ref.data$PUTATIVE_TARGET[ref.data$DRUG_ID == drug.id]) != '',
                   unique(as.character(ref.data$PUTATIVE_TARGET[ref.data$DRUG_ID == drug.id])),
                   unique(as.character(ref.data$target_or_activity_of_compound[ref.data$DRUG_ID == drug.id])))
  drug.label.title <- paste0(tissue, " - ",
                             unique(as.character(ref.data$DRUG_NAME[ref.data$DRUG_ID == drug.id]))
                             , " (", target,")")

  min.y.axis <- min(c(selected.dr$y.orig, max.conc)) - 0.5
  max.y.axis <- max(c(selected.dr$y.orig, max.conc)) + 0.5

  g <- plot_ly(data = selected.dr, type = 'box') %>%
    add_boxplot(y=~y.orig, x=~x.orig,
                marker = list(color = unname(tissue.cols[tissue])),
                line = list(color = tissue.cols[tissue]),
                hoverinfo = 'none') %>%
    layout(title = drug.label.title,
           margin = list(
             l = 50,
             r = 50,
             b = 100,
             t = 100,
             pad = 4
           ),
           showlegend = FALSE,
           xaxis = list(
             title = paste0(gene.label, ' mutation'),
             ticktext = c('Absent', 'Present'),
             tickvals = 1:2),
           yaxis = list(title = 'log10(IC50)',
                        range = c(min.y.axis, max.y.axis)),
           shapes=list(type='line',
                       x0= 0.5, x1= 2.5,
                       y0=max.conc,
                       y1=max.conc,
                       line=list(dash='dot', width=1))
    ) %>%
    add_markers(data = selected.dr, x=~x, y=~y,
                hoverinfo = 'text',
                marker = list(color = tissue.cols[tissue]),
                text = ~GDSC_cell_line)
  if (outliers > 0){
    outlier.dr <- selected.dr %>% filter(x.orig == 2) %>% arrange(y.orig)
    outlier.dr <- outlier.dr[(nrow(outlier.dr) - (outliers-1)):nrow(outlier.dr),]

    g <- g %>%
      add_markers(data = outlier.dr, x=~x, y=~y,
                  marker = list(color = 'black',
                                size = 8),
                  hoverinfo = 'none')
  }
  g$elementId <- NULL
  return(g)
}


# User interface ----
ui <- fluidPage(
  titlePanel("Drug response comparisons"),
  fluidRow(wellPanel(selectInput("dataset",
                                 label = "Which dataset should be considered?",
                                 choices = c('GDSC',
                                             'CTRP'),
                                 selected = "GDSC"))),
  tabsetPanel(
    tabPanel('All associations',
             column(4,
                    uiOutput("tissueSelection"),
                    uiOutput("geneSelection"),
                    uiOutput("drugSelection")
             ),
             column(5,
                    plotlyOutput("plot.b", height = '480px'))
    ),
    tabPanel('Putative markers',
             column(4,
                    uiOutput("putativeSelection")
             ),
             column(4,
                    plotlyOutput("plot.putative", height = '480px')),
             column(4,
                    uiOutput('putative.info'))
    ),
    tabPanel('ANOVA filtered volcano plots',
             plotlyOutput('plot.volcano')
    ),
    tabPanel('Putative volcano',
             column(6,
                    plotlyOutput('putative.volcano')),
             column(6,
                    plotlyOutput('putative.boxplot.reactive'),
                    htmlOutput('marker.description'))
    )
  )
)

# Server logic ----
server <- function(input, output) {
  ref.data <- reactive({
    switch(input$dataset,
           GDSC = all.gdsc,
           CTRP = all.ctrp)
  })
  
  ref.bem <- reactive({
    switch(input$dataset,
           GDSC = all.bems.tidy.gdsc,
           CTRP = all.bems.tidy.ctrp)
  })

  putative.assoc.data <- reactive({
    switch(input$dataset,
           GDSC = gdsc.assoc.filt,
           CTRP = ctrp.assoc.filt)
  })

  anova.sign.data <- reactive({
    switch(input$dataset,
           GDSC = gdsc.res.list$sign.anova,
           CTRP = ctrp.res.list$sign.anova)
  })
  
  fdr.tree <- reactive({
    switch(input$dataset,
           GDSC = fdr.tree.gdsc,
           CTRP = fdr.tree.ctrp)
  })

  output$tissueSelection <- renderUI({
    dat <- ref.data()
    selectInput('tissueSelection', 'Choose tissue',
                choices = all.tissues,
                selectize = TRUE)
  })

  output$geneSelection <- renderUI({
    simple.bem <- ref.bem() %>%
      filter(Tissue == input$tissueSelection)
    selectInput('geneSelection', 'Choose gene',
                choices = simple.bem$alteration,
                selectize = TRUE)
  })

  output$drugSelection <- renderUI({
    dat <- ref.data()
    drug.ids <- dat$DRUG_ID
    names(drug.ids) <- paste0(dat$DRUG_NAME,' (' , dat$DRUG_ID, ')')
    drug.ids <- drug.ids[!duplicated(drug.ids)]
    paste0(dat$DRUG_NAME,' (' , dat$DRUG_ID, ')')
    selectInput('drugSelection', 'Choose drug',
                choices = drug.ids,
                selectize = TRUE)
  })

  output$putativeSelection <- renderUI({
    assoc.dat <- putative.assoc.data()


    selectInput('putativeSelection', 'Choose association',
                choices = assoc.dat$clean.name,
                selectize = TRUE)
  })


  output$plot.b <- renderPlotly({
    dat <- ref.data()
    drug.id <- input$drugSelection

    if (length(drug.id) == 0){
      plotly_empty()
    } else {
      all.plots <- vector(mode='list', length = length(drug.id))
      for (i in 1:length(all.plots)){
        all.plots[[i]] <- plot.comp.boxplot(tissue = input$tissueSelection,
                                            drug.id = drug.id[i],
                                            gene = input$geneSelection,
                                            ref.bems = ref.bem(),
                                            ref.data = ref.data())
      }
      subplot(all.plots, shareX = FALSE, shareY = TRUE)
    }
  })

  output$plot.putative <- renderPlotly({
    assoc.dat <- putative.assoc.data()
    assoc.id <- input$putativeSelection
    chosen.assoc <- assoc.dat %>% filter(clean.name == assoc.id)

    if (length(assoc.id) == 0){
      plotly_empty()
    } else {
      plot.comp.boxplot(tissue = chosen.assoc$tissue,
                        drug.id = chosen.assoc$drug,
                        gene = chosen.assoc$alteration,
                        ref.data = ref.data(),
                        ref.bems = ref.bem(),
                        outliers = chosen.assoc$n.outliers)
    }
  })

  output$putative.info <- renderTable({
    assoc.dat <- putative.assoc.data()
    assoc.id <- input$putativeSelection
    chosen.assoc <- assoc.dat %>% filter(clean.name == assoc.id) %>%
      dplyr::select(c(outlier.p.value,
                      outlier.q.value, putative.markers)) %>%
      mutate(outlier.p.value=ifelse(outlier.p.value < 0.01,
                                    scientific(outlier.p.value),
                                    as.character(round(outlier.p.value, 2))),
             outlier.q.value=ifelse(outlier.q.value < 0.01,
                                    scientific(outlier.q.value),
                                    as.character(round(outlier.q.value, 2))))
    chosen.assoc
  })


  output$plot.volcano <- renderPlotly({
    dat <- ref.data()
    anova.out <- anova.sign.data()
    anova.out <- left_join(anova.out,
                           unique(dplyr::select(dat, DRUG_ID, DRUG_NAME)),
                           by = c('drug' = 'DRUG_ID')) %>%
      mutate(drug=as.numeric(drug),
             tissue=as.factor(tissue))

      plot_ly(data = anova.out,
              x=~effect,
              y=~corrected.p.value...q.value.,
              color=~tissue,
              colors = tissue.cols[levels(anova.out$tissue)],
              size=~n.samples,
              sizes = c(20, 50)) %>%
      add_trace(mode = 'markers',
                hoverinfo = 'text',
                opacity = 0.8,
                text = ~paste(DRUG_NAME, '-',
                              alteration, '\n',
                              tissue, '\n',
                              'FDR = ', scientific(corrected.p.value...q.value.),
                              'Effect = ', round(effect, 3))) %>%
      layout(yaxis=list(autorange = 'reversed',
                        type = 'log',
                        title = 'FDR',
                        exponentformat='E'))
  })

  output$putative.volcano <- renderPlotly({
    dat <- ref.data()
    putative.out <- putative.assoc.data()
    fdr.thres <- 0.15
    putative.out <- putative.out %>%
      mutate(tissue=ifelse(outlier.q.value > fdr.thres,
                           'Not significant', tissue),
             tissue=as.factor(tissue)) %>%
      left_join(unique(dplyr::select(dat, DRUG_ID, DRUG_NAME)),
                by = c('drug' = 'DRUG_ID')) %>%
      rownames_to_column(var = 'key.id')

    plot_ly(data = putative.out,
            x=~effect.size,
            y=~outlier.q.value,
            color=~tissue,
            colors = unname(tissue.cols[levels(putative.out$tissue)]),
            size=~n.outliers,
            sizes = c(20, 50),
            key = ~key.id,
            source = 'source',
            legendgroup = putative.out$tissue) %>%
      add_trace(mode = 'markers',
                hoverinfo = 'text',
                opacity = 0.9,
                text = ~paste0(DRUG_NAME,' (' , drug, ')\n',
                              alteration, '\n',
                              tissue, '\n',
                              'FDR = ', scientific(outlier.q.value),
                              '\nN. outliers = ', n.outliers)) %>%
      add_lines(hoverinfo = 'none',
                split=~association,
                line=list(width = 1,
                          alpha = 0.3,
                          dash = 'dot'),
                showlegend = FALSE) %>%
      layout(yaxis=list(autorange = 'reversed',
                        type = 'log',
                        title = 'Outlier P value',
                        exponentformat='E'),
             xaxis=list(title = 'Decrease in SD'))
  })

  output$putative.boxplot.reactive <- renderPlotly({
    assoc.dat <- putative.assoc.data()
    eventdata <- event_data('plotly_click', source = 'source')
    validate(need(!is.null(eventdata),
                  "Click on any point to show the corresponding boxplot"))
    assoc.id <-  as.numeric(eventdata$key)[1]
    chosen.assoc <- assoc.dat[assoc.id,]

    if (length(assoc.id) == 0){
      plotly_empty()
    } else {
      plot.comp.boxplot(tissue = chosen.assoc$tissue,
                        drug.id = chosen.assoc$drug,
                        gene = chosen.assoc$alteration,
                        ref.data = ref.data(),
                        ref.bems = ref.bem(),
                        outliers = chosen.assoc$n.outliers)
    }
  })

  output$marker.description <- renderTable({
    assoc.dat <- putative.assoc.data()
    eventdata <- event_data('plotly_click', source = 'source')
    validate(need(!is.null(eventdata),
                  "Click on any point to show the corresponding putative resistance markers"))
    assoc.id <-  as.numeric(eventdata$key)[1]
    chosen.assoc <- assoc.dat[assoc.id,]
    chosen.assoc[c('out.cell.lines', 'putative.markers')]

  })


}
# Run app ----
shinyApp(ui, server)
