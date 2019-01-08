library(shiny)
library(shinythemes)
library(shinyWidgets)
library(data.table)
library(ggplot2)
library(ggthemes)
library(gridExtra)
library(reticulate)

# ---------------------------------------------------- #
# Configure python backend so we can use reticulate
# for loading .pkl'ed dictionaries into R as lists.
# ---------------------------------------------------- #
if (!is.null(PYTHON)) {
  use_python(PYTHON)
}

pyVers <- as.numeric(py_config()$version)
if (pyVers < 3){
  stop(paste0('Must use python3 backend (you\'re using ', pyVers, '). Use --python option.'))
}

# Load python functions with reticulate package.
source_python('reticulate_helpers.py')


# ---------------------------------------------------- #
# Assess chromosome / gene manifest from 
# supplied covcurv output directory, set global vars.
# ---------------------------------------------------- #
EXONDT <- fread(file.path(DATADIR
                          , 'gene_exon_metadata.csv'))
CHROMS <- list.files(DATADIR
                     , pattern = paste0(unique(EXONDT$chr), collapse = '|'))
readsDT <- fread(file.path(DATADIR
                           , 'read_counts.csv'))
SAMPLE_IDS <- names(readsDT)[3:ncol(readsDT)]


# ---------------------------------------------------- #
# Plotting functions
# ---------------------------------------------------- #

# find union of any overlapping exons, then find splits in union.
GetExonBreaks <- function(starts, ends){
  uVec <- numeric()
  for (i in 1:length(starts)){
    uVec <- c(uVec, starts[i]:ends[i])
  }
  
  uVec <- sort(unique(uVec))
  splits <- which(diff(uVec) > 1)
  
  if (length(splits) > 0){
    exonGrid <- c(min(uVec), uVec[splits], max(uVec))
  } else {
    exonGrid <- c(min(uVec), max(uVec))
  }
  
  return(exonGrid)
}

# format the coverage matrix from covcurve for R use.
FormatCoverageData <- function(covMat){
  # ----------------------------------------------- #
  # Cast coverage matrix into long-format data.table
  # ----------------------------------------------- #
  DT <- as.data.table(t(covMat))
  names(DT) <- SAMPLE_IDS
  
  return(DT)
}

# build a coverage curve grob.
PlotGeneCoverage <- function(covMat, geneName, chromName){
  geneLen <- ncol(covMat)
  
  # ----------------------------------------------- #
  # Build top plot
  # ----------------------------------------------- #
  covDT <- FormatCoverageData(covMat)
  
  # cast from wide to long format.
  covDT <- melt(covDT
                , id.vars = NULL
                , measure.vars = names(covDT)
                , value.name = 'coverage'
                , variable.name = 'sample')
  
  # add gene base position index (1-indexed)
  out <- covDT[, genePos := 1:nrow(.SD)
               , by=sample]
  
  pltTop <- ggplot(covDT
                   , mapping = aes(x = genePos
                                   , y = coverage
                                   , color = sample)) +
    geom_line() +
    ylab('') +
    xlab('position relative to concatentated exons') +
    theme_gdocs() +
    theme_bw() + 
    scale_color_gdocs() +
    theme(legend.position = 'bottom'
          , legend.title=element_blank()) +
    ggtitle(paste0('READS COVERAGE -- CHROM: ', chromName, '  GENE: ', geneName)) +
    xlab('')
  
  # ----------------------------------------------- #
  # Build lower plot showing exon break points
  # ----------------------------------------------- #
  geneExonDT <- EXONDT[gene == geneName]
  breaks <- GetExonBreaks(geneExonDT$start
                          , ends = geneExonDT$end)
  pltBot <- ggplot(data.frame(xmin = breaks[1]
                              , xmax = breaks[length(breaks)]
                              , ymin = 0
                              , ymax = 1)
                   , aes(xmin = xmin, xmax = xmax,ymin = ymin, ymax = ymax)) +
    geom_rect(fill = 'red')
  
  if (length(breaks) > 2){
    for (i in 2:(length(breaks) - 1)){
      pltBot <- pltBot + geom_segment(x = breaks[i]
                                      , y = 0
                                      , xend = breaks[i]
                                      , yend = 1
                                      , color = 'white')
    }
  }

  pltBot <- pltBot + xlab('') + ylab('') + 
    scale_x_continuous(breaks = c(breaks[1], breaks[length(breaks)])
                       , labels = c(breaks[1], breaks[length(breaks)])) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)
          , axis.text.y = element_blank()
          , axis.ticks.y = element_blank()
          , panel.grid = element_blank()
          , panel.border = element_blank())
    
  grob <- arrangeGrob(grobs=list(pltTop, pltBot)
                      , nrow = 2
                      , ncol = 1
                      , heights = c(3.5, 1))
  return(grob)
}


# ---------------------------------------------------- #
# Begin UI/SERVER code...
# ---------------------------------------------------- #

# Define UI for app that draws a histogram ----
ui <- fluidPage(
  tags$head(
    tags$style(
      HTML(".shiny-notification {
             position:fixed;
             width: 600px;;
             top: calc(90%);;
             left: calc(30%);;
             }"
      )
    )
    # tags$script(
    #   HTML('
    #             Shiny.addCustomMessageHandler(
    #                 type = "jsCode"
    #                 ,function(message) {
    #                 Shiny.onInputChange("confirmExportAll",eval(message.value));
    #             })
    #         ')
    # )
  ),
  
  theme = shinytheme('sandstone'),

  # App title ----
  titlePanel('covcurv'),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # render chromosome choices
      selectInput('chromosome'
                  , label = 'Select chromosome'
                  , choices = CHROMS
                  , selected = CHROMS[1]
                  , selectize = TRUE
                  , multiple = FALSE)
      
      # render *initial* gene choices based on first chromsome in gtf metadata.
      , selectInput('gene'
                  , label = 'Select gene'
                  , choices = unique(EXONDT[chr == CHROMS[1], gene])
                  , selectize = TRUE
                  , multiple = FALSE)
      , tags$br()
      , tags$br()
      , tags$br()
      , tags$br()
      , actionButton('exportAll'
                     , label = 'Export all chromosomes'
                     , style = 'color: white;
                                background-color: #F85356;')
      
      , width = 3
    ),

    # Main panel for displaying outputs ----
    mainPanel(
      
      # display coverage plot.
      # uiOutput('covPlot')
      plotOutput(outputId = 'covPlot')
      
      # fluid row to render action buttons to export plot or data.
      , fluidRow(
        column(7)
        , column(5
                 , actionButton('exportPlot'
                                , label = 'Export Plot'
                                , style = 'color: white;
                                background-color: #5098F5;')
                 , actionButton('exportData'
                                , label = 'Export Data'
                                , style = 'color: white;
                                background-color: #5098F5;')
        )
      )
    )
  )
)

# Define server logic required to draw a histogram ----
server <- function(input, output, session) {

  # update gene choices based on chromosome selection and load chromosome's coverage matrices.
  observe({
    covFile <- file.path(DATADIR
                         , input$chromosome
                         , paste0('coverage_matrices_', input$chromosome, '.pkl'))
    COV_DAT <<- load_coverage_data(covFile)
    
    updateSelectInput(session
                      , inputId = 'gene'
                      , choices = unique(EXONDT[chr == input$chromosome, gene])
    )
  })
  
  # plot gene coverage once COV_DAT is loaded.
  output$covPlot <- renderPlot({
    if (input$gene %in% names(COV_DAT)){
      grid.arrange(PlotGeneCoverage(COV_DAT[[input$gene]]
                                    , geneName = input$gene
                                    , chromName = input$chromosome))
    }
  })
  
  # export coverage plot when action button hit.
  observeEvent(input$exportPlot, {
    plt <- PlotGeneCoverage(COV_DAT[[input$gene]]
                            , geneName = input$gene
                            , chromName = input$chromosome)
    pltFile <- file.path(DATADIR
                         , input$chromosome
                         , paste0(input$gene, '_coverage.png'))
    # use relative path from DATADIR in a notification to user. 
    showNotification(paste0('Saving coverage plot to '
                            , file.path(basename(DATADIR)
                                        , input$chromosome
                                        , paste0(input$gene, '_coverage.png')))
                     , duration = 1.7)
    ggsave(pltFile
           , plot = plt
           , dpi = 300
           , width = 10
           , height = 7)
  })
  
  # export coverage matrix when action button hit.
  observeEvent(input$exportData, {
    dat <- FormatCoverageData(COV_DAT[[input$gene]])
    txtFile <- file.path(DATADIR
                         , input$chromosome
                         , paste0(input$gene, '_coverage.txt'))
    # use relative path from DATADIR in a notification to user. 
    showNotification(paste0('Saving coverage matrix to '
                            , file.path(basename(DATADIR)
                                        , input$chromosome
                                        , paste0(input$gene, '_coverage.txt')))
                     , duration = 1.7)
    write.table(dat
                , file = txtFile
                , row.names = FALSE)
  })
  
  # confirmation button when "Export all matrices" button is hit - warn user that export takes a bit.
  observeEvent(input$exportAll, {
    confirmSweetAlert(session = session
                      , inputId = 'confirmExportAll'
                      , type = 'warning'
                      , title = 'Export all coverage matrices?\n (This could take several minutes)'
                      , btn_labels = c('Cancel', 'Ok')
    )
  })
  
  # if confirmation button is OK'ed, run export.
  observeEvent(input$confirmExportAll, {
    if (input$confirmExportAll){
      withProgress(message = 'Saving all coverage matrices', style = 'notification', value = 0, max = 100, {
        incSize <- floor(100 / length(CHROMS))
        
        # iterate over chromosomes.
        for (chrom in CHROMS){
          msg <- paste0('chromosome: ', chrom)
          incProgress(incSize
                      , detail = chrom)
          
          # load pkl'ed coverage matrices into a list.
          covFile <- file.path(DATADIR
                               , chrom
                               , paste0('coverage_matrices_', chrom, '.pkl'))
          covDat <- load_coverage_data(covFile)
          
          # iterate over genes, saving one coverage matrix at a time.
          for (gene in names(covDat)){
          
            # format and save coverage matrix.
            dat <- FormatCoverageData(covDat[[gene]])
            txtFile <- file.path(DATADIR
                                 , chrom
                                 , paste0(gene, '_coverage.txt'))
            write.table(dat
                        , file = txtFile
                        , row.names = FALSE)
          }
        }
      })
      
      # Tell user the export was successful.
      showNotification(paste0('Success! All coverage matrices saved to directory ', basename(DATADIR))
                       , duration = 4)
    }
  })

}


# instantiate application.
shinyApp(ui = ui, server = server)



