#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#
suppressPackageStartupMessages({
  library(shiny)
  library(shinyjs)
  library(shinydashboard)
  library(shinycssloaders)
  library(DESeq2)
  library(pasilla)
  library(vsn)
  library(tidyverse)
  library(shinyvalidate)
  library(ggplot2)
  library(tibble)
})
library(conflicted)
conflict_prefer("select", "dplyr")
conflicts_prefer(shinydashboard::box)


random_icon <- sample(c("canadian-maple-leaf", "dragon", "user", "cog", 
                        "dice-d20", "dumpster-fire", "pastafarianism"), 1)
skin_color <- sample(c("red", "yellow", "green", "blue", "purple", "black"), 1)

normalization_methods <- list( # no blank needed
  "Total Sum Scaling" = "tss",
  "Common Sum Scaling" = "css",
  "Center Log Ratio" = "clr",
  "DESeq2 VST" = "vst",
  "DESeq2 Regularized Log Transform" = "rlog"
)

# Define UI for application that draws a histogram
ui <- dashboardPage(
  skin = skin_color,
  
  dashboardHeader(title = "Basic Template", 
                  titleWidth = "calc(100% - 44px)" # puts sidebar toggle on right
  ),
  
  dashboardSidebar(
    # https://fontawesome.com/icons?d=gallery&m=free
    sidebarMenu(
      id = "tabs",
      menuItem("Demo Tab", tabName = "dataset_info_tab", icon = icon(random_icon))
    ),
    tags$a(href = "https://debruine.github.io/shinyintro/", 
           "ShinyIntro book", style="padding: 1em;")
  ),
  
  dashboardBody(
    shinyjs::useShinyjs(),
    tags$head(
      # links to files in www/
      tags$link(rel = "stylesheet", type = "text/css", href = "basic_template.css"), 
      tags$link(rel = "stylesheet", type = "text/css", href = "custom.css"), 
      tags$script(src = "custom.js")
    ),
    tabItems(
      tabItem(
        tabName = "dataset_info_tab",
        
        # Row 1: info -> controls -> info
        fluidRow(
          infoBox("Overview of the dataset",
                  textOutput("sample_info"),
                  width = 4, color = "red", icon = shiny::icon("leaf")),
          box(title = "Filter by minimum total count per gene",
              numericInput("n_min",
                           label = "Enter the minimum number of raw counts for filtering out uninformative genes",
                           value = 10, min = 1, max = 10000, step = 1),
              actionButton("update_n_min", "Update n_min!"),
              width = 4),
          infoBox(textOutput("n_min"),
                  textOutput("filtered_info"),
                  width = 4, color = "red", icon = shiny::icon("leaf"))
        ),
        
        # Row 2: two plots side-by-side and MA plot
        fluidRow(
          column(6,
            box(title = "Mean-SD Plot",
              withSpinner(plotOutput("filt_ms_plot", height = 400)), width =12)
          ),
          column(6,
            box(title = "MA Plots",
              withSpinner(plotOutput("ma_plot", height = 400)), width = 12)
          )
        ),

        # Row 4: full-width controls
        fluidRow(
          box(title = "Choose normalization method to compare",
              selectInput("normalization",
                          label   = "Normalization Method",
                          choices = normalization_methods,
                          selected = NULL),
              width = 12)
        )
      )
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    # Load pasilla data, extract raw counts and column names (samples)
    fn <- system.file("extdata", "pasilla_gene_counts.tsv",
                      package = "pasilla", mustWork = TRUE)
    counts <- read.csv(fn, sep = "\t", row.names = "gene_id")
    sample_columns <- colnames(counts)
    rs <- rowSums(counts) # cache for speed
    
    # Output basic info
    output$sample_info <- renderText(
      paste0("There are ", ncol(counts), 
             " samples and ",nrow(counts), " total genes")
    )
    
    # Filter out genes with low count.
    # First validate the n_min input.
    iv_n_min <- shinyvalidate::InputValidator$new()
    iv_n_min$add_rule("n_min", shinyvalidate::sv_integer())
    iv_n_min$add_rule("n_min", shinyvalidate::sv_between(1, 1000))
    iv_n_min$enable()
    
    # Define the reactive values
    v <- reactiveValues(
      n_min = 10L,
      title = "",
      keep = rep(TRUE, nrow(counts)),
      counts_filt = NULL
      # raw_df = NULL,
      # filt_df = NULL,
    )
    
    # init_nmin <- 10L
    # init_compute(init_nmin)
    
    # Commit on button click
    observeEvent(input$update_n_min, {
      if (!iv_n_min$is_valid()) return(NULL)
      v$title <- "Filtering and building plots..." 
      nmin <- as.integer(input$n_min)
      v$n_min <- nmin
      # Recompute derived data
      init_compute(nmin)
    }, ignoreInit = TRUE)
    
    # Number of genes after filtering output
    output$filtered_info <- renderText(v$title)
    
    # Mean-SD Plot
    output$filt_ms_plot <- renderPlot({
      req(v$counts_filt)
      meanSdPlot(as.matrix(counts_filt), ranks = FALSE)
    })
    
    # Helper for MA Plot, prepares the data.
    ma_df <- reactive({
      req(v$counts_filt)  # depends on committed filter
      cm <- v$counts_filt
      # Build MA for each column and bind
      dplyr::bind_rows(lapply(colnames(cm), MA_comparison, data_matrix = cm))
    })
    
    # MA Plots        
    output$ma_plot <- renderPlot({
      df <- ma_df()
      MA_Plot(df)
    })
    
    # Functions:
    # One-shot initializer (optional): build initial plots at startup
    init_compute <- function(n_min_val) {
      keep <- rs >= n_min_val
      raw_mat <- as.matrix(counts)
      v$counts_filt <- counts[keep,] # Keep the info of these genes
      filt_mat <- as.matrix(counts_filt)
      v$keep <- keep
      v$title <- sprintf("After filtering (n_min = %d): %d genes remain",
                         n_min_val, nrow(filt_mat))
      # v$raw_df <- tibble::tibble(
      #   mean = rowMeans(raw_mat),
      #   sd = matrixStats::rowSds(raw_mat),
      #   kept = keep
      # )
      # v$filt_df <- tibble::tibble(
      #   mean = rowMeans(filt_mat),
      #   sd = matrixStats::rowSds(filt_mat)
      # )
    }
    
}

# Function for MA-plots, which is inspired by the function in edgeR
MA_comparison <- function(ref_sample, data_matrix,
                          pseudocount=1, data_is_log=FALSE){
  ##############################################################################
  # ref_sample: the column we are creating the plot for
  # data_matrix: the complete dataset.
  # pseudocount: since we are transforming counts into log space, it is useful
  # to add a pseudocount to avoid NAs
  # data_is_log: a simple flag to indicate whether the counts in
  #‘data_matrix‘ are already in log space or not.
  ##############################################################################
  if(data_is_log){
    ref_data <- data_matrix %>% select(eval(ref_sample))
    rest_data <- data_matrix %>% select(-eval(ref_sample)) %>% rowMeans()
  }else if(data_is_log == FALSE){
    d_pseudo <- data_matrix + 1 # We add a pseudocount to the 0 values
    ref_data <- d_pseudo %>% select(eval(ref_sample)) %>% log2()
    # We collect the data for our column of interest and transform to log2
    rest_data <- d_pseudo %>% select(-eval(ref_sample)) %>% rowMeans() %>% log2()
    # We collect the average in all other samples and conver to log2
  }
  fc <- ref_data[,ref_sample] - rest_data
  # compute the fold change show on the y-axis
  out.df <- data.frame("avg.expression"=rest_data,
                       # the average gene expression show on x-axis
                       "logfc" = fc, # the fold change show on y-axis
                       "sample" = ref_sample)
  return(out.df)
}

MA_Plot <- function(df){
  ggplot(df, aes(x=avg.expression, y=logfc)) +
    geom_point(color=alpha("grey", 0.7)) +
    geom_hline(yintercept = 0, color="red") +
    geom_smooth(aes(x=avg.expression,y=logfc), method = "loess",
                color="blue", linewidth=0.5) +
    facet_wrap(~sample) +
    labs(x="Average Expression (on log scale)",
         y="log Fold Change",
         title="MA-plot of raw counts")
}

# Run the application 
shinyApp(ui = ui, server = server)

# fluidPage(
#   
#   # Application title
#   titlePanel("Old Faithful Geyser Data"),
#   
#   # Sidebar with a slider input for number of bins 
#   sidebarLayout(
#     sidebarPanel(
#       sliderInput("bins",
#                   "Number of bins:",
#                   min = 1,
#                   max = 50,
#                   value = 30)
#     ),
#     
#     # Show a plot of the generated distribution
#     mainPanel(
#       plotOutput("distPlot")
#       
#     )
#   )
# )
# 
# n_min <- eventReactive(input$update_n_min, {
#   req(iv_n_min$is_valid()) # guard the update by the validator
#   as.integer(input$n_min)
# }, ignoreInit = TRUE
# )
# counts_filt <- reactive( { 
#   req(!is.null(n_min()))
#   keep <- rowSums(counts) >= n_min() 
#   counts[keep, ]
# }
# )
# # Print the info for filtered out data with renderText
# output$filtered_info <- renderText({
#   req(!is.null(counts_filt()))
#   paste0("After filtering with respect to total counts (n_min=", 
#          n_min(), ") ",
#          nrow(counts_filt()), 
#          " genes remain for analysis")
# }
# )
# 
# output$n_min <- renderText(
#   paste0("")
# )



# # raw mean–SD with heatmap + filtered overlays
# output$raw_ms_plot <- renderPlot({
#   req(v$raw_df, v$keep)
#   ggplot(v$raw_df, ggplot2::aes(mean, sd)) +
#     geom_hex(bins = 50) +
#     scale_fill_viridis_c(name = "Count") +
#     geom_point(
#       data = subset(v$raw_df, !kept),
#       ggplot2::aes(mean, sd),
#       color = "#E4572E", size = 1.5, alpha = 0.9
#     ) +
#     labs(x = "Average Expression (all)", y = "Standard Deviation (all)", title = "Mean-SD Plot of Raw Counts") +
#     theme_minimal(base_size = 12) +
#     theme(legend.position = "top")
# })


# # filtered-only heatmap (same style)
# output$filt_ms_plot <- renderPlot({
#   req(v$filt_df)
#   ggplot(v$filt_df, ggplot2::aes(mean, sd)) +
#     geom_hex(bins = 50) +
#     scale_fill_viridis_c(name = "Count") +
#     theme_minimal(base_size = 12) +
#     theme(legend.position = "top") +
#     labs(x = "Average Expression (filtered)", y = "Standard Deviation (filtered)", title = "Mean-SD Plot of Raw Counts")
# })
