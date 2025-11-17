#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

# nolint start
suppressPackageStartupMessages({
  library(shiny)
  library(shinyjs)
  library(shinydashboard)
  library(shinycssloaders)
  library(DESeq2)
  library(pasilla)
  library(vsn)
  library(compositions)
  library(tidyverse)
  library(shinyvalidate)
  library(ggplot2)
  library(tibble)
  library(future)
  library(promises)
  library(bslib)
  library(plotly)  # ← Add for interactive PCA plots
  library(ggrepel) # ← Add for better label placement
})
library(conflicted)
conflict_prefer("select", "dplyr")
conflicts_prefer(shinydashboard::box)

future::plan(multisession)

random_icon <- sample(c("canadian-maple-leaf"))
                        # "dragon", "user", "cog", 
                        # "dice-d20", "dumpster-fire", "pastafarianism"), 1)
skin_color <- sample(c("green")) #, "yellow", "green", "blue", "purple", "black"), 1)

normalization_methods <- list( # no blank needed
  "Total Sum Scaling" = "tss",
  "Common Sum Scaling" = "css",
  "Center Log Ratio" = "clr",
  "DESeq2 Median of Ratios" = "mor",
  "DESeq2 VST" = "vst",
  "DESeq2 Regularized Log Transform" = "rlog"
)

# Define UI for application that draws a histogram
ui <- dashboardPage(
  skin = skin_color,

  # INITIALIZING DASHBOARD STRUCTURE
  dashboardHeader(
    title = "Basic Template", 
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
            numericInput("n_min_input",
              label = "Enter the minimum number of raw counts for filtering out uninformative genes",
              value = 10, min = 1, max = 10000, step = 1),
            input_task_button("update_n_min", "Update n_min!"),
            width = 4),
          infoBox(textOutput("n_min"), 
            textOutput("filtered_info"),
            width = 4, color = "red", icon = shiny::icon("leaf"))
        ),
        
        # Row 2: two plots side-by-side and MA plot
        fluidRow(
          column(6,
            box(title = "Mean-SD Plot",
              withSpinner(plotOutput("filt_ms_plot", height = 400)), width = 12)
          ),
          column(6,
            box(title = "MA Plots",
              withSpinner(plotOutput("ma_plot", height = 400)), width = 12)
          )
        ),

        # Row 4: Normalization method selection
        fluidRow(
          box(title = "Choose normalization method to compare",
              checkboxGroupInput("normalization",
                          label   = "Normalization Method",
                          choices = normalization_methods,
                          selected = NULL),
              input_task_button("normalize_button", "Normalize and plot!"),  # Initially disabled
              input_task_button("plot_button", "Generate Plots"),  # ← NEW BUTTON
              width = 12)
        ),

        # In UI, add an output area (after line 126):
        fluidRow(
          box(
            title = "Normalization Status",
            uiOutput("norm_status"),
            width = 12
          )
        ),

        # Row 5: Mean-SD and MA plots after normalization
        uiOutput("norm_plots_ui")
      )
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  # Load pasilla data, extract raw counts and column names (samples)
  fn <- system.file("extdata", "pasilla_gene_counts.tsv",
                    package = "pasilla", mustWork = TRUE)
  counts <- read.csv(fn, sep = "\t", row.names = "gene_id")
  sample_columns <- colnames(counts)
  rs <- rowSums(counts) # cache for speed
  
  # Annotation data (not used in this simple app)
  # Load Pasilla sample annotation data
  annotationFile = system.file("extdata",
    "pasilla_sample_annotation.csv",
    package = "pasilla", mustWork = TRUE)

  # Create sample annotation dataframe
  pasillaSampleAnno.df = read_csv(annotationFile) %>%
    mutate(sample = sub("fb", "", file),
    type = sub("-.*", "", type)) %>%
    column_to_rownames("sample")

  # Output basic info
  output$sample_info <- renderText(
    paste0("There are ", ncol(counts), 
          " samples and ",nrow(counts), " total genes")
  )
  
  # Filter out genes with low count.
  # First validate the n_min input.
  iv_n_min <- shinyvalidate::InputValidator$new()
  iv_n_min$add_rule("n_min_input", shinyvalidate::sv_integer())
  iv_n_min$add_rule("n_min_input", shinyvalidate::sv_between(1, 1000))
  iv_n_min$enable()
  
  # Define the reactive values
  v <- reactiveValues(
    n_min = 10L,
    title = "",
    keep = rep(TRUE, nrow(counts)),
    counts_filt = NULL
  )
  
  # Reactive event: Update n_min -> Filter wrt n_min and display message.
  observeEvent(input$update_n_min, {
    if (!iv_n_min$is_valid()) return(NULL)
    v$title <- "Filtering and building plots..." 
    v$n_min <- as.integer(input$n_min_input)
    filter_func(counts, v) # Perform filtering on the raw count data.
    req(v$counts_filt)
    counts_filt_snapshot <- v$counts_filt # snapshot for async use
    session$onFlushed(function() {
      meanSdPlotTask$invoke(counts_filt_snapshot) # Build Mean-SD plot
      MAPlotTask$invoke(counts_filt_snapshot)     # Build MA plots
    }, once = TRUE)
  }, ignoreInit = TRUE)
      
  # Number of genes after filtering output
  output$filtered_info <- renderText({
    v$title
  })
  
  # Mean-SD Plot
  output$filt_ms_plot <- renderPlot({
    meanSdPlotTask$result()
  })
  
  # MA Plots        
  output$ma_plot <- renderPlot({
    MAPlotTask$result()
  })

  # ========== NORMALIZATION SECTION ==========
  # Storage for normalization results
  v$norm_results <- reactiveValues()
  v$norm_methods_committed <- NULL

  # Create ExtendedTasks for all normalization plot methods (at startup)
  norm_plot_tasks <- list()
  for (method in c("tss", "css", "clr", "mor", "vst", "rlog")) {
    local({
      m <- method # Set local copy for closures

      # Mean-SD plot task
      norm_plot_tasks[[paste0("ms_", m)]] <<- ExtendedTask$new(
        function(norm_counts) {
          future_promise({
            meanSdPlot(as.matrix(norm_counts), ranks = FALSE)
          })
        }
      ) |> bind_task_button("plot_button")
      
      # MA plot task
      norm_plot_tasks[[paste0("ma_", m)]] <<- ExtendedTask$new(
        function(norm_counts) {
          future_promise({
            data_is_log <- m %in% c("clr", "vst", "rlog")
            ma_df <- dplyr::bind_rows(lapply(
              colnames(norm_counts),
              MA_comparison,
              data_matrix = norm_counts,
              data_is_log = data_is_log
            ))
            MA_Plot(ma_df)
          })
        }
      ) |> bind_task_button("plot_button")
      
      # Render Mean-SD plot from task result
      output[[paste0("ms_plot_", m)]] <- renderPlot({
        norm_plot_tasks[[paste0("ms_", m)]]$result()
      })
      
      # Render MA plot from task result
      output[[paste0("ma_plot_", m)]] <- renderPlot({
        norm_plot_tasks[[paste0("ma_", m)]]$result()
      })
    })
  }

  # ========== NORMALIZE DATA BUTTON ==========
  observeEvent(input$normalize_button, {
    # Get selected methods
    methods <- input$normalization
    req(length(methods) > 0, v$counts_filt) # Require at least one method and filtered counts

    # Disable button during processing
    shinyjs::disable("normalize_button")
    on.exit(shinyjs::enable("normalize_button"), add = TRUE)

    # Store committed methods in reactiveValues
    v$norm_methods_committed <- methods

    # Snapshot of filtered counts
    counts_snapshot <- v$counts_filt

    # Apply normalization for each method
    for (method in methods) {
      result <- tryCatch({
          norm_counts <- normalization_func(counts_snapshot, method, metadata.df=pasillaSampleAnno.df)
          list(norm_counts = norm_counts, error = NULL) # No error
        }, error = function(e) {
          list(norm_counts = NULL, error = conditionMessage(e)) # Capture error message
        })
      
      # Store results in reactiveValues
      v$norm_results[[method]] <- result
    }

    # Enable plot button after normalization completes
    shinyjs::enable("plot_button")
  })

  # ========== INITIALIZE PLOT BUTTON STATE ==========
  observe({
    # Disable plot button initially
    shinyjs::disable("plot_button")
  })

  # ========== GENERATE PLOT BUTTON ==========
  observeEvent(input$plot_button, {
    methods <- v$norm_methods_committed
    req(length(methods) > 0)
    
    # # Disable button during plotting
    # shinyjs::disable("plot_button")
    # on.exit(shinyjs::enable("plot_button"), add = TRUE)
    
    # Invoke plotting tasks for each normalized method (ASYNCHRONOUS)
    for (method in methods) {
      result <- v$norm_results[[method]]
      req(!is.null(result))
      
      # If successful, invoke plotting tasks
      if (is.null(result$error)) {
        local({
          m <- method
          nc <- result$norm_counts
          # Invoke tasks after reactive cycle completes
          session$onFlushed(function() {
            norm_plot_tasks[[paste0("ms_", m)]]$invoke(nc)
            norm_plot_tasks[[paste0("ma_", m)]]$invoke(nc)
          }, once = TRUE)
        })
      }
    }
  })

  # ========== DYNAMIC UI SUCCESS MESSAGE ==========
  output$norm_status <- renderUI({
    methods <- v$norm_methods_committed
    req(length(methods) > 0)
    
    # Get success/failure status
    success_index <- sapply(methods, function(m) {
      res <- v$norm_results[[m]]
      !is.null(res) && is.null(res$error)
    })
    
    # Extract method codes
    success_methods <- methods[success_index]
    failed_methods <- methods[!success_index]
    
    # Convert to human-readable labels
    get_labels <- function(method_codes) {
      sapply(method_codes, function(code) {
        names(normalization_methods)[normalization_methods == code]
      })
    }
    
    success_labels <- get_labels(success_methods)
    failed_labels <- get_labels(failed_methods)
    
    div(
      if (length(success_methods) > 0) {
        tags$p(style = "color: green; font-weight: bold;", 
              icon("check-circle"), 
              paste("Successfully normalized:", paste(success_labels, collapse = ", ")))
      },
      if (length(failed_methods) > 0) {
        tags$p(style = "color: red; font-weight: bold;",
              icon("exclamation-circle"),
              paste("Failed:", paste(failed_labels, collapse = ", ")))
      }
    )
  })

  # ========== DYNAMIC UI FOR PLOTS ==========
  output$norm_plots_ui <- renderUI({
    methods <- v$norm_methods_committed
    req(length(methods) > 0)

    do.call(tagList, lapply(methods, function(method) {
      label <- names(normalization_methods)[normalization_methods == method]
      res <- v$norm_results[[method]]
      
      # Show error message if normalization failed
      if (!is.null(res) && !is.null(res$error)) {
        fluidRow(
          box(
            title = paste("Error -", label),
            div(style = "color: red; padding: 20px;",
                icon("exclamation-triangle"),
                paste("Error:", res$error)),
            width = 12
          )
        )
      } else {
        # Show plots (with spinners while tasks are running)
        fluidRow(
          box(
            title = paste("Mean-SD Plot -", label),
            withSpinner(plotOutput(paste0("ms_plot_", method), height = 400)),
            width = 6
          ),
          box(
            title = paste("MA Plot -", label),
            withSpinner(plotOutput(paste0("ma_plot_", method), height = 400)),
            width = 6
          )
        )
      }
    }))
  }) |> bindEvent(input$plot_button)

  # FUNCTIONS AND EXTENDED TASKS:
  # One-shot initializer (optional): build initial plots at startup
  filter_func <- function(counts, v) {
    "##############################################################################
    # n_min_val: minimum total count threshold for filtering
    # returns: updates reactiveValues v with filtered counts and title
    ##############################################################################"
    n_min_val <- v$n_min
    keep <- rs >= n_min_val
    v$counts_filt <- counts[keep,] # Keep the info of these genes
    v$keep <- keep
    v$title <- sprintf("After filtering (n_min = %d): %d genes remain",
                      n_min_val, nrow(v$counts_filt))
  }

  # ExtendedTask for building Mean-SD plot
  meanSdPlotTask <- ExtendedTask$new( 
    function(counts) {
    future_promise({
      meanSdPlot(as.matrix(counts), ranks = FALSE)
    })
  }) |> bind_task_button("update_n_min")      

  # ExtendedTask for building MA plots
  MAPlotTask <- ExtendedTask$new(
    function(counts) {
    future_promise({
      # Build MA for each column and bind
      ma_df <- bind_rows(lapply(colnames(counts), 
                        MA_comparison, 
                        data_matrix = counts))
      MA_Plot(ma_df)
    })
  }) |> bind_task_button("update_n_min")
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

normalization_func <- function(data_matrix, method, pseudocount=1, metadata.df=NULL){ 
  ##############################################################################
  # data_matrix: filtered raw counts
  # method: normalization method to apply
  ##############################################################################
  # Reorder metadata to match count matrix columns
  if (!is.null(metadata.df)) {
    print("Metadata provided, checking samples...")
  # Check if all samples in data_matrix are in metadata
  if (!all(colnames(data_matrix) %in% rownames(metadata.df))) {
    stop("Some samples in count matrix are missing from metadata")
  }
  # Reorder metadata to match count matrix order
  metadata.df <- metadata.df[colnames(data_matrix), , drop = FALSE]
    print("Metadata samples match count matrix columns.")
  }
  
  if(method == "tss"){
    norm_counts <- sweep(data_matrix, 2, colSums(data_matrix), FUN="/")

  }else if(method == "css"){
    min_s <- min(colSums(data_matrix))
    scaled_data <- min_s * data_matrix
    norm_counts <- sweep(scaled_data, 2, colSums(data_matrix), FUN="/")
  
  }else if(method == "clr"){
    # Add pseudo count
    data_pseudo <- data_matrix + pseudocount
    # Calculate geometric mean of each sample and divide read counts by this
    sample_gmeans <- apply(data_pseudo, 2, function(vec){
      exp(mean(log(vec)))
    })
    scaled_data <- sweep(data_pseudo, 2, sample_gmeans, FUN="/")
    norm_counts <- log2(scaled_data)
  
  }else if(method == "mor"){
    # Create dummy metadata
    dummy_colData <- DataFrame(
      condition = rep("A", ncol(data_matrix)),
      row.names = colnames(data_matrix)
    )
  
    dds <- DESeqDataSetFromMatrix(countData = data_matrix,
                                  colData = metadata.df,
                                  design = ~1)
    dds <- estimateSizeFactors(dds) 
    norm_counts <- as.data.frame(counts(dds, normalized=TRUE))

  }else if(method == "vst"){
    # Create dummy metadata
    dummy_colData <- DataFrame(
      condition = rep("A", ncol(data_matrix)),
      row.names = colnames(data_matrix)
    )

    dds <- DESeqDataSetFromMatrix(countData = data_matrix,
                                  colData = metadata.df,
                                  design = ~1)
    dds <- estimateSizeFactors(dds)
    norm_counts <- as.data.frame(assay(vst(dds, blind=TRUE)))
  
  }else if(method == "rlog"){
    # Create dummy metadata
    dummy_colData <- DataFrame(
      condition = rep("A", ncol(data_matrix)),
      row.names = colnames(data_matrix)
    )
    
    dds <- DESeqDataSetFromMatrix(countData = data_matrix,
                                  colData = metadata.df,
                                  design = ~1)
    dds <- estimateSizeFactors(dds)
    norm_counts <- as.data.frame(assay(rlog(dds, blind=TRUE))) # returns a data matrix
  }
  
  return(norm_counts)
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


# # Helper for MA Plot, prepares the data.
# ma_df <- reactive({
#   req(v$counts_filt)  # depends on committed filter
#   cm <- v$counts_filt
#   # Build MA for each column and bind
#   dplyr::bind_rows(lapply(colnames(cm), MA_comparison, data_matrix = cm))
# })

# # MA Plots        
# output$ma_plot <- renderPlot({
#   df <- ma_df()
#   MA_Plot(df)
# })

# nolint end


  # # Reactive event: normalization results change -> generate plots
  # observe({
  #   req(v$norm_results)
  #   for (method in v$norm_methods_committed) {
  #     local({
  #       m <- method
  #       output[[paste0("ms_plot_", m)]] <- renderPlot({
  #         res <- v$norm_results[[m]]
  #         req(!is.null(res))
  #         if (!is.null(res$error)) {
  #           plot.new()
  #           text(0.5, 0.5, paste("Error:", res$error))
  #           return()
  #         }
  #         ms <- vsn::meanSdPlot(as.matrix(res$norm_counts), ranks = FALSE)
  #         print(ms$meanSdPlot)
  #       })
  #       output[[paste0("ma_plot_", m)]] <- renderPlot({
  #         res <- v$norm_results[[m]]
  #         req(!is.null(res))
  #         if (!is.null(res$error)) {
  #           plot.new()
  #           text(0.5, 0.5, paste("Error:", res$error))
  #           return()
  #         }
  #         ma_df <- dplyr::bind_rows(lapply(
  #           colnames(res$norm_counts),
  #           MA_comparison,
  #           data_matrix = res$norm_counts
  #         ))
  #         MA_Plot(ma_df)
  #       })
  #     })
  #   }
  # })


  # observe({
  #   req(v$norm_methods_committed)
  #   for (method in v$norm_methods_committed) {
  #     local({
  #       # Require the normalization results to be ready
  #       req(v$norm_results[[method]])
  #       # Set local variables for each iteration
  #       m <- method
  #       norm_counts_snapshot <- v$norm_results[[m]]$norm_counts
  #       # Generate Mean-SD
  #       output[[paste0("ms_plot_", m)]] <- renderPlot({
  #         future_promise({
  #           meanSdPlot(as.matrix(norm_counts_snapshot), ranks = FALSE)
  #         })
  #       # Generate MA plot
  #       })
  #       output[[paste0("ma_plot_", m)]] <- renderPlot({
  #         future_promise({
  #           ma_df <- bind_rows(lapply(colnames(norm_counts_snapshot), MA_comparison, data_matrix = norm_counts_snapshot))
  #           MA_Plot(ma_df)
  #         })
  #       })
  #     })
  #   }
  # })
