###most recent working code, AGR 090824 @ 4:54pm EST (Stand-alone Version)
# library(rsconnect)
# rsconnect::deployApp("~/Z-GENIE-Master")
# Sys.setenv(CXX11 = "g++")
# Sys.setenv(CXXFLAGS = "-std=c++11")
#install.packages("Rcpp", version = "1.0.6")
# Sys.setenv(CXXFLAGS = "-std=c++11")
# Sys.setenv(CXX11FLAGS = "-std=c++11")
# Function to check and install missing packages

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c("GenomicRanges", "BiocParallel", "Biostrings"))
BiocManager::install("msa", type = "binary")  # Optional, if you encounter issues with `msa`

BiocManager::install("msa", type = "binary")


library(shiny)
library(shinydashboard)
library(shinyjs)
library(plotly)
library(DT)
library(dplyr)
library(magrittr)
library(Biostrings)
library(msaR)
library(seqinr)
library(msa)
library(ape)
library(rentrez)
library(parallel)
library(doParallel)
library(reticulate)
library(ggplot2)
library(hrbrthemes)
library(tidyr)
library(viridis)
library(data.table)
library(GenomicRanges)
library(stringr)
library(BiocParallel)
library(renv)
# setwd("~/Z-GENIE-Master")
# renv::init()
# renv::snapshot()
# renv::activate()
# renv::restore()
# Set R to use CRAN and Bioconductor repositories for package installation
# options(repos = BiocManager::repositories())
# BiocManager::install("msa", type = "binary")

# Now when you install a package, R will look in both CRAN and Bioconductor
# install.packages("msa")

ui <- dashboardPage(
  dashboardHeader(title = "Z-GENIE"),
  dashboardSidebar(sidebarMenu(
    menuItem("Home", tabName = "home", icon = icon("home")),
    menuItem("Run and Process", tabName = "run_process", icon = icon("cogs")),
    menuItem(
      "Visualization",
      tabName = "visualization",
      icon = icon("chart-bar")
    ),
    menuItem("MSA and Tree", tabName = "msa_tree", icon = icon("tree"))
  )),
  
  dashboardBody(
    useShinyjs(),
    tags$head(tags$style(
      HTML(
        "
        .content-wrapper {
          background-color: #f4f6f9 !important;
        }
        .box, .main-panel {
          border-radius: 10px;
          padding: 20px;
          background-color: white !important; /* Force white background for all box and panel elements */
        }
        .main-header .logo {
          background-color: #0073b7 !important;
        }
        .main-header .navbar {
          background-color: #3c8dbc !important;
        }
        .main-sidebar .sidebar {
          background-color: #222d32 !important;
          color: white !important;
        }
        .main-sidebar .sidebar a {
          color: white !important;
        }
        .shiny-input-container {
          margin-bottom: 15px;
        }
        .shiny-bound-output {
          padding: 10px;
          background-color: white !important; /* Ensuring white background for output panels */
          border-radius: 10px;
        }
        .tab-pane {
          background-color: white !important; /* Making sure tab panes have a white background */
        }
        .sidebar, .well {
          background-color: white !important; /* Apply white to sidebars and wells */
          border: none;
          box-shadow: none;
        }
      "
      )
    )),
    
    tabItems(
      # Home Tab
      tabItem(tabName = "home", fluidRow(
        box(
          title = "Welcome to Z-GENIE",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          h4("Z-GENIE (Z-DNA GENomic Information Extractor)"),
          p(
            "This tool helps analyze Z-DNA genomic information using custom FASTA sequences or Z-Hunt output. You can either fetch sequences from NCBI or manually input data."
          )
        )
      )),
      
      # Run and Process Tab
      tabItem(tabName = "run_process", fluidRow(
        box(
          title = "Step 1: Run and Process",
          status = "info",
          solidHeader = TRUE,
          width = 12,
          collapsible = TRUE,
          
          checkboxInput("skip_steps", "Skip Fetch and Run Z-Hunt", value = FALSE),
          
          conditionalPanel(
            condition = "input.skip_steps == false",
            textInput("nucleotide_id", "Enter Nucleotide ID (e.g., U81553.1)", value = ""),
            actionButton("fetch", "Fetch and Save FASTA"),
            textInput("fasta_path", "Enter Path to FASTA File", value = "Path/to/file/yourfile.fasta"),
            textInput("params", "Z-Hunt Parameters", value = "8 6 8"),
            actionButton("run", "Run Z-Hunt"),
            numericInput("zscore_threshold", "Manual Z-Score Threshold", value = 600)
          ),
          
          conditionalPanel(
            condition = "input.skip_steps == true",
            textInput(
              "manual_fasta_path",
              "Enter Path to Original FASTA File (with >)",
              value = "Path/to/original.fasta"
            ),
            textInput("manual_zscore_path", "Enter Path to Z-SCORE File", value = "Path/to/zscore.Z-SCORE"),
            numericInput("manual_zscore_threshold", "Manual Z-Score Threshold", value = 600)
          ),
          
          actionButton("process", "Process Z-SCORE File"),
          textInput("download_name", "Enter Download File Name", value = "Filtered_Z_GENIE_output"),
          downloadButton('download', 'Download Processed Data'),
          DTOutput("data_table"),
          verbatimTextOutput("fetch_output"),
          verbatimTextOutput("output"),
          verbatimTextOutput("command_output"),
          verbatimTextOutput("final_output")
        )
      )),
      
      # Visualization Tab
      tabItem(tabName = "visualization", fluidRow(
        box(
          title = "ZFS Visualization",
          status = "info",
          solidHeader = TRUE,
          width = 12,
          collapsible = TRUE,
          sidebarLayout(
            sidebarPanel(
              fileInput("file", "Upload CSV file", buttonLabel = "Browse..."),
              hr(),
              h4("Filter DataTable"),
              hr(),
              h4("Plotly Configuration"),
              selectInput("x", "X-axis:", choices = NULL, selected = "start"),
              selectInput("y", "Y-axis:", choices = NULL, selected = "log10_ZScore"),
              selectInput("color", "Color by:", choices = NULL, selected = "ISS"),
              DTOutput("table", width = "auto"),
              downloadButton("download2", "Download Filtered Data")
            ),
            mainPanel(plotlyOutput("plot"))
          )
        )
      )),
      # Modified MSA and Tree Tab (MSA on top, Tree at the bottom)
      tabItem(tabName = "msa_tree", fluidRow(column(
        width = 12,
        # Full width
        box(
          title = "Multiple Sequence Alignment",
          status = "info",
          solidHeader = TRUE,
          width = 12,
          selectInput(
            "alignment_method",
            "Alignment Method",
            choices = c("ClustalW", "ClustalOmega", "Muscle"),
            selected = "ClustalW"
          ),
          msaROutput("msa", width = "100%", height = "auto")  # Set appropriate height for the MSA output
        )
      )), fluidRow(column(
        width = 12,
        # Full width
        box(
          title = "Phylogenetic Tree",
          status = "info",
          solidHeader = TRUE,
          width = 12,
          plotOutput("tree", height = "1500px")  # Set appropriate height for the tree plot
        )
      )))
      
    )
  )
)






# Server
server <- function(input, output, session) {
  # Function to check if a command is available and install if missing
  check_command <- function(command, install_command = NULL) {
    if (Sys.which(command) == "") {
      if (!is.null(install_command)) {
        showNotification(paste(command, "is not installed. Installing..."),
                         type = "warning")
        system(install_command)
      }
      if (Sys.which(command) == "") {
        stop(paste(
          command,
          "could not be installed. Please install it manually."
        ))
      }
    }
  }
  
  # Check for 'git' and 'gcc' availability, install if missing
  check_command("git", "sudo apt-get install git -y")
  check_command("gcc", "sudo apt-get install gcc -y")
  
  #Clone the Z-Hunt-III repository from GitHub if not already done
  if (!file.exists("zhunt")) {
    system("git clone https://github.com/Ho-Lab-Colostate/zhunt.git")
  }
  
  # Helper function to preprocess the FASTA file (remove >)
  preprocess_fasta <- function(file_path) {
    fasta_lines <- readLines(file_path)
    if (startsWith(fasta_lines[1], ">")) {
      fasta_lines <- fasta_lines[-1]
      writeLines(fasta_lines, file_path)
    }
  }
  
  # Fetch FASTA from NCBI
  fetch_and_save_fasta <- function(nucleotide_id) {
    fasta <- rentrez::entrez_fetch(db = "nucleotide",
                                   id = nucleotide_id,
                                   rettype = "fasta")
    file_name <- paste0(nucleotide_id, ".fasta")
    write(fasta, file = file_name)
    return(file_name)
  }
  
  # Remove the first line from a .fasta.Z-SCORE file if it starts with "/"
  remove_first_line <- function(file_path) {
    lines <- readLines(file_path)
    if (startsWith(lines[1], "/")) {
      lines <- lines[-1]
    }
    writeLines(lines, file_path)
  }
  
  # Run Z-Hunt
  run_zhunt <- function(input_file, params) {
    zhunt_path <- normalizePath("zhunt")
    Z_path <- file.path(zhunt_path, "bin", "zhunt")
    
    modified_fasta_file <- paste0(tools::file_path_sans_ext(input_file), "_mod.fasta")
    file.copy(input_file, modified_fasta_file, overwrite = TRUE)
    preprocess_fasta(modified_fasta_file)
    
    command <- paste(Z_path, params, modified_fasta_file)
    output <- system(command, intern = TRUE, ignore.stderr = FALSE)
    output
  }
  
  # Observe events and process based on input
  observeEvent(input$fetch, {
    req(input$nucleotide_id)
    showNotification("Fetching FASTA file...", type = "message")
    original_fasta_file <- fetch_and_save_fasta(input$nucleotide_id)
    updateTextInput(session, "fasta_path", value = normalizePath(original_fasta_file))
    showNotification("FASTA file fetched successfully!", type = "message")
  })
  
  observeEvent(input$run, {
    req(input$fasta_path)
    original_fasta_file <- normalizePath(input$fasta_path)
    params <- input$params
    command_output <- run_zhunt(original_fasta_file, params)
    
    output$output <- renderText({
      paste("Z-Hunt run complete for file:",
            original_fasta_file,
            "with parameters:",
            params)
    })
    
    output$command_output <- renderText({
      paste("Command output:\n",
            paste(command_output, collapse = "\n"))
    })
  })
  
  observeEvent(input$process, {
    if (input$skip_steps) {
      # Manual input
      original_fasta_file <- input$manual_fasta_path
      zscore_file <- input$manual_zscore_path
      threshold <- input$manual_zscore_threshold
    } else {
      # From fetched data
      original_fasta_file <- normalizePath(input$fasta_path)
      zscore_file <- paste0(tools::file_path_sans_ext(original_fasta_file),
                            "_mod.fasta.Z-SCORE")
      threshold <- input$zscore_threshold
    }
    
    remove_first_line(zscore_file)
    
    if (file.exists(zscore_file)) {
      showNotification("Processing Z-SCORE file...", type = "message") # Process the Z-SCORE file and prepare data
      processed_data <- fread(zscore_file,
                              col.names = c("V1", "V2", "zscore", "conformation"))
      processed_data$Position <- 1:nrow(processed_data)
      processed_data$max_Z_Score <- processed_data$zscore
      processed_data$OligoConformation <- processed_data$conformation
      
      threshold <- as.numeric(threshold)
      processed_data <- processed_data[order(processed_data$max_Z_Score, decreasing = TRUE), ]
      end <- tail(grep(processed_data$max_Z_Score >= threshold, pattern = TRUE),
                  n = 1)
      processed_data <- processed_data[1:end, ]
      processed_data <- processed_data[order(processed_data$Position, decreasing = FALSE), ]
      
      params <- as.numeric(unlist(strsplit(input$params, " ")))
      middle_value <- params[2]
      limit_value <- 2 * middle_value
      
      tryCatch({
        c1 <- r_to_py(processed_data$Position)
        py$c1 <- c1
      }, error = function(e) {
        print(e)
      })
      
      py_run_string(
        sprintf(
          "
import sys
tst = c1
il = []
ol = []
for k, v in enumerate(tst):
    if k > 0:
        if abs(tst[k] - tst[k-1]) < %d:
            if tst[k-1] not in il:
                il.append(tst[k-1])
            if tst[k] not in il:
                il.append(tst[k])
        else:
            ol.append(list(il))
            il = []
ol.append(list(il))
c1_ol = ol
            ",
          limit_value
        )
      )
      
      c1_ol <- py$c1_ol
      
      c1_empty <- grep("list()", c1_ol)
      c1_list <- c1_ol[-c(c1_empty)]
      
      if (length(c1_ol) == 0) {
        output$final_output <- renderText({
          "c1_ol is empty."
        })
        return()
      }
      
      if (length(c1_list) == 0) {
        output$final_output <- renderText({
          "c1_list is empty."
        })
        return()
      }
      
      non_Mountain_c1_list <- list()
      
      for (i in seq_along(c1_list)) {
        lo <- c1_list[[i]]
        non_Mountain_c1_list <- c(lo, non_Mountain_c1_list)
      }
      
      non_Mountain_c1_list <- unlist(non_Mountain_c1_list)
      non_Mountain_c1_list <- as.data.frame(non_Mountain_c1_list)
      non_Mountain_c1_list <- processed_data[-which(processed_data$Position %in% non_Mountain_c1_list$non_Mountain_c1_list), ]
      non_Mountain_c1_list <- as.list(non_Mountain_c1_list$Position)
      c1_list <- c(non_Mountain_c1_list, c1_list)
      
      # Use the original FASTA file for reading sequences
      processed_fasta <- read.fasta(original_fasta_file)
      dna <- toString(processed_fasta[[names(processed_fasta)]])
      dna <- gsub(", ", "", dna)
      dna <- toupper(dna)
      dna <- as.character(dna)
      DNA <- DNAString(dna)
      c1_fasta <- s2c(dna)
      
      if (length(c1_fasta) == 0) {
        output$final_output <- renderText({
          "c1_fasta is empty."
        })
        return()
      }
      
      Mc1 <- data.frame(
        "start",
        "end",
        "range",
        "sequence",
        "oligo_length",
        "GC_of_sequence",
        "max_Z_Score"
      )
      names(Mc1) <- c(
        "start",
        "end",
        "range",
        "sequence",
        "oligo_length",
        "GC_of_sequence",
        "max_Z_Score"
      )
      
      system.time({
        for (i in seq_along(c1_list)) {
          start <- head(c1_list[[i]], n = 1)
          last <- tail(c1_list[[i]], n = 1)
          position <- grep(processed_data$Position,
                           pattern = paste0("\\b", eval(parse(text = last)), "\\b"))
          length <- nchar(processed_data$OligoConformation[[position]])
          end <- last + length - 1
          range <- paste0(start, ":", end)
          sequence <- toString(c1_fasta[eval(parse(text = range))])
          sequence <- gsub(", ", "", sequence)
          sequence <- toupper(sequence)
          sequence <- as.character(sequence)
          oligo_length <- nchar(sequence)
          Sequence <- s2c(sequence)
          GC_of_sequence <- GC(Sequence) * 100
          GC_of_sequence <- as.numeric(formatC(
            GC_of_sequence,
            format = "f",
            digits = 2
          ))
          zscore_position_start <- grep(processed_data$Position, pattern = start)
          zscore_position_end <- position
          max_Z_Score <- max(processed_data$max_Z_Score[zscore_position_start:zscore_position_end])
          Mc1[i, ] <- c(start,
                        end,
                        range,
                        sequence,
                        oligo_length,
                        GC_of_sequence,
                        max_Z_Score)
        }
        Masterlist_Chr1_ZScoresgreaterthan600 <- Mc1
      })
      
      Masterlist_Chr1_ZScoresgreaterthan600$start <- as.numeric(Masterlist_Chr1_ZScoresgreaterthan600$start)
      Masterlist_Chr1_ZScoresgreaterthan600$end <- as.numeric(Masterlist_Chr1_ZScoresgreaterthan600$end)
      Masterlist_Chr1_ZScoresgreaterthan600$oligo_length <- as.numeric(Masterlist_Chr1_ZScoresgreaterthan600$oligo_length)
      Masterlist_Chr1_ZScoresgreaterthan600$GC_of_sequence <- as.numeric(Masterlist_Chr1_ZScoresgreaterthan600$GC_of_sequence)
      Masterlist_Chr1_ZScoresgreaterthan600$max_Z_Score <- as.numeric(Masterlist_Chr1_ZScoresgreaterthan600$max_Z_Score)
      
      df <- Masterlist_Chr1_ZScoresgreaterthan600
      df$log10_ZScore <- log(df$max_Z_Score, base = 10)
      
      calculate_at_content <- function(sequence) {
        at_count <- sum(str_count(sequence, "AT") + str_count(sequence, "TA"))
        total_count <- nchar(sequence) - 1
        at_content <- at_count / total_count
        return(at_content)
      }
      
      calculate_gc_content <- function(sequence) {
        gc_count <- sum(str_count(sequence, "GC") + str_count(sequence, "CG"))
        total_count <- nchar(sequence) - 1
        gc_content <- gc_count / total_count
        return(gc_content)
      }
      
      calculate_atgc_content <- function(sequence) {
        at_count <- sum(str_count(sequence, "A") + str_count(sequence, "T"))
        gc_count <- sum(str_count(sequence, "G") + str_count(sequence, "C"))
        atgc_content <- at_count / gc_count
        return(atgc_content)
      }
      
      df <- df %>%
        mutate(
          at_content = sapply(sequence, calculate_at_content),
          gc_content = sapply(sequence, calculate_gc_content),
          atgc_content = sapply(sequence, calculate_atgc_content)
        )
      
      df$at_content <- df$at_content * 100
      df$gc_content <- df$gc_content * 100
      
      find_purine_pattern <- function(sequence) {
        pattern <- "(A|G)(A|G)(CG)+?(C|T)(C|T)"
        matches <- gregexpr(pattern, sequence, perl = TRUE)
        if (length(matches[[1]]) == 1 && matches[[1]][1] == -1) {
          return(data.frame(
            start = integer(0),
            end = integer(0),
            sequence = character(0)
          ))
        }
        matching_sequences <- regmatches(sequence, matches)
        start_positions <- matches[[1]]
        end_positions <- start_positions + attr(matches[[1]], "match.length") - 1
        result <- data.frame(start = start_positions,
                             end = end_positions,
                             sequence = matching_sequences[[1]])
        return(result)
      }
      
      analyze_sequences <- function(df) {
        df$ISS <- NA
        df$ISSeq <- NA
        for (i in 1:nrow(df)) {
          sequence <- df$sequence[i]
          match_result <- find_purine_pattern(sequence)
          if (nrow(match_result) > 0) {
            df$ISS[i] <- "TRUE"
            df$ISSeq[i] <- paste(match_result$sequence, collapse = ",")
          } else {
            df$ISS[i] <- "FALSE"
          }
        }
        return(df)
      }
      
      # After result_df is created and populated with data
      result_df <- analyze_sequences(df)
      
      # Truncate all numerical values in the dataframe to two decimal places
      result_df <- result_df %>%
        mutate_if(is.numeric, ~ round(., 2))
      
      find_palindromes <- function(sequence) {
        dna <- DNAString(sequence)
        palindromes <- findPalindromes(dna)
        if (length(palindromes) == 0) {
          return(NULL)
        } else {
          palindrome_sequences <- as.character(Views(dna, palindromes))
          return(paste(palindrome_sequences, collapse = ", "))
        }
      }
      
      result_df$palindromic_sequence <- sapply(result_df$sequence, find_palindromes)
      result_df$palindromic_sequence <- as.character(result_df$palindromic_sequence)
      
      output$final_output <- renderText({
        paste("Processed and finalized dataset with",
              nrow(result_df),
              "rows.")
      })
      
      Zp <- (sum(result_df$oligo_length) / nchar(dna)) * 100
      
      output$final_output <- renderText({
        paste("Z-Potentiality is ", Zp, "%")
      })
      
      output$data_table <- renderDT({
        datatable(result_df,
                  options = list(pageLength = 10, searchHighlight = TRUE))
      })
      
      output$download <- downloadHandler(
        filename = function() {
          paste0(input$download_name, ".csv")  # Use the input from the user as the filename
        },
        content = function(fname) {
          write.csv(result_df,
                    fname,
                    row.names = F,
                    col.names = T)
        }
      )
    } else {
      showNotification("The .fasta.Z-SCORE file was not found.", type = "error")
    }
  })
  
  
  # Visualization plot
  data <- reactiveVal(NULL)
  
  # observeEvent(input$file, {
  #   req(input$file)
  #   df <- read.csv(input$file$datapath, header = TRUE) %>%
  #     mutate_all(~ if (is.character(.))
  #       as.factor(.)
  #       else
  #         as.numeric(.))
  #   updateSelectInput(session, "x", choices = names(df), selected = "start")
  #   updateSelectInput(session, "y", choices = names(df), selected = "log10_ZScore")
  #   updateSelectInput(session,
  #                     "color",
  #                     choices = names(df),
  #                     selected = "ISS")
  #   data(df)
  # })

  observeEvent(input$file, {
    req(input$file)
    
    # Read CSV and apply mutation
    df <- read.csv(input$file$datapath, header = TRUE) %>%
      mutate_all(~ {
        if (is.character(.)) {
          as.factor(.)
        } else {
          as.numeric(.)
        }
      })
    
    # Update select inputs based on the data
    updateSelectInput(session, "x", choices = names(df), selected = "start")
    updateSelectInput(session, "y", choices = names(df), selected = "log10_ZScore")
    updateSelectInput(session, "color", choices = names(df), selected = "ISS")
    
    # Set reactive data value
    data(df)
  })
  
  output$table <- renderDT({
    req(data())
    datatable(data(), filter = 'top', selection = 'none')
  })
  
  output$plot <- renderPlotly({
    req(data())
    
    # Get the selected values
    x_var <- input$x
    y_var <- input$y
    color_var <- input$color
    
    plot_ly(
      data = data(),
      x = ~ get(x_var),
      y = ~ get(y_var),
      color = ~ get(color_var),
      type = 'scatter',
      mode = 'markers',
      text = ~ paste(
        "Sequence: ",
        data()$sequence,
        "<br>Start: ",
        data()$start,
        "<br>End: ",
        data()$end,
        "<br>Oligo Length: ",
        data()$oligo_length,
        "<br>GC Content: ",
        data()$GC_of_sequence,
        "<br>AT Content: ",
        data()$at_content,
        "<br>Max Z-Score: ",
        data()$max_Z_Score,
        "<br>ISS: ",
        data()$ISS,
        "<br>ISSeq: ",
        data()$ISSeq,
        "<br>Palindromic Sequence: ",
        data()$palindromic_sequence
      )
    ) %>% layout(
      dragmode = 'select',
      xaxis = list(title = x_var, range = c(min(data(
        
      )[[x_var]]), max(data(
        
      )[[x_var]]))),
      yaxis = list(title = y_var, range = c(min(data(
        
      )[[y_var]]), max(data(
        
      )[[y_var]])))
    )
  })
  
  observe({
    req(data())
    if (!is.null(input$table_rows_all)) {
      filtered_data <- data()[input$table_rows_all, ]
      output$plot <- renderPlotly({
        req(filtered_data)
        
        # Get the selected values
        x_var <- input$x
        y_var <- input$y
        color_var <- input$color
        
        plot_ly(
          data = filtered_data,
          x = ~ get(x_var),
          y = ~ get(y_var),
          color = ~ get(color_var),
          type = 'scatter',
          mode = 'markers',
          text = ~ paste(
            "Sequence: ",
            filtered_data$sequence,
            "<br>Start: ",
            filtered_data$start,
            "<br>End: ",
            filtered_data$end,
            "<br>Oligo Length: ",
            filtered_data$oligo_length,
            "<br>GC Content: ",
            filtered_data$GC_of_sequence,
            "<br>AT Content: ",
            filtered_data$at_content,
            "<br>Max Z-Score: ",
            filtered_data$max_Z_Score,
            "<br>ISS: ",
            filtered_data$ISS,
            "<br>ISSeq: ",
            filtered_data$ISSeq,
            "<br>Palindromic Sequence: ",
            filtered_data$palindromic_sequence
          )
        )  %>% layout(
          dragmode = 'select',
          xaxis = list(title = x_var, range = c(
            min(filtered_data[[x_var]]), max(filtered_data[[x_var]])
          )),
          yaxis = list(title = y_var, range = c(
            min(filtered_data[[y_var]]), max(filtered_data[[y_var]])
          ))
        )
      })
    }
  })
  
  observe({
    req(data())
    if (!is.null(input$table_rows_all)) {
      filtered_data <- data()[input$table_rows_all, ]
      req(filtered_data)
      output$download2 <- downloadHandler(
        filename = function() {
          paste0(input$download_name, "_Filtered.csv")  # Use the input from the user as the filename for the filtered data
        },
        content = function(fname) {
          write.csv(filtered_data,
                    fname,
                    row.names = F,
                    col.names = T)
        }
      )
      output$msa <- renderMsaR({
        req(filtered_data)
        cdr3aa_reactive <- reactive({
          req(filtered_data)
          t <- as.data.frame(filtered_data[, 4])
          colnames(t) <- NULL
          t[, 1] <- as.data.frame(toupper(t[, 1]))
          
          # Use the selected alignment method
          alignment_method <- input$alignment_method
          
          aa_msa <- msa(t[, 1],
                        method = alignment_method,
                        type = "protein",
                        order = "input")
          AA_msa <- msaConvert(aa_msa, type = "seqinr::alignment")
          aa_msa <- AA_msa
          aa_msa$nam <- t[, 1]
          AA_msa <- as.data.frame(AA_msa[["seq"]])
          colnames(AA_msa) <- NULL
          
          AA_msa <- as.data.frame(t(AA_msa))
          colnames(AA_msa) <- AA_msa[1, ]
          write.fasta(
            sequences = AA_msa,
            names = names(AA_msa),
            file.out = "filtered_sequences.fasta"
          )
          
          output$tree <- renderPlot({
            req(aa_msa)
            d <- dist.alignment(aa_msa, "identity")
            as.matrix(d)
            hemoTree <- njs(d)
            tree <- plot(hemoTree, main = "Phylogenetic Tree of Filtered Sequences")
            
          })
          proteins <- ape::read.FASTA("filtered_sequences.fasta", type =
                                        "AA")
          return(proteins)
        })
        req(cdr3aa_reactive())
        msaR(
          cdr3aa_reactive(),
          colorscheme = "clustal",
          conservation = T,
          labelNameLength = 300,
          overviewboxWidth = "fixed",
          overviewboxHeight = "fixed",
          alignmentHeight = 500,
          rowheight = 30,
          height = 500,
          width = 1500
        )
      })
      
    }
  })
  
}


shinyApp(ui, server)
