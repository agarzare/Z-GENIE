###most recent working code, AGR 092724 @ 11:55pm EST (RShiny Version -- stand alone version)
# library(rsconnect)
# rsconnect::deployApp("~/Z-GENIE-Master")
# Sys.setenv(CXX11 = "g++")
# Sys.setenv(CXXFLAGS = "-std=c++11")
#install.packages("Rcpp", version = "1.0.6")
# Sys.setenv(CXXFLAGS = "-std=c++11")
# Sys.setenv(CXX11FLAGS = "-std=c++11")
# Function to check and install missing packages

options(install.packages.check.source = "no")

# Install BiocManager if it's not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Function to install a package only if it's not already installed
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE)
  }
}

# List of required packages
packages <- c("GenomicRanges", "BiocParallel", "Biostrings", "msa")

# Install the packages only if they are not already installed
sapply(packages, install_if_missing)

# Optionally install "msa" as a binary package if not present
if (!requireNamespace("msa", quietly = TRUE)) {
  BiocManager::install("msa", type = "binary", ask = FALSE)
}



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
library(processx)

# setwd("~/Z-GENIE-main/")
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

# Updated UI
ui <- dashboardPage(
  dashboardHeader(title = "Z-GENIE"),
  dashboardSidebar(sidebarMenu(
    menuItem("Home", tabName = "home", icon = icon("home")),
    menuItem("Run and Process", tabName = "run_process", icon = icon("cogs")),
    menuItem("Visualization", tabName = "visualization", icon = icon("chart-bar")),
    menuItem("MSA and Tree", tabName = "msa_tree", icon = icon("tree"))
  )),
  
  dashboardBody(
    useShinyjs(),
    tags$head(tags$style(
      HTML(
        "
        .content-wrapper { background-color: #f4f6f9 !important; }
        .box, .main-panel { border-radius: 10px; padding: 20px; background-color: white !important; }
        .main-header .logo { background-color: #0073b7 !important; }
        .main-header .navbar { background-color: #3c8dbc !important; }
        .main-sidebar .sidebar { background-color: #222d32 !important; color: white !important; }
        .main-sidebar .sidebar a { color: white !important; }
        .shiny-input-container { margin-bottom: 15px; }
        .shiny-bound-output { padding: 10px; background-color: white !important; border-radius: 10px; }
        .tab-pane { background-color: white !important; }
        .sidebar, .well { background-color: white !important; border: none; box-shadow: none; }
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
          p("This tool helps analyze Z-DNA genomic information using custom FASTA sequences or Z-Hunt output. You can either fetch sequences from NCBI or manually input data."),
          p("Adapted from: Ho, Pui S., et al. 'A computer aided thermodynamic approach for predicting the formation of Z‐DNA in naturally occurring sequences.' The EMBO journal 5.10 (1986): 2737-2744.")
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
          
          radioButtons("input_source", "Select Input Method:",
                       choices = list("Fetch FASTA from NCBI" = "fetch",
                                      "Upload FASTA File" = "upload",
                                      "Skip Fetch and Run Z-Hunt" = "manual"),
                       selected = "fetch"),
          
          # Conditional Panel for fetching FASTA
          conditionalPanel(
            condition = "input.input_source == 'fetch'",
            textInput("nucleotide_id", "Enter Nucleotide ID (e.g., \"U81553.1\")", value = ""),
            actionButton("fetch", "Fetch and Save FASTA"),
            downloadButton('download_fasta', 'Download Original FASTA'),
            textInput("fasta_path", "Enter Path to FASTA File", value = "Path/to/file/yourfile.fasta"),
            textInput("params", "Z-Hunt Parameters", value = "8, 6, 8"),
            actionButton("run", "Run Z-Hunt"),
            downloadButton('download_zscore', 'Download Z-SCORE Output'),
            numericInput("zscore_threshold", "Manual Z-Score Threshold", value = 600)
          ),
          
          # Conditional Panel for manually inputting paths
          conditionalPanel(
            condition = "input.input_source == 'manual'",
            fileInput("manual_fasta_path", "Upload Original FASTA File (with >)", accept = c(".fasta"), placeholder = "No file selected"),
            fileInput("manual_zscore_path", "Upload Z-SCORE File", accept = c(".Z-SCORE"), placeholder = "No file selected"),
            numericInput("manual_zscore_threshold", "Manual Z-Score Threshold", value = 600)
          ),
          
          # Conditional Panel for uploading a FASTA file and running Z-Hunt
          conditionalPanel(
            condition = "input.input_source == 'upload'",
            fileInput("upload_fasta", "Upload FASTA File", accept = c(".fasta")),
            textInput("params", "Z-Hunt Parameters", value = "8, 6, 8"),
            actionButton("run", "Run Z-Hunt"),
            downloadButton('download_zscore', 'Download Z-SCORE Output'),
            numericInput("zscore_threshold", "Manual Z-Score Threshold", value = 600)
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
      
      # MSA and Tree Tab
      tabItem(tabName = "msa_tree", fluidRow(
        column(width = 12,
               box(title = "Multiple Sequence Alignment", status = "info", solidHeader = TRUE, width = 12,
                   selectInput("alignment_method", "Alignment Method", choices = c("ClustalW", "ClustalOmega", "Muscle"), selected = "ClustalW"),
                   msaROutput("msa", width = "100%", height = "auto"))
        )),
        fluidRow(column(width = 12,
                        box(title = "Phylogenetic Tree", status = "info", solidHeader = TRUE, width = 12,
                            plotOutput("tree", height = "1500px")))
        )
      ))
  )
)


# Server
server <- function(input, output, session) {
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
                                   id = eval(expression(text = nucleotide_id)),
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
  
  # Function to run Z-Hunt
  run_zhunt <- function(input_file, params) {
    # Define the path to the zHunt binary
    Z_path <- normalizePath(file.path(getwd(), "zhunt/bin/zhunt"), mustWork = FALSE)
    
    # Check if the zHunt binary exists at the specified path
    if (!file.exists(Z_path)) {
      stop(paste("Error: The zHunt binary was not found at path:", Z_path))
    }
    
    # Ensure the zHunt binary has executable permissions
    system(paste("chmod +x", Z_path), ignore.stderr = TRUE)
    
    # Check if the zHunt binary is executable
    if (file.access(Z_path, mode = 1) != 0) {
      stop(paste("Error: The zHunt binary is not executable. Check the file permissions (chmod +x) and try again. Path:", Z_path))
    }
    
    # Prepare the modified FASTA file
    modified_fasta_file <- paste0(tools::file_path_sans_ext(input_file), "_mod.fasta")
    file.copy(input_file, modified_fasta_file, overwrite = TRUE)
    preprocess_fasta(modified_fasta_file)
    
    # Construct the command with the path to Z-Hunt and parameters
    command <- paste(Z_path, paste(params, collapse = " "), modified_fasta_file)
    
    # Debugging output: Print the constructed command to ensure it's correct
    cat("Executing command:", command, "\n")
    
    # Execute the Z-Hunt binary using processx::run instead of system for better error handling
    tryCatch({
      result <- processx::run(
        command = Z_path,
        args = c(params, modified_fasta_file),
        error_on_status = FALSE,
        echo_cmd = TRUE,
        echo = TRUE
      )
      
      # Check if there's an error in stderr output
      if (result$stderr != "") {
        return(paste("Error occurred during Z-Hunt execution:\n", result$stderr))
      }
      
      # If successful, return the result from stdout
      return(result$stdout)
      
    }, error = function(e) {
      return(paste("Error running Z-Hunt:", e$message))
    })
  }
  
  # Fetch and Save FASTA File
  observeEvent(input$fetch, {
    req(input$nucleotide_id)
    showNotification("Fetching FASTA file...", type = "message")
    original_fasta_file <- fetch_and_save_fasta(input$nucleotide_id)
    updateTextInput(session, "fasta_path", value = normalizePath(original_fasta_file))
    showNotification("FASTA file fetched successfully!", type = "message")
    
    # Save the path of the fetched FASTA file for download
    output$download_fasta <- downloadHandler(
      filename = function() {
        "fasta.fasta"
      },
      content = function(file) {
        file.copy(original_fasta_file, file)
      },
      contentType = "text/plain"
    )
  })
  
  
  
  # Run Z-Hunt
  observeEvent(input$run, {
    # Check the selected input source
    if (input$input_source == "upload") {
      req(input$upload_fasta)
      original_fasta_file <- input$upload_fasta$datapath
    } else if (input$input_source == "fetch") {
      req(input$fasta_path)
      original_fasta_file <- normalizePath(input$fasta_path)
    } else {
      # For manual input, we do not run Z-Hunt
      output$output <- renderText({
        "Z-Hunt run is not applicable for 'manual' input method."
      })
      return()
    }
    
    # Split the input string into a vector and validate parameters
    params <- as.numeric(unlist(strsplit(input$params, ",")))
    
    # Check if params has exactly 3 elements
    if (length(params) != 3) {
      output$output <- renderText({
        "Error: Please ensure 'params' contains exactly three numbers separated by commas (e.g., '16, 6, 16')"
      })
      return()
    }
    
    # Run the Z-Hunt command using the updated run_zhunt function
    command_output <- tryCatch({
      run_zhunt(original_fasta_file, params)
    }, error = function(e) {
      paste("Error running Z-Hunt:", e$message)
    })
    
    # Render the output in the Shiny UI
    output$output <- renderText({
      paste("Z-Hunt run complete for file:", original_fasta_file,
            "\nParameters used:", paste(params, collapse = " "),
            "\n\nOutput:\n", paste(command_output, collapse = "\n"))
    })
    
    # Save the Z-SCORE output to a temporary file
    output_file_path <- paste0(tools::file_path_sans_ext(original_fasta_file), "_mod.fasta.Z-SCORE")
    
    
    # Add the downloadHandler for the output
    output$download_zscore <- downloadHandler(
      filename = function() {
        paste0("fasta.Z.SCORE")
      },
      content = function(file) {
        file.copy(output_file_path, file)
      },
      contentType = "text/plain"
    )
  })
  
  
  
  # Process Z-SCORE File
  observeEvent(input$process, {
    req(input$input_source)  # Ensure input source is selected
    
    # Handling different input sources
    if (input$input_source == "manual") {
      req(input$manual_fasta_path, input$manual_zscore_path)
      original_fasta_file <- input$manual_fasta_path$datapath
      zscore_file <- input$manual_zscore_path$datapath
      threshold <- input$manual_zscore_threshold
      
      # Remove the first line from the Z-SCORE file if necessary
      remove_first_line(zscore_file)
      
    } else if (input$input_source == "upload") {
      req(input$upload_fasta)  # Ensure the FASTA file is uploaded
      original_fasta_file <- input$upload_fasta$datapath
      zscore_file <- paste0(tools::file_path_sans_ext(original_fasta_file), "_mod.fasta.Z-SCORE")
      threshold <- input$zscore_threshold
      
      # Remove the first line from the Z-SCORE file if necessary
      remove_first_line(zscore_file)
      
    } else if (input$input_source == "fetch") {
      req(input$fasta_path)
      original_fasta_file <- normalizePath(input$fasta_path)
      zscore_file <- paste0(tools::file_path_sans_ext(original_fasta_file), "_mod.fasta.Z-SCORE")
      threshold <- input$zscore_threshold
      
      # Remove the first line from the Z-SCORE file if necessary
      remove_first_line(zscore_file)
    }
    
    # Check if the Z-SCORE file exists
    if (file.exists(zscore_file)) {
      showNotification("Processing Z-SCORE file...", type = "message")
      
      
      # Continue with the processing steps
      processed_data <- fread(zscore_file, col.names = c("V1", "V2", "zscore", "conformation"))
      processed_data$Position <- 1:nrow(processed_data)
      processed_data$max_Z_Score <- processed_data$zscore
      processed_data$OligoConformation <- processed_data$conformation
      
      threshold <- as.numeric(threshold)
      processed_data <- processed_data[order(processed_data$max_Z_Score, decreasing = TRUE), ]
      
      # Ensure to use `which` instead of `grep` to find rows that match the threshold criteria
      end <- tail(which(processed_data$max_Z_Score >= threshold), n = 1)
      
      if (length(end) == 0) {
        showNotification("No scores meet the threshold criteria.", type = "error")
        return()
      }
      
      processed_data <- processed_data[1:end, ]
      processed_data <- processed_data[order(processed_data$Position, decreasing = FALSE), ]
      
      
      # params <- as.numeric(unlist(strsplit(input$params, " ")))
      # middle_value <- params[2]
      # limit_value <- 2 * middle_value
      # Correct the params handling
      params <- as.numeric(trimws(unlist(strsplit(input$params, ","))))
      
      # Now params should be a numeric vector: c(16, 6, 16)
      # Extract middle_value and calculate limit_value
      if (length(params) == 3) {
        middle_value <- params[2]
        limit_value <- 2 * middle_value
      } else {
        stop("Error: 'params' should contain exactly three numbers separated by commas (e.g., '16, 6, 16')")
      }
      
      
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
      
      ###include p53 motif
      # Function to find 5′-RRRCWWGYYY-3′ motif in a sequence
      find_motif_pattern <- function(sequence) {
        # Define the pattern for RRRCWWGYYY (R = A|G, W = A|T, Y = C|T)
        pattern <- "(A|G)(A|G)(A|G)(C)(A|T)(A|T)(G)(C|T)(C|T)(C|T)"
        
        # Find matches of the updated pattern in the DNA sequence
        matches <- gregexpr(pattern, sequence, perl = TRUE)
        
        # Handle cases where no matches are found
        if (length(matches[[1]]) == 1 && matches[[1]][1] == -1) {
          return(data.frame(
            start = integer(0),
            end = integer(0),
            sequence = character(0)
          ))
        }
        
        # Extract matching sequences
        matching_sequences <- regmatches(sequence, matches)
        start_positions <- matches[[1]]
        end_positions <- start_positions + attr(matches[[1]], "match.length") - 1
        
        # Create a data frame to hold the results
        result <- data.frame(start = start_positions,
                             end = end_positions,
                             sequence = matching_sequences[[1]])
        return(result)
      }
      
      # Function to analyze a dataframe of sequences for the motif pattern
      analyze_sequences2 <- function(df) {
        df$p53motif <- NA  # Column to store TRUE/FALSE if the motif is found
        df$p53motifeq <- NA  # Column to store the matching sequences
        
        # Loop through each sequence in the dataframe
        for (i in 1:nrow(df)) {
          sequence <- df$sequence[i]
          match_result <- find_motif_pattern(sequence)
          
          # Check if any matches were found
          if (nrow(match_result) > 0) {
            df$p53motif[i] <- "TRUE"
            df$p53motifeq[i] <- paste(match_result$sequence, collapse = ",")
          } else {
            df$p53motif[i] <- "FALSE"
          }
        }
        
        return(df)
      }
      
      
      # After result_df is created and populated with data
      result_df <- analyze_sequences(df)
      result_df <- analyze_sequences2(result_df)
      
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
  
  observeEvent(input$file, {
    req(input$file)
    df <- read.csv(input$file$datapath, header = TRUE) %>%
      mutate_all(~ if (is.character(.))
        as.factor(.)
        else
          as.numeric(.))
    updateSelectInput(session, "x", choices = names(df), selected = "start")
    updateSelectInput(session, "y", choices = names(df), selected = "log10_ZScore")
    updateSelectInput(session,
                      "color",
                      choices = names(df),
                      selected = "ISS")
    data(df)
  })
  
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