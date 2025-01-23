###most recent working code, AGR 012225 @ 12:37 AM EST (RShiny Version -- stand alone version)
# library(rsconnect)
# rsconnect::deployApp("~/Z-GENIE-Master")
# Sys.setenv(CXX11 = "g++")
# Sys.setenv(CXXFLAGS = "-std=c++11")
# install.packages("Rcpp", version = "1.0.6")
# Sys.setenv(CXXFLAGS = "-std=c++11")
# Sys.setenv(CXX11FLAGS = "-std=c++11")

options(install.packages.check.source = "no")


# OPTIONAL: Increase max upload size if large files are needed
options(shiny.maxRequestSize = 1024 * 1024^2) # 1 GB
# options(shiny.maxRequestSize = 100 * 1024^2)  # 100 MB, etc.

# Install BiocManager if not present
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Function to install BiocManager packages if not already installed
install_if_missing <- function(pkg){
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE)
  }
}

# List of BiocManager required packages
packages <- c("GenomicRanges","BiocParallel","Biostrings","msa")
sapply(packages, install_if_missing)

# OPTIONAL: install "msa" if not present
if (!requireNamespace("msa", quietly = TRUE)) {
  BiocManager::install("msa", type="binary", ask=FALSE)
}

# list of packages needed
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

# ------------------------------------------------------------------------------
# UI Creation
ui <- dashboardPage(
  dashboardHeader(title = "Z-GENIE"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Home", tabName="home", icon=icon("home")),
      menuItem("Run and Process", tabName="run_process", icon=icon("cogs")),
      menuItem("Visualization", tabName="visualization", icon=icon("chart-bar")),
      menuItem("MSA and Tree", tabName="msa_tree", icon=icon("tree"))
    )
  ),
  
  dashboardBody(
    useShinyjs(),
    tags$head(tags$style(HTML("
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
    "))),
    
    tabItems(
      # Home
      tabItem(tabName="home",
              fluidRow(
                box(
                  title="Welcome to Z-GENIE", status="primary", solidHeader=TRUE, width=12,
                  h4("Z-GENIE (Z-DNA GENomic Information Extractor)"),
                  p("This tool helps analyze Z-DNA genomic information using custom FASTA sequences or Z-Hunt output. You can either fetch sequences from NCBI or manually input data."),
                  p("Adapted from: Ho, Pui S., et al. 'A computer aided thermodynamic approach for predicting the formation of Z‐DNA in naturally occurring sequences.' The EMBO journal 5.10 (1986): 2737-2744.")
                )
              )
      ),
      
      # Run and Process
      tabItem(tabName="run_process",
              fluidRow(
                box(
                  title="Step 1: Run and Process", status="info", solidHeader=TRUE, width=12, collapsible=TRUE,
                  
                  radioButtons("input_source", "Select Input Method:",
                               choices=list("Fetch FASTA from NCBI"="fetch",
                                            "Upload FASTA File"="upload",
                                            "Skip Fetch and Run Z-Hunt"="manual",
                                            "Manual DNA Sequence"="manual_seq"),
                               selected="fetch"
                  ),
                  
                  # FETCH
                  conditionalPanel(
                    condition = "input.input_source=='fetch'",
                    textInput("nucleotide_id","Enter Nucleotide ID (e.g., \"U81553.1\")",value=""),
                    actionButton("fetch","Fetch and Save FASTA"),
                    downloadButton("download_fasta","Download Original FASTA"),
                    textInput("fasta_path","Enter Path to FASTA File","Path/to/file/yourfile.fasta"),
                    textInput("params","Z-Hunt Parameters","8, 6, 8"),
                    actionButton("run","Run Z-Hunt"),
                    numericInput("zscore_threshold","Manual Z-Score Threshold",600)
                  ),
                  
                  # MANUAL
                  conditionalPanel(
                    condition="input.input_source=='manual'",
                    fileInput("manual_fasta_path","Upload Original FASTA File (with >)",accept=c(".fasta"),placeholder="No file selected"),
                    fileInput("manual_zscore_path","Upload Z-SCORE File",accept=c(".Z-SCORE"),placeholder="No file selected"),
                    numericInput("manual_zscore_threshold","Manual Z-Score Threshold",600)
                  ),
                  
                  # MANUAL SEQ
                  conditionalPanel(
                    condition="input.input_source=='manual_seq'",
                    textInput("organism_name","Organism Name for FASTA header","Organism"),
                    textAreaInput("dna_sequence","Enter DNA Sequence (A,T,C,G,N only)",value="",rows=5,placeholder="e.g. ATCGNNN..."),
                    actionButton("create_fasta","Create FASTA from typed sequence (for Z-Hunt)"),
                    # NEW TEXT FIELD to store or show path for Manual DNA:
                    textInput("fasta_path2","Enter Path to FASTA File","Path/to/file/yourfile.fasta"),
                    textInput("params_seq","Z-Hunt Parameters","8, 6, 8"),
                    actionButton("run","Run Z-Hunt"),
                    numericInput("zscore_threshold_seq","Manual Z-Score Threshold",600)
                  ),
                  
                  # UPLOAD
                  conditionalPanel(
                    condition="input.input_source=='upload'",
                    radioButtons("upload_mode","FASTA Input Style:",choices=c("File Upload"="file","Local File Path"="path"),selected="file"),
                    
                    conditionalPanel(
                      condition="input.upload_mode=='file'",
                      fileInput("upload_fasta","Upload FASTA File",accept=c(".fasta"))
                    ),
                    conditionalPanel(
                      condition="input.upload_mode=='path'",
                      textInput("upload_fasta_path","Enter Path to FASTA File","/path/to/large_file.fasta")
                    ),
                    textInput("params","Z-Hunt Parameters","8, 6, 8"),
                    actionButton("run","Run Z-Hunt"),
                    numericInput("zscore_threshold","Manual Z-Score Threshold",600)
                  ),
                  actionButton("process","Process Z-SCORE File"),
                  textInput("download_name","Enter Download File Name","Filtered_Z_GENIE_output"),
                  downloadButton("download","Download Processed Data"),
                  
                  # Additional naming
                  textInput("download_zscore_name","Z-SCORE Download Name","fasta.Z.SCORE"),
                  downloadButton("download_zscore","Download Z-SCORE Output"),
                  textInput("download_bed_name","BED Download Name","Z-GENIE_Minimal.bed"),
                  
                  actionButton("generate_bed","Generate BED File (Minimal)"),
                  downloadButton("download_bed","Download BED File"),
                  
                  DTOutput("data_table"),
                  verbatimTextOutput("fetch_output"),
                  verbatimTextOutput("output"),
                  verbatimTextOutput("command_output"),
                  verbatimTextOutput("final_output")
                )
              )
      ),
      
      # Visualization
      tabItem(tabName="visualization",
              fluidRow(
                box(
                  title="ZFS Visualization", status="info", solidHeader=TRUE, width=12, collapsible=TRUE,
                  sidebarLayout(
                    sidebarPanel(
                      fileInput("file","Upload CSV file",buttonLabel="Browse..."),
                      hr(),
                      h4("Filter DataTable"),
                      hr(),
                      h4("Plotly Configuration"),
                      selectInput("x","X-axis:",choices=NULL,selected="start"),
                      selectInput("y","Y-axis:",choices=NULL,selected="log10_ZScore"),
                      selectInput("color","Color by:",choices=NULL,selected="ISS"),
                      DTOutput("table",width="auto"),
                      downloadButton("download2","Download Filtered Data")
                    ),
                    mainPanel(plotlyOutput("plot"))
                  )
                )
              )
      ),
      
      # MSA and Tree
      tabItem(tabName="msa_tree",
              fluidRow(
                column(width=12,
                       box(title="Multiple Sequence Alignment",status="info",solidHeader=TRUE,width=12,
                           selectInput("alignment_method","Alignment Method",choices=c("ClustalW","ClustalOmega","Muscle"),selected="ClustalW"),
                           msaROutput("msa",width="100%",height="auto")
                       )
                )
              ),
              fluidRow(
                column(width=12,
                       box(title="Phylogenetic Tree",status="info",solidHeader=TRUE,width=12,
                           plotOutput("tree",height="1500px"))
                )
              )
      )
    )
  )
)

# -------------------------------------------------------------------
# SERVER

server <- function(input, output, session){
  
  processed_data_reactive <- reactiveVal(NULL)
  bed_data_reactive <- reactiveVal(NULL)
  
  # A reactive that stores the path to .Z-SCORE file
  zscore_path_reactive <- reactiveVal(NULL)
  
  # A reactive that stores the path to the manually created FASTA
  manual_fasta_path <- reactiveVal(NULL)
  
  # Helper function: Preprocess FASTA (remove >)
  preprocess_fasta <- function(file_path){
    fasta_lines <- readLines(file_path)
    if (startsWith(fasta_lines[1], ">")){
      fasta_lines <- fasta_lines[-1]
      writeLines(fasta_lines,file_path)
    }
  }
  
  # Fetch FASTA from NCBI
  fetch_and_save_fasta <- function(nucleotide_id){
    fasta <- rentrez::entrez_fetch(db="nucleotide",
                                   id=eval(expression(text=nucleotide_id)),
                                   rettype="fasta")
    file_name <- paste0(nucleotide_id,".fasta")
    write(fasta,file=file_name)
    return(file_name)
  }
  
  # Remove first line if starts with "/"
  remove_first_line <- function(file_path){
    lines <- readLines(file_path)
    if (startsWith(lines[1],"/")){
      lines <- lines[-1]
    }
    writeLines(lines,file_path)
  }
  
  # Run Z-Hunt
  run_zhunt <- function(input_file, params){
    Z_path <- normalizePath(file.path(getwd(),"zhunt/bin/zhunt"), mustWork=FALSE)
    
    if (!file.exists(Z_path)){
      stop("Error: The zHunt binary was not found at path:", Z_path)
    }
    system(paste("chmod +x",Z_path),ignore.stderr=TRUE)
    if (file.access(Z_path, mode=1)!=0){
      stop("Error: The zHunt binary is not executable. Check file permissions. Path:",Z_path)
    }
    
    modified_fasta_file <- paste0(tools::file_path_sans_ext(input_file),"_mod.fasta")
    file.copy(input_file,modified_fasta_file,overwrite=TRUE)
    preprocess_fasta(modified_fasta_file)
    
    command <- paste(Z_path, paste(params,collapse=" "),modified_fasta_file)
    cat("Executing command:", command, "\n")
    
    result <- tryCatch({
      processx::run(
        command=Z_path,
        args=c(params,modified_fasta_file),
        error_on_status=FALSE,
        echo_cmd=TRUE,
        echo=TRUE
      )
    }, error=function(e){
      return(paste("Error running Z-Hunt:",e$message))
    })
    
    if (!is.null(result$stderr) && nzchar(result$stderr)){
      return(paste("Error occurred during Z-Hunt execution:\n",result$stderr))
    }
    return(result$stdout)
  }
  
  # ObserveEvent: create FASTA from manual DNA
  observeEvent(input$create_fasta,{
    withProgress(message="Creating FASTA from typed sequence...", value=0, {
      time_start <- Sys.time()
      incProgress(0.3, detail="Validating sequence")
      
      req(input$organism_name)
      req(input$dna_sequence)
      seq_raw <- gsub("\\s+","",input$dna_sequence)
      if (nchar(seq_raw)<1){
        showNotification("Please enter a DNA sequence.", type="error")
        return()
      }
      if (grepl("[^atcgnATCGN]", seq_raw)){
        showNotification("Invalid characters found. Only A,T,C,G,N allowed.", type="error")
        return()
      }
      
      incProgress(0.6, detail="Writing to temporary file")
      manual_file <- tempfile(pattern="ManualDNA_",fileext=".fasta")
      cat(">",input$organism_name,"\n",seq_raw,"\n",file=manual_file,sep="")
      manual_fasta_path(manual_file)
      
      incProgress(0.9, detail="Updating path field")
      # We update the separate 'fasta_path2' for manual_seq
      updateTextInput(session,"fasta_path2",value=normalizePath(manual_file))
      
      time_end <- Sys.time()
      total_time <- round(as.numeric(difftime(time_end,time_start,units="secs")),2)
      showNotification(paste("FASTA created in", total_time,"seconds at:",manual_file), type="message")
    })
  })
  
  # ObserveEvent: fetch
  observeEvent(input$fetch,{
    withProgress(message="Fetching FASTA from NCBI...", value=0, {
      time_start <- Sys.time()
      incProgress(0.2, detail="Checking Nucleotide ID")
      req(input$nucleotide_id)
      
      incProgress(0.4, detail="Fetching from NCBI")
      original_fasta_file <- fetch_and_save_fasta(input$nucleotide_id)
      updateTextInput(session,"fasta_path",value=normalizePath(original_fasta_file))
      incProgress(0.8, detail="Finishing up")
      
      showNotification("FASTA file fetched successfully!",type="message")
      
      output$download_fasta <- downloadHandler(
        filename=function(){ "fasta.fasta" },
        content=function(file){ file.copy(original_fasta_file, file) },
        contentType="text/plain"
      )
      time_end <- Sys.time()
      total_time <- round(as.numeric(difftime(time_end,time_start,units="secs")),2)
      showNotification(paste("Fetch completed in", total_time,"seconds."), type="message")
    })
  })
  
  # ObserveEvent: run Z-Hunt
  observeEvent(input$run,{
    withProgress(message="Running Z-Hunt...", value=0, {
      time_start <- Sys.time()
      incProgress(0.2, detail="Determining input source")
      
      # We handle logic for each input_source
      if (input$input_source=="upload"){
        if (input$upload_mode=="file"){
          req(input$upload_fasta)
          original_fasta_file <- input$upload_fasta$datapath
        } else {
          req(input$upload_fasta_path)
          original_fasta_file <- input$upload_fasta_path
        }
        params <- as.numeric(unlist(strsplit(input$params,",")))
        
      } else if (input$input_source=="fetch"){
        req(input$fasta_path)
        original_fasta_file <- normalizePath(input$fasta_path)
        params <- as.numeric(unlist(strsplit(input$params,",")))
        
      } else if (input$input_source=="manual_seq"){
        # If manual_seq, we can check manual_fasta_path(); if not null, that means user used create_fasta
        # Otherwise, fallback to input$fasta_path2
        if (!is.null(manual_fasta_path())){
          original_fasta_file <- manual_fasta_path()
        } else {
          original_fasta_file <- input$fasta_path2
        }
        params <- as.numeric(unlist(strsplit(input$params_seq,",")))
        
      } else {
        # for 'manual' we skip running Z-hunt
        output$output <- renderText("Z-Hunt run not applicable for 'manual' input method.")
        return()
      }
      
      incProgress(0.4, detail="Validating parameters")
      if (length(params)!=3){
        output$output <- renderText("Error: 'params' must have exactly three numbers.")
        return()
      }
      
      incProgress(0.6, detail="Invoking run_zhunt()")
      command_output <- tryCatch({
        run_zhunt(original_fasta_file, params)
      }, error=function(e){
        paste("Error running Z-Hunt:", e$message)
      })
      
      incProgress(0.8, detail="Finalizing results")
      output$output <- renderText({
        paste("Z-Hunt run complete for file:", original_fasta_file,
              "\nParameters used:", paste(params,collapse=" "),
              "\n\nOutput:\n", command_output)
      })
      
      output_file_path <- paste0(tools::file_path_sans_ext(original_fasta_file), "_mod.fasta.Z-SCORE")
      zscore_path_reactive(output_file_path)
      
      time_end <- Sys.time()
      total_time <- round(as.numeric(difftime(time_end,time_start,units="secs")),2)
      showNotification(paste("Z-Hunt completed in", total_time,"seconds."), type="message")
    })
  })
  
  # Download Z-SCORE
  output$download_zscore <- downloadHandler(
    filename=function(){
      validate(need(nzchar(input$download_zscore_name),"Please enter a file name for Z-SCORE."))
      input$download_zscore_name
    },
    content=function(file){
      req(zscore_path_reactive())
      real_path <- zscore_path_reactive()
      if (!file.exists(real_path)){
        stop("Z-SCORE file not found. Please run Z-Hunt first.")
      }
      file.copy(real_path, file)
    },
    contentType="text/plain"
  )
  
  # ---------------------------
  # Full "Process Z-SCORE File"
  # ---------------------------
  observeEvent(input$process,{
    withProgress(message="Processing Z-SCORE File (FULL)", value=0, {
      time_start <- Sys.time()
      incProgress(0.2, detail="Identifying input source")
      
      req(input$input_source)
      zscore_file <- NULL
      
      if (input$input_source=="manual"){
        req(input$manual_fasta_path, input$manual_zscore_path)
        original_fasta_file <- input$manual_fasta_path$datapath
        zscore_file <- input$manual_zscore_path$datapath
        threshold <- input$manual_zscore_threshold
        
        if (!file.exists(zscore_file)){
          showNotification("The .fasta.Z-SCORE file was not found.", type="error")
          return()
        }
        remove_first_line(zscore_file)
        
      } else if (input$input_source=="upload"){
        if (input$upload_mode=="file"){
          req(input$upload_fasta)
          original_fasta_file <- input$upload_fasta$datapath
        } else {
          req(input$upload_fasta_path)
          original_fasta_file <- input$upload_fasta_path
        }
        zscore_file <- paste0(tools::file_path_sans_ext(original_fasta_file), "_mod.fasta.Z-SCORE")
        threshold <- input$zscore_threshold
        
        if (!file.exists(zscore_file)){
          showNotification("The .fasta.Z-SCORE file was not found.", type="error")
          return()
        }
        remove_first_line(zscore_file)
        
      } else if (input$input_source=="fetch"){
        req(input$fasta_path)
        original_fasta_file <- normalizePath(input$fasta_path)
        zscore_file <- paste0(tools::file_path_sans_ext(original_fasta_file), "_mod.fasta.Z-SCORE")
        threshold <- input$zscore_threshold
        
        if (!file.exists(zscore_file)){
          showNotification("The .fasta.Z-SCORE file was not found.", type="error")
          return()
        }
        remove_first_line(zscore_file)
        
      } else if (input$input_source=="manual_seq"){
        # If manual_seq, again check if manual_fasta_path() is not null, else fallback
        if (!is.null(manual_fasta_path())){
          original_fasta_file <- manual_fasta_path()
        } else {
          original_fasta_file <- input$fasta_path2
        }
        
        zscore_file <- paste0(tools::file_path_sans_ext(original_fasta_file), "_mod.fasta.Z-SCORE")
        threshold <- input$zscore_threshold_seq
        
        if (!file.exists(zscore_file)){
          showNotification("The .fasta.Z-SCORE file was not found for manual DNA sequence. Please run Z-Hunt first.", type="error")
          return()
        }
        remove_first_line(zscore_file)
      }
      
      incProgress(0.4, detail="Loading Z-SCORE data")
      processed_data <- fread(zscore_file, col.names=c("V1","V2","zscore","conformation"))
      processed_data$Position <- 1:nrow(processed_data)
      processed_data$max_Z_Score <- processed_data$zscore
      processed_data$OligoConformation <- processed_data$conformation
      
      threshold <- as.numeric(threshold)
      processed_data <- processed_data[order(processed_data$max_Z_Score, decreasing=TRUE),]
      end <- tail(which(processed_data$max_Z_Score>=threshold),1)
      if (length(end)==0){
        showNotification("No scores meet the threshold criteria.", type="error")
        return()
      }
      processed_data <- processed_data[1:end,]
      processed_data <- processed_data[order(processed_data$Position),]
      
      # param
      if (input$input_source=="manual_seq"){
        param_vec <- as.numeric(unlist(strsplit(input$params_seq,",")))
      } else {
        param_vec <- as.numeric(unlist(strsplit(input$params,",")))
      }
      if (length(param_vec)==3){
        middle_value <- param_vec[2]
        limit_value <- 2*middle_value
      } else {
        stop("Error: 'params' should have exactly three numbers.")
      }
      
      incProgress(0.6, detail="Grouping consecutive positions in R")
      positions <- processed_data$Position
      if (length(positions)<1){
        showNotification("No positions after threshold filter.", type="warning")
        return()
      }
      
      groups <- list()
      temp <- c(positions[1])
      for (i in seq(2,length(positions))){
        if (abs(positions[i] - positions[i-1]) < limit_value){
          temp <- c(temp, positions[i])
        } else {
          groups[[length(groups)+1]] <- temp
          temp <- c(positions[i])
        }
      }
      groups[[length(groups)+1]] <- temp
      
      incProgress(0.7, detail="Building final data frame")
      processed_fasta <- read.fasta(original_fasta_file)
      dna <- toString(processed_fasta[[names(processed_fasta)]])
      dna <- gsub(", ","", dna)
      dna <- toupper(dna)
      dna <- as.character(dna)
      c1_fasta <- s2c(dna)
      
      MasterDF <- data.frame(
        "start"=numeric(0),
        "end"=numeric(0),
        "range"=character(0),
        "sequence"=character(0),
        "oligo_length"=numeric(0),
        "GC_of_sequence"=numeric(0),
        "max_Z_Score"=numeric(0)
      )
      
      idx <- 1
      for (grp in groups){
        start_pos <- head(grp,1)
        last_pos  <- tail(grp,1)
        last_row_idx <- which(processed_data$Position==last_pos)
        if (length(last_row_idx)<1){
          showNotification(paste("No row found for last_pos=",last_pos),type="warning")
          next
        }
        length_oligo <- nchar(processed_data$OligoConformation[last_row_idx])
        end_pos <- last_pos + length_oligo - 1
        
        rng <- paste0(start_pos,":",end_pos)
        seq_sub <- toString(c1_fasta[ eval(parse(text=rng)) ])
        seq_sub <- gsub(", ","", seq_sub)
        seq_sub <- toupper(seq_sub)
        
        seq_len <- nchar(seq_sub)
        seq_s2c <- s2c(seq_sub)
        GCpct <- GC(seq_s2c)*100
        
        row_indices <- which(processed_data$Position %in% grp)
        local_max_z <- max(processed_data$max_Z_Score[row_indices])
        
        MasterDF[idx, ] <- c(
          start_pos,
          end_pos,
          rng,
          seq_sub,
          seq_len,
          round(GCpct,2),
          local_max_z
        )
        idx <- idx+1
      }
      
      MasterDF$start          <- as.numeric(MasterDF$start)
      MasterDF$end            <- as.numeric(MasterDF$end)
      MasterDF$oligo_length   <- as.numeric(MasterDF$oligo_length)
      MasterDF$GC_of_sequence <- as.numeric(MasterDF$GC_of_sequence)
      MasterDF$max_Z_Score    <- as.numeric(MasterDF$max_Z_Score)
      MasterDF$log10_ZScore   <- log(MasterDF$max_Z_Score, base=10)
      
      # Additional calculations
      calculate_at_content <- function(sequence){
        at_count <- sum(str_count(sequence,"AT") + str_count(sequence,"TA"))
        total_count <- nchar(sequence) - 1
        if (total_count<1) return(0)
        at_content <- at_count/total_count
        return(at_content)
      }
      calculate_gc_content <- function(sequence){
        gc_count <- sum(str_count(sequence,"GC") + str_count(sequence,"CG"))
        total_count <- nchar(sequence) - 1
        if (total_count<1) return(0)
        gc_content <- gc_count/total_count
        return(gc_content)
      }
      calculate_atgc_content <- function(sequence){
        at_count <- sum(str_count(sequence,"A") + str_count(sequence,"T"))
        gc_count <- sum(str_count(sequence,"G") + str_count(sequence,"C"))
        if (gc_count==0) return(Inf)
        atgc_content <- at_count/gc_count
        return(atgc_content)
      }
      
      MasterDF$sequence <- as.character(MasterDF$sequence)
      MasterDF <- MasterDF %>%
        mutate(
          at_content  = sapply(sequence, calculate_at_content)*100,
          gc_content  = sapply(sequence, calculate_gc_content)*100,
          atgc_content= sapply(sequence, calculate_atgc_content)
        )
      
      find_purine_pattern <- function(sequence){
        pattern <- "(A|G)(A|G)(CG)+?(C|T)(C|T)"
        matches <- gregexpr(pattern,sequence,perl=TRUE)
        if (length(matches[[1]])==1 && matches[[1]][1]==-1){
          return(data.frame(start=integer(0),end=integer(0),sequence=character(0)))
        }
        matching_sequences <- regmatches(sequence,matches)
        start_positions <- matches[[1]]
        end_positions <- start_positions + attr(matches[[1]],"match.length") -1
        data.frame(start=start_positions,end=end_positions,sequence=matching_sequences[[1]])
      }
      
      analyze_sequences <- function(df){
        df$ISS <- NA
        df$ISSeq <- NA
        for (i in seq_len(nrow(df))){
          sequence <- df$sequence[i]
          match_result <- find_purine_pattern(sequence)
          if (nrow(match_result)>0){
            df$ISS[i] <- "TRUE"
            df$ISSeq[i] <- paste(match_result$sequence,collapse=",")
          } else {
            df$ISS[i] <- "FALSE"
          }
        }
        df
      }
      
      find_motif_pattern <- function(sequence){
        pattern <- "(A|G)(A|G)(A|G)(C)(A|T)(A|T)(G)(C|T)(C|T)(C|T)"
        matches <- gregexpr(pattern,sequence,perl=TRUE)
        if (length(matches[[1]])==1 && matches[[1]][1]==-1){
          return(data.frame(start=integer(0),end=integer(0),sequence=character(0)))
        }
        matching_sequences <- regmatches(sequence,matches)
        start_positions <- matches[[1]]
        end_positions <- start_positions + attr(matches[[1]],"match.length") -1
        data.frame(start=start_positions,end=end_positions,sequence=matching_sequences[[1]])
      }
      
      analyze_sequences2 <- function(df){
        df$p53motif <- NA
        df$p53motifeq <- NA
        for (i in seq_len(nrow(df))){
          sequence <- df$sequence[i]
          match_result <- find_motif_pattern(sequence)
          if (nrow(match_result)>0){
            df$p53motif[i] <- "TRUE"
            df$p53motifeq[i] <- paste(match_result$sequence,collapse=",")
          } else {
            df$p53motif[i] <- "FALSE"
          }
        }
        df
      }
      
      find_palindromes <- function(sequence){
        dna <- DNAString(sequence)
        palindromes <- findPalindromes(dna)
        if (length(palindromes)==0){
          return(NULL)
        } else {
          palindrome_sequences <- as.character(Views(dna,palindromes))
          return(paste(palindrome_sequences,collapse=", "))
        }
      }
      
      incProgress(0.8, detail="Applying additional calculations")
      result_df <- analyze_sequences(MasterDF)
      result_df <- analyze_sequences2(result_df)
      result_df$palindromic_sequence <- sapply(result_df$sequence, find_palindromes)
      result_df$palindromic_sequence <- as.character(result_df$palindromic_sequence)
      result_df <- result_df %>% mutate_if(is.numeric, ~ round(.,2))
      
      dna_length <- nchar(dna)
      Zp <- (sum(result_df$oligo_length)/dna_length)*100
      output$final_output <- renderText({
        paste("Z-Potentiality is", round(Zp,2), "%")
      })
      
      output$data_table <- renderDT({
        datatable(result_df, options=list(pageLength=10, searchHighlight=TRUE))
      })
      
      processed_data_reactive(result_df)
      
      output$download <- downloadHandler(
        filename=function(){ paste0(input$download_name,".csv") },
        content=function(fname){ write.csv(result_df,fname,row.names=FALSE) }
      )
      
      incProgress(1.0, detail="Done")
      time_end <- Sys.time()
      total_time <- round(as.numeric(difftime(time_end,time_start,units="secs")),2)
      showNotification(paste("Processing (FULL) complete in", total_time,"seconds. Table updated."),type="message")
    })
  })
  
  # ---------------------------
  # Generate BED File (Minimal)
  # ---------------------------
  observeEvent(input$generate_bed,{
    withProgress(message="Generating BED File (Minimal)", value=0, {
      time_start <- Sys.time()
      incProgress(0.2, detail="Identifying input source")
      
      if (input$input_source=="manual"){
        req(input$manual_zscore_path)
        zscore_file <- input$manual_zscore_path$datapath
        threshold <- input$manual_zscore_threshold
        if (!file.exists(zscore_file)){
          showNotification("The .fasta.Z-SCORE file was not found. Cannot generate BED.", type="error")
          return()
        }
        remove_first_line(zscore_file)
        
      } else if (input$input_source=="upload"){
        if (input$upload_mode=="file"){
          req(input$upload_fasta)
          original_fasta_file <- input$upload_fasta$datapath
        } else {
          req(input$upload_fasta_path)
          original_fasta_file <- input$upload_fasta_path
        }
        zscore_file <- paste0(tools::file_path_sans_ext(original_fasta_file), "_mod.fasta.Z-SCORE")
        threshold <- input$zscore_threshold
        if (!file.exists(zscore_file)){
          showNotification("The .fasta.Z-SCORE file was not found. Cannot generate BED.", type="error")
          return()
        }
        remove_first_line(zscore_file)
        
      } else if (input$input_source=="fetch"){
        req(input$fasta_path)
        original_fasta_file <- normalizePath(input$fasta_path)
        zscore_file <- paste0(tools::file_path_sans_ext(original_fasta_file), "_mod.fasta.Z-SCORE")
        threshold <- input$zscore_threshold
        if (!file.exists(zscore_file)){
          showNotification("The .fasta.Z-SCORE file was not found. Cannot generate BED.", type="error")
          return()
        }
        remove_first_line(zscore_file)
        
      } else if (input$input_source=="manual_seq"){
        if (!is.null(manual_fasta_path())){
          original_fasta_file <- manual_fasta_path()
        } else {
          original_fasta_file <- input$fasta_path2
        }
        zscore_file <- paste0(tools::file_path_sans_ext(original_fasta_file), "_mod.fasta.Z-SCORE")
        threshold <- input$zscore_threshold_seq
        if (!file.exists(zscore_file)){
          showNotification("The .fasta.Z-SCORE file was not found for manual DNA sequence. Please run Z-Hunt first.",type="error")
          return()
        }
        remove_first_line(zscore_file)
        
      } else {
        showNotification("No valid input source for Z-SCORE found.",type="error")
        return()
      }
      
      incProgress(0.4, detail="Reading Z-SCORE data")
      bed_data <- fread(zscore_file, col.names=c("V1","V2","zscore","conformation"))
      bed_data$Position <- 1:nrow(bed_data)
      bed_data$max_Z_Score <- bed_data$zscore
      bed_data$OligoConformation <- bed_data$conformation
      
      threshold <- as.numeric(threshold)
      bed_data <- bed_data[order(bed_data$max_Z_Score, decreasing=TRUE),]
      end_idx <- tail(which(bed_data$max_Z_Score>=threshold),1)
      if (length(end_idx)==0){
        showNotification("No scores meet threshold for BED generation.", type="error")
        return()
      }
      bed_data <- bed_data[1:end_idx,]
      bed_data <- bed_data[order(bed_data$Position),]
      
      if (input$input_source=="manual_seq"){
        param_vec <- as.numeric(unlist(strsplit(input$params_seq,",")))
      } else {
        param_vec <- as.numeric(unlist(strsplit(input$params,",")))
      }
      if (length(param_vec)==3){
        middle_value <- param_vec[2]
        limit_value <- 2*middle_value
      } else {
        showNotification("Params must have exactly three values. Cannot generate BED.", type="error")
        return()
      }
      
      incProgress(0.6, detail="Performing pure R grouping for BED")
      positions <- bed_data$Position
      if (length(positions)<1){
        showNotification("No positions found after threshold filter for BED.", type="warning")
        return()
      }
      groups <- list()
      temp <- c(positions[1])
      for (i in seq(2,length(positions))){
        if (abs(positions[i]-positions[i-1])<limit_value){
          temp <- c(temp, positions[i])
        } else {
          groups[[length(groups)+1]] <- temp
          temp <- c(positions[i])
        }
      }
      groups[[length(groups)+1]] <- temp
      
      incProgress(0.8, detail="Building minimal BED DataFrame")
      minimal_bed_df <- data.frame(start=numeric(0), end=numeric(0), max_Z_Score=numeric(0))
      idx <- 1
      for (grp in groups){
        start_pos <- head(grp,1)
        last_pos  <- tail(grp,1)
        
        row_idx <- which(bed_data$Position==last_pos)
        if (length(row_idx)<1){
          showNotification(paste("No row found for last_pos=",last_pos,". Skipping group."), type="warning")
          next
        }
        length_oligo <- nchar(bed_data$OligoConformation[row_idx])
        end_pos <- last_pos + length_oligo -1
        
        grp_rows <- which(bed_data$Position %in% grp)
        local_max_zscore <- max(bed_data$max_Z_Score[grp_rows])
        
        minimal_bed_df[idx, ] <- c(start_pos, end_pos, local_max_zscore)
        idx <- idx+1
      }
      
      chrom_name <- if (!is.null(input$nucleotide_id) && nzchar(input$nucleotide_id)){
        input$nucleotide_id
      } else {
        "chrom"
      }
      
      bed_final <- data.frame(
        chrom=chrom_name,
        start=minimal_bed_df$start,
        end=minimal_bed_df$end,
        score=minimal_bed_df$max_Z_Score
      )
      
      bed_data_reactive(bed_final)
      
      output$download_bed <- downloadHandler(
        filename=function(){
          validate(need(nzchar(input$download_bed_name),"Please enter a file name for BED."))
          input$download_bed_name
        },
        content=function(file){
          write.table(bed_final, file, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
        }
      )
      
      incProgress(1.0, detail="Done")
      time_end <- Sys.time()
      total_time <- round(as.numeric(difftime(time_end,time_start,units="secs")),2)
      showNotification(paste("Minimal BED file generated in", total_time,"seconds. Click 'Download BED File'."), type="message")
    })
  })
  
  # Visualization
  data <- reactiveVal(NULL)
  
  observeEvent(input$file,{
    req(input$file)
    df <- read.csv(input$file$datapath, header=TRUE) %>%
      mutate_all(~ if(is.character(.)) as.factor(.) else as.numeric(.))
  # mutate_all(~ ifelse(is.character(.), as.factor(.),  as.numeric(.)))
                 
    updateSelectInput(session,"x",choices=names(df),selected="start")
    updateSelectInput(session,"y",choices=names(df),selected="log10_ZScore")
    updateSelectInput(session,"color",choices=names(df),selected="ISS")
    data(df)
  })
  
  output$table <- renderDT({
    req(data())
    datatable(data(), filter='top', selection='none')
  })
  
  output$plot <- renderPlotly({
    req(data())
    x_var <- input$x
    y_var <- input$y
    color_var <- input$color
    
    plot_ly(
      data=data(),
      x=~ get(x_var),
      y=~ get(y_var),
      color=~ get(color_var),
      type='scatter',
      mode='markers',
      text=~ paste(
        "Sequence: ", data()$sequence,
        "<br>Start: ", data()$start,
        "<br>End: ", data()$end,
        "<br>Oligo Length: ", data()$oligo_length,
        "<br>GC Content: ", data()$GC_of_sequence,
        "<br>AT Content: ", data()$at_content,
        "<br>Max Z-Score: ", data()$max_Z_Score,
        "<br>ISS: ", data()$ISS,
        "<br>ISSeq: ", data()$ISSeq,
        "<br>Palindromic Sequence: ", data()$palindromic_sequence
      )
    ) %>%
      layout(
        dragmode='select',
        xaxis=list(title=x_var),
        yaxis=list(title=y_var)
      )
  })
  
  observe({
    req(data())
    if (!is.null(input$table_rows_all)){
      filtered_data <- data()[ input$table_rows_all, ]
      output$plot <- renderPlotly({
        req(filtered_data)
        x_var <- input$x
        y_var <- input$y
        color_var <- input$color
        
        plot_ly(
          data=filtered_data,
          x=~ get(x_var),
          y=~ get(y_var),
          color=~ get(color_var),
          type='scatter',
          mode='markers',
          text=~ paste(
            "Sequence: ", filtered_data$sequence,
            "<br>Start: ", filtered_data$start,
            "<br>End: ", filtered_data$end,
            "<br>Oligo Length: ", filtered_data$oligo_length,
            "<br>GC Content: ", filtered_data$GC_of_sequence,
            "<br>AT Content: ", filtered_data$at_content,
            "<br>Max Z-Score: ", filtered_data$max_Z_Score,
            "<br>ISS: ", filtered_data$ISS,
            "<br>ISSeq: ", filtered_data$ISSeq,
            "<br>Palindromic Sequence: ", filtered_data$palindromic_sequence
          )
        ) %>%
          layout(
            dragmode='select',
            xaxis=list(title=x_var),
            yaxis=list(title=y_var)
          )
      })
    }
  })
  
  observe({
    req(data())
    if (!is.null(input$table_rows_all)){
      filtered_data <- data()[ input$table_rows_all, ]
      req(filtered_data)
      
      output$download2 <- downloadHandler(
        filename=function(){ paste0(input$download_name,"_Filtered.csv") },
        content=function(fname){ write.csv(filtered_data,fname,row.names=FALSE) }
      )
      
      output$msa <- renderMsaR({
        req(filtered_data)
        cdr3aa_reactive <- reactive({
          req(filtered_data)
          t <- as.data.frame(filtered_data[,4])
          colnames(t)<-NULL
          t[,1]<-as.data.frame(toupper(t[,1]))
          alignment_method <- input$alignment_method
          aa_msa <- msa(t[,1],method=alignment_method,type="protein",order="input")
          AA_msa <- msaConvert(aa_msa, type="seqinr::alignment")
          aa_msa <- AA_msa
          aa_msa$nam <- t[,1]
          AA_msa <- as.data.frame(AA_msa[["seq"]])
          colnames(AA_msa)<-NULL
          AA_msa <- as.data.frame(t(AA_msa))
          colnames(AA_msa)<-AA_msa[1,]
          write.fasta(sequences=AA_msa, names=names(AA_msa), file.out="filtered_sequences.fasta")
          
          output$tree <- renderPlot({
            req(aa_msa)
            d <- dist.alignment(aa_msa,"identity")
            as.matrix(d)
            hemoTree <- njs(d)
            plot(hemoTree, main="Phylogenetic Tree of Filtered Sequences")
          })
          
          proteins <- ape::read.FASTA("filtered_sequences.fasta", type="AA")
          return(proteins)
        })
        req(cdr3aa_reactive())
        msaR(cdr3aa_reactive(),
             colorscheme="clustal",
             conservation=TRUE,
             labelNameLength=300,
             overviewboxWidth="fixed",
             overviewboxHeight="fixed",
             alignmentHeight=500,
             rowheight=30,
             height=500,
             width=1500)
      })
    }
  })
}

shinyApp(ui, server)



