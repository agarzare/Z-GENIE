library(shinytest2)
library(profvis)
library(callr)  # used to run the app in the background, otherwise it stalls

# first start the local app
# shiny::runApp("inst/shinyapp", host = "127.0.0.1", port = 4186, launch.browser = TRUE)
# callr::r_bg()
check_if_done <- function(nuc_id, my_csv_path){
  if (!file.exists(my_csv_path)) {
    stop("The file does not exist at the specified path.")
  }
  # Define the valid steps
  valid_steps <- c("fetch", "zhunt", "zscoreproc", "bedfile")
  
  # Load the CSV into a data frame
  data <- read.csv(my_csv_path, header = TRUE)
  
  # Filter the data for the given nuc_id
  nuc_data <- subset(data, nuc_ID == nuc_id)
  
  # Check if all steps are present for the given nuc_id
  missing_steps <- setdiff(valid_steps, nuc_data$step)
  
  if (length(missing_steps) == 0) {
    cat("The nuc_id has all required steps.\n")
    return(TRUE)
  } else {
    cat("The nuc_id is missing the following steps:", paste(missing_steps, collapse = ", "), "\n")
    return(FALSE)
  }
  
}


get_max_z <- function(nuc_id){
  # step 1: get the nuc_id+ "_mod.fasta.Z-SCORE into a dataframe
  nuc_z_file <- paste0("inst/shinyapp/", nuc_id, "_mod.fasta.Z-SCORE")
  data <- read.table(nuc_z_file, header = FALSE)
  
  # Assign column names for better readability
  colnames(data) <- c("Column1", "Column2", "Z_Score", "Sequence")
  # Step 2: get the largest value in the 3rd column
  max_z <- max(data$Z_Score)
  rounded_down <- floor(max_z / 100) * 100
  # step 3: put that value into the zscore threshold from 600 to that
  # basically wait again and proceed to next step as normal 
  return(rounded_down)
}
  
  

step_done <- function(nuc_id, step, my_csv_path){
  if (!file.exists(my_csv_path)) {
    stop("The file does not exist at the specified path.")
  }
  cat("FILE EXISTS ")
  
  # Load the CSV into a data frame
  data <- read.csv(my_csv_path, header = TRUE)
  
  # The nuc_id and step to check
  nuc_id_to_check <- nuc_id
  step_to_check <- step
  
  # Check if the pair exists in the data frame
  exists <- any(data$nuc_ID == nuc_id_to_check & data$step == step_to_check)
  # Print the result
  if (exists) {
    cat("The nuc_id and step pair exists in the CSV.\n")
    return(TRUE)
  } else {
    cat("The nuc_id and step pair does not exist in the CSV.\n")
    return(FALSE)
  }
}


test_this <- function(nuc_id, my_csv_path, run_type){
  # Create a dynamic profiling output path for each nuc_id
  prof_output_path <- paste0("~/zgenie/ZGENIE/inst/shinyapp/benchmarks/bench_", nuc_id, run_type, ".Rprof")
  # OG STUFF
  shiny_process <- r_bg(function() {
    shiny::runApp("inst/shinyapp", host = "127.0.0.1", port = 4186)
  })
  
  # Give the app some time to start
  Sys.sleep(5)
  
  # Initialize test driver for the app
  app <- AppDriver$new(
    "http://127.0.0.1:4186",  # Change if your app runs on a different port
    load_timeout = 10000  # Allow up to 10s for the app to load
  )
  cat("App started successfully\n")
  
  
  
  # go to teh right tab
  # app$set_inputs(tabName = "run_process")
  app$click(selector = "a[data-value='run_process']")
  Sys.sleep(2)
  cat("clicked the right tab \n")
  
  # setting the ID, no need to wait
  app$set_inputs(nucleotide_id = nuc_id, wait_ = FALSE)
  
  # TIMINGTEST 1: Timing the fetch, save to dataframe
  app$click("fetch")
  app$wait_for_js("Array.from(document.querySelectorAll('.shiny-notification'))
      .some(el => el.innerText.includes('FASTA file fetched successfully!'))",
                  timeout = 5000000)
  Sys.sleep(2)
  cat("Timing test 1 done\n")
  
  
  #TIMINGTEST 2: timing z hunt run
  app$click("run", wait_ = FALSE)
  app$wait_for_js("Array.from(document.querySelectorAll('.shiny-notification'))
      .some(el => el.innerText.includes('Z-Hunt completed in'))",
                  timeout = 5000000)
  Sys.sleep(2)
  cat("Timing test 2 done\n")
  
  
  
  #TIMINGTEST 3: timing the time taken for processing z score
  app$click("process", wait_ = FALSE)
  # app$wait_for_js("Array.from(document.querySelectorAll('.shiny-notification'))
  #     .some(el => el.innerText.includes('Z-Hunt completed in'))",
  #                 timeout = 5000000)
  
  app$wait_for_js(
    "document.querySelectorAll('.shiny-notification').length > 0",
    timeout = 5000000
  )
  cat("Proceed to check if step is in CSV \n")
  Sys.sleep(2)
  
  z_hunt_failed <- !step_done(nuc_id, "zscoreproc", my_csv_path)
  
  if(z_hunt_failed){
    # getting the max value 
    max_z <- get_max_z(nuc_id)
    cat("MAX Z IS BELOW ________________________")
    cat(max_z)
    if(max_z <= 0){
      cat("No available threshold! \n")
      part3 <- data.frame(nuc_id = nuc_id, step = "zscoreproc", time = NA, run_type = "", z_thresh = -1, stringsAsFactors = FALSE)
      write.table(part3, my_csv_path, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
      part4 <- data.frame(nuc_id = nuc_id, step = "bedfile", time = NA, run_type = "", z_thresh = -1, stringsAsFactors = FALSE)
      app$stop()
      return("done: no threshold")
      }
    else{
      cat("lowering threshold")
      app$set_inputs(zscore_threshold = max_z, wait_ = FALSE)
      app$click("process", wait_ = FALSE)
      # basically wait again and proceed to next step as normal
      app$wait_for_js(
        "Array.from(document.querySelectorAll('.shiny-notification'))
     .some(el => el.innerText.includes('Z-Hunt completed in'))",
        timeout = 5000000
      )
    }
    cat("Successfully lowered threshold!!!!!!!!!\n")
  }

  Sys.sleep(2)
  cat("Timing test 3 done\n")
  
  
  #TIMINGTEST 4: How long does this take BEDFILES
  app$click("generate_bed", wait_ = FALSE)
  
  app$wait_for_js("Array.from(document.querySelectorAll('.shiny-notification'))
      .some(el => el.innerText.includes('Minimal BED file generated in'))",
                  timeout = 5000000)
  Sys.sleep(2)
  cat("Timing test 4 done\n")

  # Stop the app test session after getting the last message
  app$stop()
  return("done")
}

# iterates the same ID multiple times
iterate_ids <- function(nuc_id, num_iterations, run_type){
  
  # Define the number of iterations
  num_iterations <- 4
  
  # Loop to run the test multiple times
  for (i in 1:num_iterations) {
    cat("Running test iteration:", i, "\n")
    
    #source("gettingRecLogs.R")
    tryCatch({
      test_this(nuc_id, my_csv_path, run_type)  # Replace with the actual path to getRecLogs.R
      cat("Iteration", i, "completed successfully.\n")
    }, error = function(e) {
      cat("Error in iteration", i, ":", e$message, "\n")
    })
  }
  
  cat("All test iterations completed.\n")
}

iterate_csv <- function(csv_path, my_csv_path, run_type){
  # Read the CSV file
  data <- read.csv(csv_path)
  print(names(data))
  
  # Iterate through each row and extract the "nuc_id" column
  for (i in 1:nrow(data)) {
    
    
    
    nuc_id <- data[i, "RefSeq"]
    print(nuc_id)
    if(check_if_done(nuc_id, my_csv_path)){
      # trying each one
      cat("NucID already done:", i, "\n")
    }
    else{
      tryCatch({
        test_this(nuc_id, my_csv_path, run_type)  # Replace with the actual path to getRecLogs.R
        cat("Iteration", i, "completed successfully.\n")
      }, error = function(e) {
        cat("Error in iteration", i, ":", e$message, "\n")
      })

    }

    
  }
  
}

#make sure its the same path as the app.R
#zp_path <- "zp_vals_try_3.csv"
#iterate_ids("U81553.1", 5)
csv_to_testsite <- file.path("~/zgenie/ZGENIE/inst/shinyapp", "testing_both_try_6.csv")
timing_df <- data.frame(nuc_ID = character(), step = character(), time = character(), run_type = character(), z_thresh = character(), stringsAsFactors = FALSE)
if (!file.exists(csv_to_testsite)) {
  write.csv(timing_df, csv_to_testsite, row.names = FALSE)
} 


iterate_csv("~/zgenie/ZGENIE/inst/shinyapp/sample_nucs.csv", csv_to_testsite, run_type = "both")

