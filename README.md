# Z-GENIE
Z-GENIE (Z-DNA GENomic Information Extractor)

####Installation Instructions

#Clone the repository

git clone https://github.com/agarzare/Z-GENIE.git

#Install the package locally

install.packages("devtools")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install("Biostrings")

BiocManager::install("msa", type = "binary")

devtools::install("/path/to/cloned/ZGENIE")

#Run Z-GENIE

library(ZGENIE)

library(shiny)

runApp(system.file("shinyapp", package = "ZGENIE"))
