# Z-GENIE
Z-GENIE (Z-DNA GENomic Information Extractor)

####Installation Instructions

#Clone the repository

git clone https://github.com/agarzare/Z-GENIE.git

#Install the package locally

devtools::install("/path/to/cloned/ZGENIE")

#Run Z-GENIE

library(ZGENIE)

runApp(system.file("shinyapp", package = "ZGENIE"))
