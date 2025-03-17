# Z-GENIE (Z-DNA GENomic Information Extractor)

Z-GENIE is a Shiny application designed for interactive bioinformatics analysis, including sequence alignment and visualization.

## Installation Instructions
### Install R and RStudio
If not running in browser on the shiny server, you will need to have R and RStudio installed on your device. 
Below are links on how to install each component.
Installing R: https://cran.r-project.org/
Installing RStudio: https://www.rstudio.com/



### Clone the Repository

Ways to clone:

1. HTTPS 
Input the following into terminal or command line:
```bash
git clone https://github.com/agarzare/Z-GENIE.git
```
You will be prompted to log in to your personal account.

2. SSH 
*You will need to have a public SSH key connected to your account
Input the following into terminal or command line:
```bash
git clone git@github.com:agarzare/Z-GENIE.git
```

Continue in R/RStudio from the terminal:
```r
# Install devtools
install.packages("devtools")

# Install the ZGENIE package from the cloned repository
devtools::install("/path/to/cloned/ZGENIE")

# Load the ZGENIE library
library(ZGENIE)

# Load the Shiny package
library(shiny)

# Run the Shiny app
runApp(system.file("shinyapp", package = "ZGENIE"))
```
