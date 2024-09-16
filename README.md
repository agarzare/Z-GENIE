# Z-GENIE (Z-DNA GENomic Information Extractor)

Z-GENIE is a Shiny application designed for interactive bioinformatics analysis, including sequence alignment and visualization.

## Installation Instructions

### Clone the Repository

To get started, clone the repository from GitHub:

```bash
git clone https://github.com/agarzare/Z-GENIE.git
```

Continue in R/RStudio:
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
