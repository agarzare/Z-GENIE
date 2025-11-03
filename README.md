# Z-GENIE (Z-DNA GENomic Information Extractor)

Z-GENIE is an R/Shiny application for interactive prediction and visualization of Z-DNA–forming regions in DNA. It wraps the classic Z-Hunt algorithm in a user-friendly graphical interface so you can explore Z-DNA motifs without needing to write code.

Z-GENIE v1 is described in:

> Garza Reyna A, et al. **Z-GENIE: a user-friendly R/Shiny resource for predicting Z-DNA forming regions in DNA.** *BMC Genomics*. 2025.  
> https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-025-12148-x

---

## 1. Easiest Option: Use Z-GENIE in Your Browser

Most users can use Z-GENIE directly online with no installation.

### Launch Z-GENIE Online

<p align="center">
  <a href="https://i5h2oo-angel0i0garza0reyna.shinyapps.io/zgenie_shinyapp/" target="_blank">
    <img src="https://img.shields.io/badge/Launch%20Z--GENIE%20App-CLICK%20HERE-blue?style=for-the-badge" alt="Launch Z-GENIE">
  </a>
</p>

If the image above doesn’t load, you can also use this plain link:

👉 **https://i5h2oo-angel0i0garza0reyna.shinyapps.io/zgenie_shinyapp/**

- Works in any modern web browser.
- No R, RStudio, or GitHub experience required.
- Note: The current hosting plan has usage limits (free/academic tier). If the app is unavailable or slow, please try again later or run it locally (see below).

---

## 2. Run Z-GENIE Locally (R/RStudio)

If you prefer to run the app on your own machine (for performance, privacy, or reliability), follow these steps.

### 2.1 Install R and RStudio

1. Install **R**: https://cran.r-project.org/  
2. Install **RStudio Desktop**: https://www.rstudio.com/

Open RStudio after installation.

---

### 2.2 Get the Z-GENIE Code

#### Option A (simplest): Download ZIP

1. Go to the GitHub repository:  
   https://github.com/agarzare/Z-GENIE
2. Click the green **“Code”** button.
3. Click **“Download ZIP”**.
4. Unzip the file somewhere easy to find, for example:  
   `Documents/Z-GENIE/`

#### Option B (for Git users): Clone the repository

**HTTPS (most common):**
```bash
git clone https://github.com/agarzare/Z-GENIE.git
```

**SSH (requires an SSH key configured with GitHub):**
```bash
git clone git@github.com:agarzare/Z-GENIE.git
```

---

### 2.3 Install the ZGENIE R Package

In RStudio, install and load `devtools`:

```r
install.packages("devtools")  # run this once
library(devtools)
```

Then install the package from the downloaded/cloned folder.  
Replace the path below with the path on your computer (for example `"~/Documents/Z-GENIE/ZGENIE"`):

```r
devtools::install("/path/to/ZGENIE")
```

Load the package:

```r
library(ZGENIE)
```

---

### 2.4 Launch the Shiny App Locally

Once the package is installed:

```r
library(shiny)
runApp(system.file("shinyapp", package = "ZGENIE"))
```

This should open Z-GENIE in your default web browser, running locally on your machine.

---

## 3. Advanced: Compiling the Z-Hunt Backend

Z-GENIE uses a compiled Z-Hunt executable under the hood. Most users running the packaged app do **not** need to compile this manually. These instructions are for developers or for deploying on a new server/architecture.

### 3.1 Docker-based compilation (Linux x86_64 target)

```bash
docker run --rm --platform linux/amd64   -v ~/zhunt-master:/src   -w /src gcc:10   /bin/bash -c "gcc -march=x86-64 -mtune=generic -o bin/zhunt zhun3.c -lm"
```

- Replace `~/zhunt-master` with the path to your `zhunt-master` directory.
- The resulting binary will be at `bin/zhunt` inside that directory.

### 3.2 Local compilation with gcc (Linux x86_64)

```bash
cd /path/to/zhunt-master
gcc -march=x86-64 -mtune=generic -o bin/zhunt zhun3.c -lm
```

Ensure that `bin/zhunt` is executable and that your deployment knows where to find it (for example, via a configuration option or by adding it to the system `PATH`).

---

## 4. Citing Z-GENIE

If you use Z-GENIE in your work, please cite:

> Garza Reyna A, et al. **Z-GENIE: a user-friendly R/Shiny resource for predicting Z-DNA forming regions in DNA.** *BMC Genomics*. 2025.  
> https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-025-12148-x

---

## 5. Support, Hosting, and Contact

The public Shiny server instance is currently hosted on a limited plan. If you are interested in:

- Collaborating on applications of Z-GENIE, or  
- Supporting more robust hosting (for example, upgrading the shinyapps.io plan or alternative hosting),

please contact the maintainer using the corresponding author information listed in the BMC Genomics article.

For bug reports and feature requests, please open an issue on GitHub:  
https://github.com/agarzare/Z-GENIE/issues
