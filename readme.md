# lipopepID 🔬

**lipopepID** is an R package designed for the automated identification and scoring of lipopeptide compounds from mass spectrometry data. It is specifically optimized for high-resolution mass spectrometers (e.g., Orbitrap Astral).

## Installation

You can install the development version of **lipopepID** from GitHub with:

```R
# install.packages("devtools")
devtools::install_github("honglab-meta/lipopepID")
```

You can also install the development version of **lipopepID** from a local directory with:

```R
devtools::install_local("./lipopepID-main.zip")
```

## Quick start: Testing with test data

### 1. load R package
```R
library(lipopepID)
```

### 2. set file path
```R
rdata_path <- "path/to/LP_Scoring_Workspace.RData"
```

### 3. load parameters
load CONFIG.R file to your R environment.

### 4. doReScore
rdata_path represents your workspace path.
```R
doReScore("LP_Workspace.RData", CONFIG)
```




## From raw data to statistical analysis
```R
library(lipopepID)
```

### 1. Prepare theoretical library
The lib_file is a CSV containing peptide sequences and mass info.
This CSV can be obtained from the authors.
```R
lib <- prepLib("./final_lib.csv", CONFIG)
```

### 2. Feature extraction and identification
Automatically process all mzML files in the target directory to extract features and match them against the library.
```R
processed_list <- doProcess(lib, CONFIG)
```

### 3. Weight optimization
Optimize the scoring parameters based on the identified features to improve confidence and sensitivity.
```R
weight_opt_res <- doWeight(processed_list, lib, CONFIG)
```

### 4. Final scoring & false discovery rate control
Apply the optimized weights to calculate final scores and perform 1% False Discovery Rate (FDR) filtering for the final summary report.
```R
doScore(processed_list, lib, CONFIG, optimized_weights = weight_opt_res$best_weights)
```


