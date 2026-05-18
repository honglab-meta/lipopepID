# lipopepID 🔬

**lipopepID** is an R package designed for the automated identification and scoring of lipopeptide compounds from mass spectrometry data. It is specifically optimized for high-resolution mass spectrometers (e.g., Orbitrap Astral).

## Installation

You can install the development version of **lipopepID** from GitHub with:

```R
# install.packages("MSnbase")
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("MSnbase")
# install.packages("devtools")
devtools::install_github("honglab-meta/lipopepID")
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
load CONFIG.R file to your R environment

### 4. rerun doScore
rdata_path represents your workspace path
```R
doScore(rdata_path,CONFIG)
```

### 5. rerun doStat
```R
doStat()
```



## From raw data to statistical analysis
```R
library(lipopepID)
```

### 1. Prepare Theoretical Library
The lib_file is a CSV containing peptide sequences and mass info.
This CSV can be obtained from the authors.
```R
lib <- prepLib("path/to/your_library.csv", CONFIG)
```

### 2. Execute Batch Processing
Automatically searches all mzML files in the directory.
```R
doMain(lib, CONFIG)
```

### 3. Scoring & FDR Filtering
Apply 1% False Discovery Rate (FDR) thresholds.
Note: Ensure your RData workspace is in the current directory.
```R
doScore(rdata_path,CONFIG)
```

### 4. Statistical Aggregation
Generate final summary reports across all samples.
```R
doStat()
```


