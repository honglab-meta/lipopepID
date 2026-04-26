# lipopepID 🔬

**lipopepID** is an R package designed for the automated identification and scoring of lipopeptide compounds from mass spectrometry data. It is specifically optimized for high-resolution mass spectrometers (e.g., Orbitrap Astral).

## Installation

You can install the development version of **lipopepID** from GitHub with:

```R
# install.packages("devtools")
devtools::install_github("honglab/lipopepID")
```

# Quick start
library(lipopepID)

# 1. Prepare Theoretical Library
# The lib_file is a CSV containing peptide sequences and mass info.
# This CSV can be obtained from the authors.
lib <- prepLib("path/to/your_library.csv", CONFIG)

# 2. Execute Batch Processing
# Automatically searches all mzML files in the directory.
doMain(lib, CONFIG)

# 3. Scoring & FDR Filtering
# Apply 1% False Discovery Rate (FDR) thresholds.
# Note: Ensure your RData workspace is in the current directory.
doScore(CONFIG)

# 4. Statistical Aggregation
# Generate final summary reports across all samples.
doStat()


# Adjust config parameters and rerun doScore and doStat in the RData file.

# 1. load R package
library(lipopepID)

# 2. set file path
rdata_path <- "path/to/LP_Scoring_Workspace.RData"

# 3. load parameters
# load CONFIG.R file to your R environment

# 4. rerun doScore
doScore(rdata_path,CONFIG)

# 5. rerun doStat
doStat()


