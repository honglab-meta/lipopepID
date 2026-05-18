.onAttach <- function(libname, pkgname) {
  # List of core dependency packages to be loaded into the global environment
  dependencies <- c("dplyr", "purrr", "stringr", "MSnbase", "tidyr", "tibble", "readr",
                    "tidyverse","stringi","broom","magrittr")

  # Attempt to load each package silently
  for (pkg in dependencies) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      packageStartupMessage(paste("Dependency missing. Attempting to install:", pkg))
      install.packages(pkg)
    }
    # Force attach to search path to ensure functions like bind_rows are available
    library(pkg, character.only = TRUE, warn.conflicts = FALSE)
  }

  packageStartupMessage(">>> lipopepID loaded successfully. Environment initialized.")
}
