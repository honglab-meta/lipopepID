#' Weight Sensitivity Analysis for Lipopeptide Scoring
#'
#' @description
#' This function performs a sensitivity analysis of the scoring weights using a Monte Carlo
#' simulation approach. It generates random weight combinations to evaluate the robustness
#' of the current configuration and its ability to discriminate between Target and Decoy
#' identifications using AUC (Area Under the Curve).
#'
#' @param scoring_data A data frame containing the individual scores for each metric
#' (Score_Coverage, Score_Continuity, etc.) and a 'db' column indicating "target" or "decoy".
#' @param current_config_weights A named list of the current weights used in CONFIG$weights.
#' @param n_iterations Integer. The number of random weight perturbations to perform. Default is 500.
#' @param seed Integer. Random seed for reproducibility. Default is 42.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{baseline_auc}: The AUC achieved with the current configuration.
#'   \item \code{comparison_table}: A summary table comparing the current config against
#'         equal weights and single-metric scenarios.
#'   \item \code{random_results}: A data frame containing weights and AUCs for all iterations.
#'   \item \code{plot}: A ggplot2 histogram showing the distribution of AUCs vs. the baseline.
#' }
#'
#' @export
#'
#' @importFrom pROC roc auc
#' @importFrom ggplot2 ggplot aes geom_histogram geom_vline annotate labs theme_minimal
#' @importFrom dplyr mutate select filter arrange
#' @importFrom magrittr %>%
#'
#' @examples
#' \dontrun{
#' # Assuming analysis_base and CONFIG$weights are available in the environment:
#' results <- perform_weight_sensitivity_analysis(analysis_base, CONFIG$weights)
#' print(results$comparison_table)
#' plot(results$plot)
#' }
perform_weight_sensitivity_analysis <- function(scoring_data,
                                                current_config_weights,
                                                n_iterations = 500,
                                                seed = 42) {

  # Ensure required namespaces are available without attaching them globally
  if (!requireNamespace("pROC", quietly = TRUE)) stop("Package 'pROC' is required.")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")

  # --- Internal Helper: Calculate AUC for a specific weight vector ---
  evaluate_weights <- function(w_vec, data) {
    # Calculate weighted score
    # Expected order: coverage, continuity, precision, intensity, end_match, adduct
    temp_scores <- (data$Score_Coverage   * w_vec[1]) +
      (data$Score_Continuity * w_vec[2]) +
      (data$Score_Precision  * w_vec[3]) +
      (data$Score_Intensity  * w_vec[4]) +
      (data$Score_End_Match  * w_vec[5]) +
      (data$Score_Adduct     * w_vec[6])

    temp_scores <- pmin(1.0, temp_scores)

    # Calculate AUC using pROC
    roc_obj <- pROC::roc(data$db == "target", temp_scores, quiet = TRUE)
    return(as.numeric(pROC::auc(roc_obj)))
  }

  # --- 1. Setup and Baseline Calculation ---
  weight_names <- c("coverage", "continuity", "precision", "intensity", "end_match", "adduct")
  current_w_vec <- unlist(current_config_weights)[weight_names]

  baseline_auc <- evaluate_weights(current_w_vec, scoring_data)

  # --- 2. Monte Carlo Simulation (Random Perturbations) ---
  set.seed(seed)
  message(paste("Simulating", n_iterations, "random weight combinations..."))

  sensitivity_results <- replicate(n_iterations, {
    # Generate random weights and normalize to the same sum as current weights
    random_w <- runif(6)
    random_w <- (random_w / sum(random_w)) * sum(current_w_vec)

    auc_val <- evaluate_weights(random_w, scoring_data)
    return(c(random_w, auc = auc_val))
  }) %>% t() %>% as.data.frame()

  colnames(sensitivity_results) <- c(weight_names, "AUC")

  # --- 3. Scenario Comparison Table ---
  message("Generating comparison scenarios...")
  scenarios <- list(
    "Current Config"  = current_w_vec,
    "Equal Weights"   = rep(1/6, 6),
    "Coverage Only"   = c(1,0,0,0,0,0),
    "Continuity Only" = c(0,1,0,0,0,0),
    "Precision Only"  = c(0,0,1,0,0,0),
    "Intensity Only"  = c(0,0,0,1,0,0),
    "End_Match Only"  = c(0,0,0,0,1,0),
    "Adduct Only"     = c(0,0,0,0,0,1)
  )

  comp_df <- data.frame(
    Scenario = names(scenarios),
    AUC = sapply(scenarios, evaluate_weights, data = scoring_data),
    stringsAsFactors = FALSE
  )

  # --- 4. Visualization ---
  p <- ggplot2::ggplot(sensitivity_results, ggplot2::aes(x = AUC)) +
    ggplot2::geom_histogram(bins = 30, fill = "skyblue", color = "white", alpha = 0.7) +
    ggplot2::geom_vline(xintercept = baseline_auc, color = "red", linetype = "dashed", size = 1) +
    ggplot2::annotate("text", x = baseline_auc, y = n_iterations/100,
                      label = paste("Current (", round(baseline_auc, 4), ")", sep=""),
                      color = "red", angle = 90, vjust = -0.5) +
    ggplot2::labs(title = "Weight Sensitivity Analysis",
                  subtitle = paste("Current Config vs.", n_iterations, "Random Samples"),
                  x = "Discrimination Power (AUC)",
                  y = "Frequency") +
    ggplot2::theme_minimal()

  print(head(sensitivity_results[order(-sensitivity_results$AUC), ], 5))

  # --- 5. Return List ---
  return(list(
    baseline_auc = baseline_auc,
    comparison_table = comp_df,
    random_results = sensitivity_results,
    plot = p
  ))
}


