#' Prepare Theoretical Lipopeptide Library
#'
#' This function processes a CSV file containing theoretical lipopeptide data,
#' expands it based on adduct offsets, handles unsaturation, and converts the
#' data into a format compatible with MSnbase Spectrum2 objects.
#'
#' @param lib_path Character. Path to the theoretical library CSV file.
#' @param CONFIG List. Global configuration containing adduct offsets and scoring parameters.
#' @export
prepLib <- function(lib_path = "./final_lib.csv", CONFIG) {
  cat("--- preparing ---\n")

  df_theory <- read.csv(lib_path, row.names = 1)
  colnames(df_theory)[9] <- "compound_name"
  colnames(df_theory)[1] <- "mz_pre"
  colnames(df_theory)[8] <- "mz"

  unsat_offset <- -2.0156

  df_theory_expanded <- purrr::map_df(names(CONFIG$adduct_offset), function(adduct_name) {
    offset <- CONFIG$adduct_offset[[adduct_name]]
    df_theory %>%
      dplyr::mutate(
        mz_pre = as.numeric(mz_pre) + offset,
        mz = as.numeric(mz) + offset,
        adduct_label = adduct_name
      )
  })

  df_theory_expanded <- df_theory_expanded %>%
    dplyr::bind_rows(
      df_theory_expanded %>%
        dplyr::filter(grepl("LF|CF", compound_name)) %>%
        dplyr::mutate(
          mz_pre = mz_pre + unsat_offset,
          mz = mz + unsat_offset,
          compound_name = paste0(compound_name, "_unsat")
        )
    )

  df_theory <- df_theory_expanded[!is.na(df_theory_expanded$mz), ]
  df_theory$intensity <- 100

  theory_list <- split(df_theory, df_theory$compound_name)
  spec_library <- lapply(seq_along(theory_list), function(i) {
    data <- theory_list[[i]][order(theory_list[[i]]$mz), ]
    if(nrow(data) == 0) return(NULL)
    new("Spectrum2", mz = as.numeric(data$mz), intensity = as.numeric(data$intensity),
        msLevel = as.integer(2), centroided = TRUE, precScanNum = as.integer(i))
  })
  spec_library <- spec_library[!sapply(spec_library, is.null)]
  names(spec_library) <- names(theory_list)

  theory_ion_counts <- df_theory %>%
    dplyr::filter(ion_type %in% c("b", "y")) %>%
    dplyr::group_by(compound_name) %>%
    dplyr::summarise(Total_Theoretical_by = n_distinct(paste0(ion_type, ion_num)), .groups = 'drop')

  cat("Done\n")
  return(list(spec_library = spec_library, theory_ion_counts = theory_ion_counts, df_theory = df_theory))
}


#' Process Mass Spectrometry Data for Lipopeptides
#'
#' @description
#' First stage of the pipeline. Reads mzML files, performs peak picking and cleaning,
#' matches against the theoretical library, annotates fragments, and writes the MS2 detail CSV.
#'
#' @param lib List. The output from the \code{prepLib} function.
#' @param CONFIG List. Global configuration containing mass tolerances.
#'
#' @return A nested list containing processed data for each file, required for down-stream scoring.
#' @export
doProcess <- function(lib, CONFIG) {
  require(stats)
  spec_library <- lib$spec_library
  df_theory <- lib$df_theory

  mzML_files <- list.files(pattern = "\\.mzML$", full.names = FALSE)
  n_files <- length(mzML_files)
  tol_val <- CONFIG$ms2_tolerance_ppm * 1e-6

  cat(sprintf("\n--- Starting doProcess Workflow. Found %d files. ---\n", n_files))

  process_results <- list()

  for (f_idx in seq_along(mzML_files)) {
    file_name <- mzML_files[f_idx]
    file_prefix <- stringr::str_remove(file_name, "\\.mzML$")
    cat(sprintf("\n[%d/%d] processing: %s\n", f_idx, n_files, file_name))

    raw_exp <- MSnbase::readMSData(file_name, msLevel = 2, mode = "onDisk")
    exp_specs <- MSnbase::spectra(raw_exp)
    names(exp_specs) <- paste0("Scan_", MSnbase::fData(raw_exp)$acquisitionNum)

    exp_specs_clean <- lapply(exp_specs, function(s) {
      s <- MSnbase::pickPeaks(s) %>% clean(all = TRUE)
      max_int <- max(MSnbase::intensity(s), na.rm = TRUE)
      if(max_int > 0) s <- removePeaks(s, t = max_int * 0.001)
      return(s)
    })

    common_peaks_matrix <- matrix(NA, nrow = length(exp_specs_clean), ncol = length(spec_library),
                                  dimnames = list(names(exp_specs_clean), names(spec_library)))

    cat("    -> calculating common peaks matrix: \n")
    pb <- txtProgressBar(min = 0, max = length(exp_specs_clean), style = 3, width = 50)
    for (i in seq_along(exp_specs_clean)) {
      s_i <- exp_specs_clean[[i]]
      for (j in seq_along(spec_library)) {
        common_peaks_matrix[i, j] <- MSnbase::compareSpectra(s_i, spec_library[[j]],
                                                             fun = "common", relative = TRUE, tolerance = tol_val)
      }
      setTxtProgressBar(pb, i)
    }
    close(pb)

    matched_results <- as.data.frame(common_peaks_matrix) %>%
      tibble::rownames_to_column("Experimental_Scan") %>%
      tidyr::pivot_longer(cols = -Experimental_Scan, names_to = "Theoretical_Lipopeptide", values_to = "Common_Peaks") %>%
      dplyr::filter(Common_Peaks > 7)

    if(nrow(matched_results) == 0) {
      cat("    [Skip] No matches found for this file\n")
      next
    }

    exp_metadata <- MSnbase::fData(raw_exp) %>% tibble::rownames_to_column("Scan_ID") %>%
      dplyr::mutate(Experimental_Scan = paste0("Scan_", acquisitionNum),
                    precursorMZ = round(as.numeric(precursorMZ), 4)) %>%
      dplyr::select(Experimental_Scan, precursorMZ)

    matched_results <- matched_results %>% dplyr::left_join(exp_metadata, by = "Experimental_Scan")

    filtered_mz_results <- matched_results %>%
      dplyr::mutate(base_theoretical_mz = as.numeric(sub("^([0-9.]+)_.*", "\\1", Theoretical_Lipopeptide))) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        errors_ppm = list(abs(precursorMZ - (base_theoretical_mz + CONFIG$adduct_offset)) /
                            (base_theoretical_mz + CONFIG$adduct_offset) * 1e6),
        matched_adduct = list(names(CONFIG$adduct_offset)[unlist(errors_ppm) < CONFIG$ms1_match_ppm])
      ) %>%
      dplyr::filter(length(matched_adduct) > 0) %>%
      dplyr::mutate(matched_adduct = paste(unlist(matched_adduct), collapse = "/")) %>%
      dplyr::ungroup() %>% dplyr::filter(Common_Peaks > 7)

    final_identified_lps <- filtered_mz_results %>%
      dplyr::group_by(Theoretical_Lipopeptide) %>% dplyr::slice_max(Common_Peaks, n = 1, with_ties = FALSE) %>% dplyr::ungroup()

    final_identified_lps_output <- final_identified_lps %>%
      dplyr::mutate(parts = stringr::str_split(Theoretical_Lipopeptide, "_"),
                    LP_Type = sapply(parts, function(x) if (length(x) >= 4 && (x[4] == "CF" || x[4] == "LF")) x[4] else x[3]),
                    AA_Sequence = sapply(parts, function(x) x[2]),
                    FA_Details = sapply(parts, function(x) {
                      if (length(x) >= 4 && (x[4] == "CF" || x[4] == "LF")) {
                        valid_idx <- c(3, 5, 6); return(paste(x[valid_idx[valid_idx <= length(x)]], collapse = "_"))
                      } else {
                        if (length(x) >= 5 && (x[3] == "CF" || x[3] == "LF")) paste(x[4], x[5], sep = "_") else x[4]
                      }
                    })) %>%
      dplyr::left_join(df_theory %>% dplyr::select(compound_name, db) %>% distinct(), by = c("Theoretical_Lipopeptide" = "compound_name")) %>%
      dplyr::select(LP_Type, FA_Details, AA_Sequence, everything(), -parts)

    detailed_matches_list <- list()
    for (i in 1:nrow(final_identified_lps)) {
      row_data <- final_identified_lps[i, ]; exp_id <- row_data$Experimental_Scan; theo_id <- row_data$Theoretical_Lipopeptide
      s_exp <- exp_specs_clean[[exp_id]]; s_theo <- spec_library[[theo_id]]
      mz_exp <- MSnbase::mz(s_exp); mz_theo <- MSnbase::mz(s_theo)
      matches <- lapply(mz_exp, function(m) {
        diffs <- abs(mz_theo - m) / m
        best <- which.min(diffs)
        if (length(best) > 0 && diffs[best] <= CONFIG$ms2_tolerance_ppm * 1e-6) {
          return(data.frame(Matched_mz_exp = m, Matched_mz_theo = mz_theo[best],
                            Intensity_exp = MSnbase::intensity(s_exp)[which(mz_exp == m)[1]], ppm_error = diffs[best] * 1e6))
        }
        return(NULL)
      })
      curr_df <- do.call(rbind, matches)
      if (!is.null(curr_df)) {
        curr_df$Experimental_Scan <- exp_id; curr_df$Theoretical_Lipopeptide <- theo_id
        detailed_matches_list[[i]] <- curr_df
      }
    }

    if (length(detailed_matches_list) == 0) next

    matched_ions_detail <- do.call(rbind, detailed_matches_list)
    lookup_table <- df_theory %>%
      dplyr::mutate(
        ion_label = paste0(as.character(ion_type), as.character(ion_num)),
        add_label = as.character(adduct_label)
      ) %>%
      dplyr::select(
        Theoretical_Lipopeptide = compound_name,
        Matched_mz_theo = mz,
        ion_label,
        add_label
      ) %>%
      dplyr::distinct(Theoretical_Lipopeptide, Matched_mz_theo, .keep_all = TRUE)

    feature_lookup <- final_identified_lps_output %>%
      dplyr::select(Experimental_Scan, Theoretical_Lipopeptide, LP_Type, FA_Details, AA_Sequence, matched_adduct, db)

    final_output_annotated <- matched_ions_detail %>%
      dplyr::left_join(exp_metadata, by = "Experimental_Scan") %>%
      dplyr::left_join(lookup_table, by = c("Theoretical_Lipopeptide", "Matched_mz_theo")) %>%
      dplyr::left_join(feature_lookup, by = c("Experimental_Scan", "Theoretical_Lipopeptide")) %>%
      dplyr::rename(Ion_Annotation = ion_label) %>%
      dplyr::select(Experimental_Scan, LP_Type, precursorMZ, Theoretical_Lipopeptide, matched_adduct,
                    Matched_mz_exp, Matched_mz_theo, ppm_error, Intensity_exp, add_label, Ion_Annotation, db)


    final_output_annotated_out <- final_output_annotated %>% dplyr::filter(db == "target") %>%
      dplyr::select(Experimental_Scan, LP_Type, precursorMZ, Theoretical_Lipopeptide, matched_adduct,
                    Matched_mz_exp, Matched_mz_theo, ppm_error, Intensity_exp, add_label, Ion_Annotation)

    write.csv(final_output_annotated_out, paste0(file_prefix, "_MS2_detail.csv"), row.names = FALSE)


    process_results[[file_name]] <- list(
      final_output_annotated = final_output_annotated,
      exp_specs_clean = exp_specs_clean,
      file_prefix = file_prefix
    )
  }

  cat("\n--- doProcess Completed. ---\n")
  return(process_results)
}


#' Optimize Weights Using Monte Carlo Simulation on Processed Lipopeptide Data
#'
#' @description
#' Second stage of the pipeline. Dynamically calculates diagnostic scores (Coverage, Continuity,
#' Precision, etc.) for both target and decoy hits, then runs a Monte Carlo simulation to
#' find the optimal weight combination that maximizes discrimination power (AUC).
#'
#' @param process_list List. The output object returned by \code{doProcess}.
#' @param lib List. The output from the \code{prepLib} function.
#' @param CONFIG List. Global configuration containing tolerances and baseline weights.
#' @param n_iterations Integer. Number of Monte Carlo random weight samples. Default is 500.
#' @param seed Integer. Random seed for reproducibility. Default is 42.
#'
#' @return A list containing the optimization baseline AUC, the best weights discovered, and the full random trial data.
#' @export
doWeight <- function(process_list, lib, CONFIG, n_iterations = 500, seed = 42) {
  if (!requireNamespace("pROC", quietly = TRUE)) stop("Package 'pROC' is required.")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")

  theory_ion_counts <- lib$theory_ion_counts
  all_file_features <- list()

  cat("\n--- Starting doWeight: Extracting Scoring Features from Processed Data... ---\n")


  for (file_name in names(process_list)) {
    cache <- process_list[[file_name]]
    final_output_annotated <- cache$final_output_annotated
    exp_specs_clean <- cache$exp_specs_clean

    scoring_details <- final_output_annotated %>%
      dplyr::left_join(theory_ion_counts, by = c("Theoretical_Lipopeptide" = "compound_name")) %>%
      dplyr::group_by(Experimental_Scan, Theoretical_Lipopeptide) %>%
      dplyr::summarise(
        db = dplyr::first(db),
        Score_Precision = pmax(0, 1 - (stats::median(abs(ppm_error), na.rm = TRUE) / CONFIG$ms2_tolerance_ppm)),
        Score_Continuity = {
          b_nums <- sort(as.numeric(stringr::str_extract(Ion_Annotation[grepl("^b", Ion_Annotation)], "\\d+")))
          y_nums <- sort(as.numeric(stringr::str_extract(Ion_Annotation[grepl("^y", Ion_Annotation)], "\\d+")))
          count_cont <- function(nums) if(length(nums) < 2) 0 else sum(diff(nums) == 1)
          total_cont_pairs <- count_cont(b_nums) + count_cont(y_nums)
          theoretical_max_pairs <- pmax(1, dplyr::first(Total_Theoretical_by) - 2)
          pmin(1.0, total_cont_pairs / theoretical_max_pairs)
        },
        Score_Intensity = {
          all_intensities <- MSnbase::intensity(exp_specs_clean[[dplyr::first(Experimental_Scan)]])
          base_peak_int <- max(all_intensities, na.rm = TRUE)
          significant_mask <- all_intensities > (base_peak_int * CONFIG$int_base_peak_ratio)
          significant_peaks_sum <- sum(all_intensities[significant_mask], na.rm = TRUE)
          matched_significant_int <- sum(Intensity_exp[Intensity_exp > (base_peak_int * CONFIG$int_base_peak_ratio)], na.rm = TRUE)
          tic_ratio <- if(significant_peaks_sum > 0) pmin(1.0, matched_significant_int / significant_peaks_sum) else 0
          avg_rel_int <- mean(Intensity_exp / base_peak_int, na.rm = TRUE)
          intensity_penalty <- dplyr::if_else(avg_rel_int > CONFIG$int_avg_rel_threshold, 1.0, (avg_rel_int / CONFIG$int_penalty_coef))
          pmin(1.0, tic_ratio * intensity_penalty)
        },
        Matched_Unique_Ions = dplyr::n_distinct(Ion_Annotation),
        Total_Theory_Ions = dplyr::first(Total_Theoretical_by),
        Score_Coverage = dplyr::if_else(Matched_Unique_Ions >= CONFIG$cov_min_unique_ions,
                                        Matched_Unique_Ions / Total_Theory_Ions,
                                        (Matched_Unique_Ions / Total_Theory_Ions) * CONFIG$cov_penalty_factor),
        Score_End_Match = {
          current_type <- dplyr::first(LP_Type)
          target_y <- if(current_type %in% c("LF", "CF")) "y10" else "y7"
          has_b1 = any(Ion_Annotation == "b1", na.rm = TRUE); has_target_y = any(Ion_Annotation == target_y, na.rm = TRUE)
          dplyr::case_when(has_b1 & has_target_y ~ 1.0, has_b1 | has_target_y ~ 0.3, TRUE ~ 0.0)
        },
        Score_Adduct = {
          curr_adduct <- as.character(dplyr::first(matched_adduct))
          val <- CONFIG$adduct_scores[curr_adduct]
          if(is.na(val)) CONFIG$adduct_scores["Default"] else val
        }, .groups = 'drop'
      )
    all_file_features[[file_name]] <- scoring_details
  }


  scoring_data <- do.call(rbind, all_file_features)

  if(is.null(scoring_data) || nrow(scoring_data) == 0) {
    stop("No score matrices extracted. Please verify if matches exist in doProcess output.")
  }


  evaluate_weights <- function(w_vec, data) {
    temp_scores <- (data$Score_Coverage   * w_vec[1]) +
      (data$Score_Continuity * w_vec[2]) +
      (data$Score_Precision  * w_vec[3]) +
      (data$Score_Intensity  * w_vec[4]) +
      (data$Score_End_Match  * w_vec[5]) +
      (data$Score_Adduct     * w_vec[6])

    temp_scores <- pmin(1.0, temp_scores)
    roc_obj <- pROC::roc(data$db == "target", temp_scores, quiet = TRUE)
    return(as.numeric(pROC::auc(roc_obj)))
  }


  weight_names <- c("coverage", "continuity", "precision", "intensity", "end_match", "adduct")
  current_w_vec <- unlist(CONFIG$weights)[weight_names]
  baseline_auc <- evaluate_weights(current_w_vec, scoring_data)


  if (!is.null(seed)) set.seed(seed)
  message(paste("Simulating", n_iterations, "random weight combinations..."))

  sensitivity_results <- replicate(n_iterations, {
    random_w <- runif(6)
    random_w <- (random_w / sum(random_w)) * sum(current_w_vec)
    auc_val <- evaluate_weights(random_w, scoring_data)
    return(c(random_w, auc = auc_val))
  }) %>% t() %>% as.data.frame()

  colnames(sensitivity_results) <- c(weight_names, "AUC")


  best_row <- sensitivity_results[which.max(sensitivity_results$AUC), ]
  best_weights <- unlist(best_row[weight_names])
  best_auc <- best_row$AUC

  message(sprintf("\nOptimization complete.\nBaseline AUC: %s\nMaximized AUC: %s",
                  round(baseline_auc, 4), round(best_auc, 4)))


  message("\nTop 5 optimized weight configurations discovered:")
  print(head(sensitivity_results[order(-sensitivity_results$AUC), ], 5))

  return(list(
    baseline_auc = baseline_auc,
    best_auc = best_auc,
    best_weights = as.list(best_weights),
    random_results = sensitivity_results
  ))
}


#' Scoring and FDR Calculation with Memory-based Processed Data
#'
#' @param process_list List. The output object returned by doProcess.
#' @param lib List. The output from the prepLib function.
#' @param CONFIG List. Global configuration.
#' @param optimized_weights List or Vector. Optional optimized weights from doWeight.
#' @export
doScore <- function(process_list, lib, CONFIG, optimized_weights = NULL) {

  cat("\n--- Starting doScore Workflow ---\n")


  theory_ion_counts <- lib$theory_ion_counts
  df_theory <- lib$df_theory


  target_weights <- CONFIG$weights
  if (!is.null(optimized_weights)) {
    cat("    -> Applying optimized weights discovered by Monte Carlo simulation.\n")
    target_weights <- as.list(optimized_weights)
    CONFIG$weights <- target_weights
  }


  for (file_name in names(process_list)) {
    cache <- process_list[[file_name]]
    final_output_annotated <- cache$final_output_annotated
    exp_specs_clean <- cache$exp_specs_clean
    file_prefix <- cache$file_prefix

    cat(sprintf("Processing Scores for: %s\n", file_prefix))


    scoring_details <- final_output_annotated %>%
      dplyr::left_join(theory_ion_counts, by = c("Theoretical_Lipopeptide" = "compound_name")) %>%
      dplyr::group_by(Experimental_Scan, Theoretical_Lipopeptide) %>%
      dplyr::summarise(
        db = dplyr::first(db),
        LP_Type = dplyr::first(LP_Type),
        matched_adduct = dplyr::first(matched_adduct),

        Score_Precision = pmax(0, 1 - (stats::median(abs(ppm_error), na.rm = TRUE) / CONFIG$ms2_tolerance_ppm)),

        Score_Continuity = {
          b_nums <- sort(as.numeric(stringr::str_extract(Ion_Annotation[grepl("^b", Ion_Annotation)], "\\d+")))
          y_nums <- sort(as.numeric(stringr::str_extract(Ion_Annotation[grepl("^y", Ion_Annotation)], "\\d+")))
          count_cont <- function(nums) if(length(nums) < 2) 0 else sum(diff(nums) == 1)
          total_cont_pairs <- count_cont(b_nums) + count_cont(y_nums)
          theoretical_max_pairs <- pmax(1, dplyr::first(Total_Theoretical_by) - 2)
          pmin(1.0, total_cont_pairs / theoretical_max_pairs)
        },

        Score_Intensity = {
          s <- exp_specs_clean[[dplyr::first(Experimental_Scan)]]
          base_peak_int <- max(MSnbase::intensity(s), na.rm = TRUE)
          significant_peaks_sum <- sum(MSnbase::intensity(s)[MSnbase::intensity(s) > (base_peak_int * CONFIG$int_base_peak_ratio)], na.rm = TRUE)
          matched_significant_int <- sum(Intensity_exp[Intensity_exp > (base_peak_int * CONFIG$int_base_peak_ratio)], na.rm = TRUE)
          tic_ratio <- if(significant_peaks_sum > 0) pmin(1.0, matched_significant_int / significant_peaks_sum) else 0
          avg_rel_int <- mean(Intensity_exp / base_peak_int, na.rm = TRUE)
          intensity_penalty <- dplyr::if_else(avg_rel_int > CONFIG$int_avg_rel_threshold, 1.0, (avg_rel_int / CONFIG$int_penalty_coef))
          pmin(1.0, tic_ratio * intensity_penalty)
        },
        Matched_Unique_Ions = dplyr::n_distinct(Ion_Annotation),
        Total_Theory_Ions = dplyr::first(Total_Theoretical_by),

        Score_Coverage = dplyr::if_else(Matched_Unique_Ions >= CONFIG$cov_min_unique_ions,
                                        Matched_Unique_Ions / Total_Theory_Ions,
                                        (Matched_Unique_Ions / Total_Theory_Ions) * CONFIG$cov_penalty_factor),

        Score_End_Match = {
          current_type <- dplyr::first(LP_Type)
          target_y <- if(current_type %in% c("LF", "CF")) "y10" else "y7"
          has_b1 = any(Ion_Annotation == "b1", na.rm = TRUE); has_target_y = any(Ion_Annotation == target_y, na.rm = TRUE)
          dplyr::case_when(has_b1 & has_target_y ~ 1.0, has_b1 | has_target_y ~ 0.3, TRUE ~ 0.0)
        },

        Score_Adduct = {
          curr_adduct <- as.character(dplyr::first(matched_adduct))
          val <- CONFIG$adduct_scores[curr_adduct]
          if(is.na(val)) CONFIG$adduct_scores["Default"] else val
        }, .groups = 'drop'
      )

    final_report <- scoring_details %>%
      dplyr::mutate(LipopepID_Score =
                      (Score_Coverage   * target_weights$coverage) +
                      (Score_Continuity * target_weights$continuity) +
                      (Score_Precision  * target_weights$precision) +
                      (Score_Intensity  * target_weights$intensity) +
                      (Score_End_Match  * target_weights$end_match) +
                      (Score_Adduct     * target_weights$adduct)) %>%
      dplyr::mutate(LipopepID_Score = pmin(1.0, round(LipopepID_Score, 4))) %>%
      dplyr::arrange(dplyr::desc(LipopepID_Score)) %>%
      dplyr::mutate(cum_target = cumsum(db == "target"),
                    cum_decoy = cumsum(db == "decoy"),
                    current_FDR = (2 * cum_decoy) / pmax(1, (cum_target + cum_decoy)))

    cutoff_1 <- final_report %>% dplyr::filter(current_FDR <= 0.01) %>% dplyr::slice_tail(n = 1) %>% dplyr::pull(LipopepID_Score)
    cutoff_1 <- if(length(cutoff_1) == 0) max(final_report$LipopepID_Score) else cutoff_1

    final_report_cleaned <- final_report %>%
      dplyr::filter(!(LP_Type == "CS" & matched_adduct == "H2O_add")) %>%
      dplyr::filter(!(LP_Type == "LS" & matched_adduct == "H2O_loss")) %>%
      dplyr::filter(!(LP_Type == "CF" & matched_adduct == "H2O_add")) %>%
      dplyr::filter(!(LP_Type == "LF" & matched_adduct == "H2O_loss"))

    write.csv(final_report_cleaned %>% dplyr::filter(db == "target"),
              paste0(file_prefix, "_MS1_Results.csv"), row.names = FALSE)

    final_conf <- final_report_cleaned %>%
      dplyr::filter(db == "target") %>%
      dplyr::filter(LipopepID_Score >= cutoff_1)

    write.csv(final_conf, paste0(file_prefix, "_MS1_results_1%FDR.csv"), row.names = FALSE)

    # intensity_calculation_support
    intensity_calculation_support <- purrr::map_df(names(exp_specs_clean), function(scan_id) {
      s <- exp_specs_clean[[scan_id]]
      base_int <- max(MSnbase::intensity(s), na.rm = TRUE)
      sig_tic <- sum(MSnbase::intensity(s)[MSnbase::intensity(s) > (base_int * CONFIG$int_base_peak_ratio)], na.rm = TRUE)
      data.frame(Experimental_Scan = scan_id, Base_Peak_Intensity = base_int, TIC_Significant_Ref = sig_tic)
    })

    final_identified_lps_output <- final_report
    rdata_file <- paste0(file_prefix, "_Scoring_Workspace.RData")

    save(
      final_output_annotated,
      final_identified_lps_output,
      theory_ion_counts,
      intensity_calculation_support,
      CONFIG,
      file_prefix,
      file = rdata_file
    )

    cat(sprintf("   [OK] %s: Identifications at 1%% FDR: %d. Workspace saved.\n", file_prefix, nrow(final_conf)))
  }
}


#' Rescoring based on saved workspace RData
#'
#' @param rdata_path Character. Path to the _Scoring_Workspace.RData file.
#' @param CONFIG List. Global configuration containing new weights.
#' @export
doReScore <- function(rdata_path, CONFIG) {

  if (!file.exists(rdata_path)) stop("no RData file: ", rdata_path)

  tmp_env <- new.env()
  load(rdata_path, envir = tmp_env)

  file_prefix <- tmp_env$file_prefix
  annotated_data <- tmp_env$final_output_annotated
  theory_counts  <- tmp_env$theory_ion_counts
  int_support    <- tmp_env$intensity_calculation_support

  cat(sprintf("\n--- Rescoring: %s ---\n", file_prefix))

  cat("    -> Recalculating scoring features from annotated data...\n")

  scoring_details <- annotated_data %>%

    dplyr::left_join(theory_counts, by = c("Theoretical_Lipopeptide" = "compound_name")) %>%
    dplyr::group_by(Experimental_Scan, Theoretical_Lipopeptide) %>%
    dplyr::summarise(
      db = dplyr::first(db),
      LP_Type = dplyr::first(LP_Type),
      matched_adduct = dplyr::first(matched_adduct),

      # [1] Precision Score
      Score_Precision = pmax(0, 1 - (stats::median(abs(ppm_error), na.rm = TRUE) / CONFIG$ms2_tolerance_ppm)),

      # [2] Continuity Score
      Score_Continuity = {
        b_nums <- sort(as.numeric(stringr::str_extract(Ion_Annotation[grepl("^b", Ion_Annotation)], "\\d+")))
        y_nums <- sort(as.numeric(stringr::str_extract(Ion_Annotation[grepl("^y", Ion_Annotation)], "\\d+")))
        count_cont <- function(nums) if(length(nums) < 2) 0 else sum(diff(nums) == 1)
        total_cont_pairs <- count_cont(b_nums) + count_cont(y_nums)
        theoretical_max_pairs <- pmax(1, dplyr::first(Total_Theoretical_by) - 2)
        pmin(1.0, total_cont_pairs / theoretical_max_pairs)
      },

      # [3] Intensity Score
      Score_Intensity = {

        scan_info <- int_support[int_support$Experimental_Scan == dplyr::first(Experimental_Scan), ]
        if(nrow(scan_info) == 0) {
          0
        } else {
          base_peak_int <- scan_info$Base_Peak_Intensity
          sig_tic_ref <- scan_info$TIC_Significant_Ref

          matched_significant_int <- sum(Intensity_exp[Intensity_exp > (base_peak_int * CONFIG$int_base_peak_ratio)], na.rm = TRUE)
          tic_ratio <- if(sig_tic_ref > 0) pmin(1.0, matched_significant_int / sig_tic_ref) else 0

          avg_rel_int <- mean(Intensity_exp / base_peak_int, na.rm = TRUE)
          intensity_penalty <- dplyr::if_else(avg_rel_int > CONFIG$int_avg_rel_threshold, 1.0, (avg_rel_int / CONFIG$int_penalty_coef))
          pmin(1.0, tic_ratio * intensity_penalty)
        }
      },

      # [4] Coverage Score
      Matched_Unique_Ions = dplyr::n_distinct(Ion_Annotation),
      Total_Theory_Ions = dplyr::first(Total_Theoretical_by),
      Score_Coverage = dplyr::if_else(Matched_Unique_Ions >= CONFIG$cov_min_unique_ions,
                                      Matched_Unique_Ions / Total_Theory_Ions,
                                      (Matched_Unique_Ions / Total_Theory_Ions) * CONFIG$cov_penalty_factor),

      # [5] End Match Score
      Score_End_Match = {
        current_type <- dplyr::first(LP_Type)
        target_y <- if(current_type %in% c("LF", "CF")) "y10" else "y7"
        has_b1 = any(Ion_Annotation == "b1", na.rm = TRUE)
        has_target_y = any(Ion_Annotation == target_y, na.rm = TRUE)
        dplyr::case_when(has_b1 & has_target_y ~ 1.0, has_b1 | has_target_y ~ 0.3, TRUE ~ 0.0)
      },

      # [6] Adduct Score
      Score_Adduct = {
        curr_adduct <- as.character(dplyr::first(matched_adduct))
        val <- CONFIG$adduct_scores[curr_adduct]
        if(is.na(val)) CONFIG$adduct_scores["Default"] else val
      }, .groups = 'drop'
    )

  target_weights <- CONFIG$weights
  cat("    -> Applying weight configuration and calculating total scores...\n")

  final_report <- scoring_details %>%
    dplyr::mutate(LipopepID_Score =
                    (Score_Coverage   * target_weights$coverage) +
                    (Score_Continuity * target_weights$continuity) +
                    (Score_Precision  * target_weights$precision) +
                    (Score_Intensity  * target_weights$intensity) +
                    (Score_End_Match  * target_weights$end_match) +
                    (Score_Adduct     * target_weights$adduct)) %>%
    dplyr::mutate(LipopepID_Score = pmin(1.0, round(LipopepID_Score, 4))) %>%
    dplyr::arrange(dplyr::desc(LipopepID_Score))

  final_report <- final_report %>%
    dplyr::mutate(cum_target = cumsum(db == "target"),
                  cum_decoy = cumsum(db == "decoy"),
                  current_FDR = (2 * cum_decoy) / pmax(1, (cum_target + cum_decoy)))

  cutoff_1 <- final_report %>%
    dplyr::filter(current_FDR <= 0.01) %>%
    dplyr::slice_tail(n = 1) %>%
    dplyr::pull(LipopepID_Score)

  cutoff_1 <- if(length(cutoff_1) == 0) max(final_report$LipopepID_Score) else cutoff_1

  write.csv(final_report %>% dplyr::filter(db == "target"),
            paste0(file_prefix, "_MS1_Results_Rescored.csv"), row.names = FALSE)

  final_conf <- final_report %>% dplyr::filter(db == "target" & LipopepID_Score >= cutoff_1)
  write.csv(final_conf, paste0(file_prefix, "_MS1_results_1%FDR_Rescored.csv"), row.names = FALSE)

  cat(sprintf("   [OK] Rescoring complete. High confidence hits (1%% FDR): %d\n", nrow(final_conf)))
}
