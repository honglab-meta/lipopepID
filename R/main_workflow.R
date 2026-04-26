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

  df_theory_expanded <- map_df(names(CONFIG$adduct_offset), function(adduct_name) {
    offset <- CONFIG$adduct_offset[[adduct_name]]
    df_theory %>%
      mutate(
        mz_pre = as.numeric(mz_pre) + offset,
        mz = as.numeric(mz) + offset,
        adduct_label = adduct_name
      )
  })

  df_theory_expanded <- df_theory_expanded %>%
    bind_rows(
      df_theory_expanded %>%
        filter(grepl("LF|CF", compound_name)) %>%
        mutate(
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
    filter(ion_type %in% c("b", "y")) %>%
    group_by(compound_name) %>%
    summarise(Total_Theoretical_by = n_distinct(paste0(ion_type, ion_num)), .groups = 'drop')

  cat("Done\n")
  return(list(spec_library = spec_library, theory_ion_counts = theory_ion_counts, df_theory = df_theory))
}


#' Core Batch Processing Workflow
#'
#' The primary engine of the package. It iterates through mzML files, performs
#' spectrum cleaning, identifies common peaks against the theoretical library,
#' filters by MS1/MS2 tolerance, and generates initial identification reports.
#'
#' @param lib_data List. The output from the \code{prepLib} function.
#' @param CONFIG List. Global configuration containing mass tolerances and scoring weights.
#' @export
doMain <- function(lib_data, CONFIG) {
  spec_library <- lib_data$spec_library
  theory_ion_counts <- lib_data$theory_ion_counts
  df_theory <- lib_data$df_theory

  mzML_files <- list.files(pattern = "\\.mzML$", full.names = FALSE)
  n_files <- length(mzML_files)
  tol_val <- CONFIG$ms2_tolerance_ppm * 1e-6

  cat(sprintf("\n--- Starting batch processing workflow. Found %d files in total. ---\n", n_files))

  for (f_idx in seq_along(mzML_files)) {
    file_name <- mzML_files[f_idx]
    file_prefix <- str_remove(file_name, "\\.mzML$")
    cat(sprintf("\n[%d/%d] processing: %s\n", f_idx, n_files, file_name))

    raw_exp <- readMSData(file_name, msLevel = 2, mode = "onDisk")
    exp_specs <- spectra(raw_exp)
    names(exp_specs) <- paste0("Scan_", fData(raw_exp)$acquisitionNum)

    exp_specs_clean <- lapply(exp_specs, function(s) {
      s <- pickPeaks(s) %>% clean(all = TRUE)
      max_int <- max(intensity(s), na.rm = TRUE)
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
        common_peaks_matrix[i, j] <- compareSpectra(s_i, spec_library[[j]],
                                                    fun = "common", relative = TRUE, tolerance = tol_val)
      }
      setTxtProgressBar(pb, i)
    }
    close(pb)

    matched_results <- as.data.frame(common_peaks_matrix) %>%
      rownames_to_column("Experimental_Scan") %>%
      pivot_longer(cols = -Experimental_Scan, names_to = "Theoretical_Lipopeptide", values_to = "Common_Peaks") %>%
      filter(Common_Peaks > 7)

    if(nrow(matched_results) == 0) {
      cat("    [Skip] No matches found for this file\n")
      next
    }

    exp_metadata <- fData(raw_exp) %>% rownames_to_column("Scan_ID") %>%
      mutate(Experimental_Scan = paste0("Scan_", acquisitionNum),
             precursorMZ = round(as.numeric(precursorMZ), 4)) %>%
      select(Experimental_Scan, precursorMZ)

    matched_results <- matched_results %>% left_join(exp_metadata, by = "Experimental_Scan")

    filtered_mz_results <- matched_results %>%
      mutate(base_theoretical_mz = as.numeric(sub("^([0-9.]+)_.*", "\\1", Theoretical_Lipopeptide))) %>%
      rowwise() %>%
      mutate(
        errors_ppm = list(abs(precursorMZ - (base_theoretical_mz + CONFIG$adduct_offset)) /
                            (base_theoretical_mz + CONFIG$adduct_offset) * 1e6),
        matched_adduct = list(names(CONFIG$adduct_offset)[unlist(errors_ppm) < CONFIG$ms1_match_ppm])
      ) %>%
      filter(length(matched_adduct) > 0) %>%
      mutate(matched_adduct = paste(unlist(matched_adduct), collapse = "/")) %>%
      ungroup() %>% filter(Common_Peaks > 7)

    final_identified_lps <- filtered_mz_results %>%
      group_by(Theoretical_Lipopeptide) %>% slice_max(Common_Peaks, n = 1, with_ties = FALSE) %>% ungroup()

    final_identified_lps_output <- final_identified_lps %>%
      mutate(parts = str_split(Theoretical_Lipopeptide, "_"),
             LP_Type = sapply(parts, function(x) if (length(x) >= 4 && (x[4] == "CF" || x[4] == "LF")) x[4] else x[3]),
             AA_Sequence = sapply(parts, function(x) x[2]),
             FA_Details = sapply(parts, function(x) {
               if (length(x) >= 4 && (x[4] == "CF" || x[4] == "LF")) {
                 valid_idx <- c(3, 5, 6); return(paste(x[valid_idx[valid_idx <= length(x)]], collapse = "_"))
               } else {
                 if (length(x) >= 5 && (x[3] == "CF" || x[3] == "LF")) paste(x[4], x[5], sep = "_") else x[4]
               }
             })) %>%
      left_join(df_theory %>% select(compound_name, db) %>% distinct(), by = c("Theoretical_Lipopeptide" = "compound_name")) %>%
      select(LP_Type, FA_Details, AA_Sequence, everything(), -parts)

    detailed_matches_list <- list()
    for (i in 1:nrow(final_identified_lps)) {
      row_data <- final_identified_lps[i, ]; exp_id <- row_data$Experimental_Scan; theo_id <- row_data$Theoretical_Lipopeptide
      s_exp <- exp_specs_clean[[exp_id]]; s_theo <- spec_library[[theo_id]]
      mz_exp <- mz(s_exp); mz_theo <- mz(s_theo)
      matches <- lapply(mz_exp, function(m) {
        diffs <- abs(mz_theo - m) / m
        best <- which.min(diffs)
        if (length(best) > 0 && diffs[best] <= CONFIG$ms2_tolerance_ppm * 1e-6) {
          return(data.frame(Matched_mz_exp = m, Matched_mz_theo = mz_theo[best],
                            Intensity_exp = intensity(s_exp)[which(mz_exp == m)[1]], ppm_error = diffs[best] * 1e6))
        }
        return(NULL)
      })
      curr_df <- do.call(rbind, matches)
      if (!is.null(curr_df)) {
        curr_df$Experimental_Scan <- exp_id; curr_df$Theoretical_Lipopeptide <- theo_id
        detailed_matches_list[[i]] <- curr_df
      }
    }
    matched_ions_detail <- do.call(rbind, detailed_matches_list)
    lookup_table <- df_theory %>%
      mutate(ion_label = paste0(ion_type, ion_num), add_label = adduct_label) %>%
      select(Theoretical_Lipopeptide = compound_name, Matched_mz_theo = mz, ion_label, add_label) %>%
      distinct(Theoretical_Lipopeptide, Matched_mz_theo, .keep_all = TRUE)

    feature_lookup <- final_identified_lps_output %>%
      select(Experimental_Scan, Theoretical_Lipopeptide, LP_Type, FA_Details, AA_Sequence, matched_adduct, db)

    final_output_annotated <- matched_ions_detail %>%
      left_join(exp_metadata, by = "Experimental_Scan") %>%
      left_join(lookup_table, by = c("Theoretical_Lipopeptide", "Matched_mz_theo")) %>%
      left_join(feature_lookup, by = c("Experimental_Scan", "Theoretical_Lipopeptide")) %>%
      rename(Ion_Annotation = ion_label) %>%
      select(Experimental_Scan, LP_Type, precursorMZ, Theoretical_Lipopeptide, matched_adduct,
             Matched_mz_exp, Matched_mz_theo, ppm_error, Intensity_exp, add_label,Ion_Annotation,db)

    final_output_annotated_out <- final_output_annotated %>% filter(db == "target") %>%
      select(Experimental_Scan, LP_Type, precursorMZ, Theoretical_Lipopeptide, matched_adduct,
             Matched_mz_exp, Matched_mz_theo, ppm_error, Intensity_exp, add_label,Ion_Annotation)

    write.csv(final_output_annotated_out, paste0(file_prefix, "_MS2_detail.csv"), row.names = FALSE)

    scoring_details <- final_output_annotated %>%
      left_join(theory_ion_counts, by = c("Theoretical_Lipopeptide" = "compound_name")) %>%
      group_by(Experimental_Scan, Theoretical_Lipopeptide) %>%
      summarise(
        Score_Precision = pmax(0, 1 - (median(abs(ppm_error), na.rm = TRUE) / CONFIG$ms2_tolerance_ppm)),
        Score_Continuity = {
          b_nums <- sort(as.numeric(str_extract(Ion_Annotation[grepl("^b", Ion_Annotation)], "\\d+")))
          y_nums <- sort(as.numeric(str_extract(Ion_Annotation[grepl("^y", Ion_Annotation)], "\\d+")))
          count_cont <- function(nums) if(length(nums) < 2) 0 else sum(diff(nums) == 1)
          total_cont_pairs <- count_cont(b_nums) + count_cont(y_nums)
          theoretical_max_pairs <- pmax(1, first(Total_Theoretical_by) - 2)
          pmin(1.0, total_cont_pairs / theoretical_max_pairs)
        },
        Score_Intensity = {
          all_intensities <- intensity(exp_specs_clean[[first(Experimental_Scan)]])
          base_peak_int <- max(all_intensities, na.rm = TRUE)
          significant_mask <- all_intensities > (base_peak_int * CONFIG$int_base_peak_ratio)
          significant_peaks_sum <- sum(all_intensities[significant_mask], na.rm = TRUE)
          matched_significant_int <- sum(Intensity_exp[Intensity_exp > (base_peak_int * CONFIG$int_base_peak_ratio)], na.rm = TRUE)
          tic_ratio <- if(significant_peaks_sum > 0) pmin(1.0, matched_significant_int / significant_peaks_sum) else 0
          avg_rel_int <- mean(Intensity_exp / base_peak_int, na.rm = TRUE)
          intensity_penalty <- if_else(avg_rel_int > CONFIG$int_avg_rel_threshold, 1.0, (avg_rel_int / CONFIG$int_penalty_coef))
          pmin(1.0, tic_ratio * intensity_penalty)
        },
        Matched_Unique_Ions = n_distinct(Ion_Annotation),
        Total_Theory_Ions = first(Total_Theoretical_by),
        Score_Coverage = if_else(Matched_Unique_Ions >= CONFIG$cov_min_unique_ions,
                                 Matched_Unique_Ions / Total_Theory_Ions,
                                 (Matched_Unique_Ions / Total_Theory_Ions) * CONFIG$cov_penalty_factor),
        Score_End_Match = {
          current_type <- first(LP_Type)
          target_y <- if(current_type %in% c("LF", "CF")) "y10" else "y7"
          has_b1 = any(Ion_Annotation == "b1", na.rm = TRUE); has_target_y = any(Ion_Annotation == target_y, na.rm = TRUE)
          case_when(has_b1 & has_target_y ~ 1.0, has_b1 | has_target_y ~ 0.3, TRUE ~ 0.0)
        },
        Score_Adduct = {
          curr_adduct <- as.character(first(matched_adduct))
          val <- CONFIG$adduct_scores[curr_adduct]
          if(is.na(val)) CONFIG$adduct_scores["Default"] else val
        }, .groups = 'drop'
      )

    final_report_v2 <- final_identified_lps_output %>%
      left_join(scoring_details, by = c("Experimental_Scan", "Theoretical_Lipopeptide")) %>%
      mutate(LipopepID_Score = (Score_Coverage * CONFIG$weights$coverage) + (Score_Continuity * CONFIG$weights$continuity) +
               (Score_Precision * CONFIG$weights$precision) + (Score_Intensity * CONFIG$weights$intensity) +
               (Score_End_Match * CONFIG$weights$end_match) + (Score_Adduct * CONFIG$weights$adduct)) %>%
      mutate(LipopepID_Score = pmin(1.0, round(LipopepID_Score, 4))) %>%
      left_join(fData(raw_exp) %>% rownames_to_column("Scan_ID") %>%
                  mutate(Experimental_Scan = paste0("Scan_", acquisitionNum)) %>%
                  select(Experimental_Scan, precursorIntensity), by = "Experimental_Scan") %>%
      rename(MS1_Intensity = precursorIntensity)

    fdr_calc <- final_report_v2 %>% arrange(desc(LipopepID_Score)) %>%
      mutate(cum_target = cumsum(db == "target"), cum_decoy = cumsum(db == "decoy"),
             current_FDR = (2 * cum_decoy) / pmax(1, (cum_target + cum_decoy)))

    cutoff_1 <- fdr_calc %>% filter(current_FDR <= 0.01) %>% slice_tail(n = 1) %>% pull(LipopepID_Score)
    cutoff_1 <- if(length(cutoff_1) == 0) max(fdr_calc$LipopepID_Score) else cutoff_1

    final_report_v2_export <- final_report_v2 %>% filter(db == "target") %>%
      mutate(across(where(is.list), ~ map_chr(.x, ~ paste(as.character(.x), collapse = "; "))))

    final_report_v2_export <- final_report_v2_export %>%
      mutate(FA_processing = if_else(LP_Type %in% c("LF", "CF"), FA_Details, NA_character_)) %>%
      mutate(parts = str_split(FA_processing, "_"),
             FA_Details = if_else(LP_Type %in% c("LF", "CF"), map_chr(parts, ~ .x[2]), FA_Details),
             Additional_Info = if_else(LP_Type %in% c("LF", "CF"),
                                       map_chr(parts, function(x) if(length(x) >= 3) paste(x[3:length(x)], collapse = "_") else NA_character_), NA_character_)) %>%
      select(-FA_processing, -parts)

    write.csv(final_report_v2_export, paste0(file_prefix, "_MS1_Results.csv"), row.names = FALSE)

    final_report_v2_confidence <- final_report_v2[final_report_v2$db == "target" & final_report_v2$LipopepID_Score >= cutoff_1, ]
    if(nrow(final_report_v2_confidence) > 0) {
      conf_export <- final_report_v2_confidence %>% mutate(across(where(is.list), ~ map_chr(.x, ~ paste(as.character(.x), collapse = "; "))))
      write.csv(conf_export, paste0(file_prefix, "_MS1_results_1%FDR.csv"), row.names = FALSE)
    }

    intensity_calculation_support <- map_df(exp_specs_clean, function(s) {
      base_int <- max(intensity(s), na.rm = TRUE)
      sig_tic <- if(base_int > 0) sum(intensity(s)[intensity(s) > (base_int * CONFIG$int_base_peak_ratio)], na.rm = TRUE) else 0
      data.frame(Experimental_Scan = paste0("Scan_", acquisitionNum(s)), Base_Peak_Intensity = base_int, TIC_Significant_Ref = sig_tic, stringsAsFactors = FALSE)
    })

    rdata_file <- paste0(file_prefix, "_Scoring_Workspace.RData")
    save(final_output_annotated, final_identified_lps_output, theory_ion_counts, intensity_calculation_support, CONFIG, raw_exp, file_prefix, file = rdata_file)
  }
  cat("\n--- Done ---\n")
}

#' Statistical Summary and Result Aggregation
#'
#' Consolidates results across multiple processed samples. It performs
#' data cleaning (e.g., removing contradictory adducts), calculates abundance
#' statistics, and organizes 1% FDR filtered results into dedicated directories.
#'
#' @export
doStat <- function() {
  cat("\n--- Starting result organization and statistical aggregation workflow ---\n")


  process_ms1_file <- function(file_path) {
    df <- read_csv(file_path, show_col_types = FALSE)
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df_cleaned <- df %>%
      filter(!(LP_Type == "CS" & matched_adduct == "H2O_add")) %>%
      filter(!(LP_Type == "LS" & matched_adduct == "H2O_loss")) %>%
      filter(!(LP_Type == "CF" & matched_adduct == "H2O_add")) %>%
      filter(!(LP_Type == "LF" & matched_adduct == "H2O_loss")) %>%
      mutate(FA_Details = case_when(
        LP_Type %in% c("LF", "CF") & str_detect(FA_Details, "_") ~ str_extract(FA_Details, "(?<=_)\\d+$"),
        TRUE ~ as.character(FA_Details)
      ))
    write_csv(df_cleaned, file_path)
  }

  result_files <- list.files(pattern = "_MS1_.*\\.csv$", full.names = FALSE)
  invisible(lapply(result_files, process_ms1_file))


  res_files <- list.files(pattern = "_MS1_Results.csv$", full.names = FALSE)
  for (res_file in res_files) {
    file_prefix <- str_remove(res_file, "_MS1_Results.csv$")
    df <- read.csv(res_file) %>% filter(db == "target")
    if(nrow(df) == 0) next

    summary_lp_type <- df %>% group_by(LP_Type) %>% summarise(Count = n(), Total_Abundance = sum(MS1_Intensity, na.rm = TRUE), Abundance_Pct = sum(MS1_Intensity, na.rm = TRUE) / sum(df$MS1_Intensity, na.rm = TRUE), .groups = 'drop')
    summary_fa_details <- df %>% group_by(LP_Type, FA_Details, AA_Sequence) %>% summarise(Observed_MZ = mean(precursorMZ, na.rm = TRUE), Scan_Count = n(), Max_Intensity = max(MS1_Intensity, na.rm = TRUE), Total_Abundance = sum(MS1_Intensity, na.rm = TRUE), .groups = 'drop') %>% group_by(LP_Type) %>% mutate(Abundance_in_Type = Total_Abundance / sum(Total_Abundance)) %>% ungroup()

    write.csv(summary_lp_type, paste0(file_prefix, "_Stat_LP_Type.csv"), row.names = FALSE)
    write.csv(summary_fa_details, paste0(file_prefix, "_Stat_FA_Details.csv"), row.names = FALSE)
  }


  lp_stat_files <- list.files(pattern = "_Stat_LP_Type.csv$", full.names = FALSE)
  fa_stat_files <- list.files(pattern = "_Stat_FA_Details.csv$", full.names = FALSE)
  all_samples_lp_summary <- map_df(lp_stat_files, function(f) { sample_name <- str_remove(f, "_Stat_LP_Type.csv$"); read.csv(f) %>% mutate(Sample_Source = sample_name) %>% relocate(Sample_Source) })
  all_samples_fa_summary <- map_df(fa_stat_files, function(f) { sample_name <- str_remove(f, "_Stat_FA_Details.csv$"); read.csv(f) %>% mutate(Sample_Source = sample_name) %>% relocate(Sample_Source) })
  write.csv(all_samples_lp_summary, "All_Samples_LP_Type.csv", row.names = FALSE)
  write.csv(all_samples_fa_summary, "All_Samples_FA_Details.csv", row.names = FALSE)


  fdr_res_files <- list.files(pattern = "FDR.*\\.csv$", full.names = FALSE)
  output_dir <- "result_FDR"
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  for (res_file in fdr_res_files) {
    file_prefix <- str_extract(res_file, "^.*(?=_MS1)")
    df <- read.csv(res_file) %>% filter(db == "target")
    if(nrow(df) == 0) next
    summary_lp_type <- df %>% group_by(LP_Type) %>% summarise(Count = n(), Total_Abundance = sum(MS1_Intensity, na.rm = TRUE), Abundance_Pct = sum(MS1_Intensity, na.rm = TRUE) / sum(df$MS1_Intensity, na.rm = TRUE), .groups = 'drop')
    summary_fa_details <- df %>% group_by(LP_Type, FA_Details, AA_Sequence) %>% summarise(Observed_MZ = mean(precursorMZ, na.rm = TRUE), Scan_Count = n(), Max_Intensity = max(MS1_Intensity, na.rm = TRUE), Total_Abundance = sum(MS1_Intensity, na.rm = TRUE), .groups = 'drop') %>% group_by(LP_Type) %>% mutate(Abundance_in_Type = Total_Abundance / sum(Total_Abundance)) %>% ungroup()
    write.csv(summary_lp_type, file = file.path(output_dir, paste0(file_prefix, "_Stat_LP_Type_1%FDR.csv")), row.names = FALSE)
    write.csv(summary_fa_details, file = file.path(output_dir, paste0(file_prefix, "_Stat_FA_Details_1%FDR.csv")), row.names = FALSE)
  }

  lp_fdr_files <- list.files(path = output_dir, pattern = "_Stat_LP_Type_1%FDR.csv$", full.names = TRUE)
  fa_fdr_files <- list.files(path = output_dir, pattern = "_Stat_FA_Details_1%FDR.csv$", full.names = TRUE)
  all_fdr_lp <- map_df(lp_fdr_files, function(f) { sample_name <- str_remove(basename(f), "_Stat_LP_Type_1%FDR.csv$"); read.csv(f) %>% mutate(Sample_Source = sample_name) %>% relocate(Sample_Source) })
  all_fdr_fa <- map_df(fa_fdr_files, function(f) { sample_name <- str_remove(basename(f), "_Stat_FA_Details_1%FDR.csv$"); read_csv(f, col_types = cols(.default = "c")) %>% mutate(Sample_Source = sample_name) %>% relocate(Sample_Source) })
  write.csv(all_fdr_fa, file.path(output_dir, "All_Samples_FA_Details_FDR.csv"), row.names = FALSE)
  write.csv(all_fdr_lp, file.path(output_dir, "All_Samples_LP_Type_FDR.csv"), row.names = FALSE)

  target_dir <- "result_FDR"
  if (!dir.exists(target_dir)) dir.create(target_dir)
  fdr_items <- list.files(pattern = "FDR", full.names = TRUE)
  fdr_files <- fdr_items[!dir.exists(fdr_items)]
  if (length(fdr_files) > 0) {
    file.copy(fdr_files, file.path(target_dir, basename(fdr_files)))
    file.remove(fdr_files)
  }
  cat("--- Statistical summary complete. ---\n")
}


#' Re-scoring and FDR Calculation
#'
#' Performs a multi-dimensional scoring analysis on previously identified features
#' using saved intermediate data. It calculates the Final LipopepID Score and
#' determines the 1% False Discovery Rate (FDR) threshold using a Target-Decoy approach.
#'
#' @param rdata_path Character. Path to the intermediate workspace file (e.g., "*_Scoring_Workspace.RData").
#' @param CONFIG List. Global configuration list. Weights provided here will override those used in the original \code{doMain} run.
#' @export
doScore <- function(rdata_path, CONFIG) {

  # attention：here will load final_output_annotated, final_identified_lps_output,
  # theory_ion_counts, intensity_calculation_support, file_prefix, raw_exp
  load(rdata_path)
  cat(sprintf("\n--- Rescoring: %s ---\n", file_prefix))


  scoring_details <- final_output_annotated %>%
    left_join(theory_ion_counts, by = c("Theoretical_Lipopeptide" = "compound_name")) %>%
    group_by(Experimental_Scan, Theoretical_Lipopeptide) %>%
    summarise(

      Score_Precision = pmax(0, 1 - (median(abs(ppm_error), na.rm = TRUE) / CONFIG$ms2_tolerance_ppm)),


      Score_Continuity = {
        b_nums <- sort(as.numeric(str_extract(Ion_Annotation[grepl("^b", Ion_Annotation)], "\\d+")))
        y_nums <- sort(as.numeric(str_extract(Ion_Annotation[grepl("^y", Ion_Annotation)], "\\d+")))
        count_cont <- function(nums) if(length(nums) < 2) 0 else sum(diff(nums) == 1)
        total_cont_pairs <- count_cont(b_nums) + count_cont(y_nums)
        theoretical_max_pairs <- pmax(1, first(Total_Theoretical_by) - 2)
        pmin(1.0, total_cont_pairs / theoretical_max_pairs)
      },


      Score_Intensity = {
        stats <- intensity_calculation_support[intensity_calculation_support$Experimental_Scan == first(Experimental_Scan), ]
        base_peak_int <- stats$Base_Peak_Intensity
        significant_peaks_sum <- stats$TIC_Significant_Ref

        matched_significant_int <- sum(Intensity_exp[Intensity_exp > (base_peak_int * CONFIG$int_base_peak_ratio)], na.rm = TRUE)
        tic_ratio <- if(significant_peaks_sum > 0) pmin(1.0, matched_significant_int / significant_peaks_sum) else 0

        avg_rel_int <- mean(Intensity_exp / base_peak_int, na.rm = TRUE)
        intensity_penalty <- if_else(avg_rel_int > CONFIG$int_avg_rel_threshold, 1.0, (avg_rel_int / CONFIG$int_penalty_coef))
        pmin(1.0, tic_ratio * intensity_penalty)
      },


      Matched_Unique_Ions = n_distinct(Ion_Annotation),
      Total_Theory_Ions = first(Total_Theoretical_by),
      Score_Coverage = if_else(Matched_Unique_Ions >= CONFIG$cov_min_unique_ions,
                               Matched_Unique_Ions / Total_Theory_Ions,
                               (Matched_Unique_Ions / Total_Theory_Ions) * CONFIG$cov_penalty_factor),


      Score_End_Match = {
        current_type <- first(LP_Type)
        target_y <- if(current_type %in% c("LF", "CF")) "y10" else "y7"
        has_b1 = any(Ion_Annotation == "b1", na.rm = TRUE)
        has_target_y = any(Ion_Annotation == target_y, na.rm = TRUE)
        case_when(has_b1 & has_target_y ~ 1.0, has_b1 | has_target_y ~ 0.3, TRUE ~ 0.0)
      },


      Score_Adduct = {
        curr_adduct <- as.character(first(matched_adduct))
        val <- CONFIG$adduct_scores[curr_adduct]
        if(is.na(val)) CONFIG$adduct_scores["Default"] else val
      }, .groups = 'drop'
    )


  final_report_v2 <- final_identified_lps_output %>%
    left_join(scoring_details, by = c("Experimental_Scan", "Theoretical_Lipopeptide")) %>%
    mutate(LipopepID_Score =
             (Score_Coverage   * CONFIG$weights$coverage) +
             (Score_Continuity * CONFIG$weights$continuity) +
             (Score_Precision  * CONFIG$weights$precision) +
             (Score_Intensity  * CONFIG$weights$intensity) +
             (Score_End_Match  * CONFIG$weights$end_match) +
             (Score_Adduct     * CONFIG$weights$adduct)) %>%
    mutate(LipopepID_Score = pmin(1.0, round(LipopepID_Score, 4))) %>%
    left_join(fData(raw_exp) %>% rownames_to_column("Scan_ID") %>%
                mutate(Experimental_Scan = paste0("Scan_", acquisitionNum)) %>%
                select(Experimental_Scan, precursorIntensity), by = "Experimental_Scan") %>%
    rename(MS1_Intensity = precursorIntensity)


  fdr_calc <- final_report_v2 %>% arrange(desc(LipopepID_Score)) %>%
    mutate(cum_target = cumsum(db == "target"), cum_decoy = cumsum(db == "decoy"),
           current_FDR = (2 * cum_decoy) / pmax(1, (cum_target + cum_decoy)))

  cutoff_1 <- fdr_calc %>% filter(current_FDR <= 0.01) %>% slice_tail(n = 1) %>% pull(LipopepID_Score)
  cutoff_1 <- if(length(cutoff_1) == 0) max(fdr_calc$LipopepID_Score) else cutoff_1


  final_report_v2_export <- final_report_v2 %>%
    mutate(across(where(is.list), ~ map_chr(.x, ~ paste(as.character(.x), collapse = "; "))))

  final_report_v2_export <- final_report_v2_export %>%
    mutate(FA_processing = if_else(LP_Type %in% c("LF", "CF"), FA_Details, NA_character_)) %>%
    mutate(
      parts = str_split(FA_processing, "_"),
      FA_Details = if_else(LP_Type %in% c("LF", "CF"), map_chr(parts, ~ .x[2]), FA_Details),
      Additional_Info = if_else(LP_Type %in% c("LF", "CF"),
                                map_chr(parts, function(x) {
                                  if(length(x) >= 3) paste(x[3:length(x)], collapse = "_")
                                  else NA_character_
                                }), NA_character_)
    ) %>%
    select(-FA_processing, -parts)

  write.csv(final_report_v2_export, paste0(file_prefix, "_MS1_Results.csv"), row.names = FALSE)

  final_report_v2_confidence <- final_report_v2[final_report_v2$db == "target" & final_report_v2$LipopepID_Score >= cutoff_1, ]
  if(nrow(final_report_v2_confidence) > 0) {
    conf_export <- final_report_v2_confidence %>%
      mutate(across(where(is.list), ~ map_chr(.x, ~ paste(as.character(.x), collapse = "; "))))
    write.csv(conf_export, paste0(file_prefix, "_MS1_results_1%FDR.csv"), row.names = FALSE)
  }

  cat(sprintf("   [OK] Re-scoring complete. Identifications at 1%% FDR: %d\n", nrow(final_report_v2_confidence)))
}
