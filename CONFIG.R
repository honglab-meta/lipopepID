# =========================================================
# Global Configuration Center (Adjust parameters here)
# =========================================================
CONFIG <- list(
  # 1. Mass Tolerance Settings 
  ms2_tolerance_ppm = 20,          # Used for compareSpectra and Score_Precision
  ms1_match_ppm = 10,              # Precursor mass match tolerance (PPM)
  
  # 2. Adduct Offset Settings
  # format: M(adduct) - M(H+) (Relative to protonated form)
  adduct_offset = c(
    "H"        = 0,
    "Na"       = 22.9898 - 1.0078,
    "NH4"      = 18.0344 - 1.0078,
    "H2O_loss" = -18.0106,
    "H2O_add"  = 18.0106),
  
  # 3. Intensity Scoring Parameters
  int_base_peak_ratio = 0.05,      # Threshold to define "significant peaks"
  int_avg_rel_threshold = 0.15,    # Average relative intensity threshold
  int_penalty_coef = 0.3,          # Denominator for the penalty coefficient
  
  # 4. Coverage Scoring Parameters
  cov_min_unique_ions = 3,         # Minimum threshold for unique ion count
  cov_penalty_factor = 0.5,        # Penalty factor when below the threshold
  
  # 5. Adduct Preference Scores (Corresponds to Score_Adduct)
  adduct_scores = c(
    "H" = 1.0, "NH4" = 0.7, "Na" = 0.8, "H2O_add" = 0.9, "H2O_loss" = 0.9, "Default" = 0.5),
  
  # 6. Total Score Weights (Sum should be 1.0)
  weights = list(
    coverage   = 0.25,
    continuity = 0.25,
    precision  = 0.25,
    intensity  = 0.15,
    end_match  = 0.05,
    adduct     = 0.05
  )
)