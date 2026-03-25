#!/usr/bin/env Rscript

# =============================================================================
# NDHS 2024 Caesarean Section Risk Pipeline -- HARD-FIX FINAL SCRIPT
# -----------------------------------------------------------------------------
# Stable BR-only production pipeline for NDHS 2024 Birth Recode (NGBR8AFL.dta)
#
# Architecture
#   1) Survey-weighted descriptive epidemiology for manuscript tables/figures
#   2) BR-only phenotype engineering
#   3) Weighted ridge logistic regression (glmnet) for prediction
#   4) Optional reduced weighted logistic refit for interpretable OR display
#      with automatic fallback to ridge-coefficient display
#   5) PSU-within-stratum development/validation split
#   6) Five figures (PNG, 650 dpi) + five tables (CSV + HTML)
# =============================================================================

options(
  stringsAsFactors = FALSE,
  scipen = 999,
  width = 170,
  survey.lonely.psu = "adjust"
)

set.seed(20260309)

# -----------------------------------------------------------------------------
# Package management
# -----------------------------------------------------------------------------
required_pkgs <- c(
  "haven", "dplyr", "tidyr", "stringr", "forcats", "purrr", "readr",
  "tibble", "ggplot2", "survey", "gt", "htmltools", "scales",
  "Matrix", "glmnet"
)

install_if_missing <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss) > 0) {
    install.packages(miss, repos = "https://cloud.r-project.org")
  }
  invisible(lapply(pkgs, library, character.only = TRUE))
}
suppressPackageStartupMessages(install_if_missing(required_pkgs))

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x)) || identical(x, "")) y else x
}

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

find_existing_file <- function(candidates, label) {
  candidates <- unique(candidates[!is.na(candidates) & nzchar(candidates)])
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0) {
    stop(sprintf("Could not locate %s. Checked:\n%s", label, paste(candidates, collapse = "\n")), call. = FALSE)
  }
  hit[[1]]
}

br_path <- find_existing_file(
  c(
    args[1] %||% "",
    Sys.getenv("NDHS_BR_PATH", unset = ""),
    "NGBR8AFL.dta",
    file.path("data", "NGBR8AFL.dta"),
    file.path(getwd(), "NGBR8AFL.dta")
  ),
  "NDHS BR dataset (NGBR8AFL.dta)"
)

analysis_root <- file.path(getwd(), "ndhs2024_cs_risk_final_outputs")
dir.create(analysis_root, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(analysis_root, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(analysis_root, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(analysis_root, "models"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(analysis_root, "logs"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(analysis_root, "supplementary"), recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Formatting helpers
# -----------------------------------------------------------------------------
fmt_n <- function(x, digits = 0) {
  ifelse(is.na(x), NA_character_, scales::number(x, accuracy = 10^-digits, big.mark = ",", trim = TRUE))
}

fmt_pct <- function(x, digits = 1) {
  ifelse(is.na(x), NA_character_, paste0(scales::number(100 * x, accuracy = 10^-digits, trim = TRUE), "%"))
}

fmt_num <- function(x, digits = 2) {
  ifelse(is.na(x), NA_character_, scales::number(x, accuracy = 10^-digits, trim = TRUE))
}

fmt_p <- function(p) {
  ifelse(is.na(p), NA_character_, ifelse(p < 1e-4, "<0.0001", sprintf("%.4f", p)))
}

clip_prob <- function(p) {
  p <- as.numeric(p)
  p[!is.finite(p)] <- NA_real_
  pmin(pmax(p, 1e-6), 1 - 1e-6)
}

safe_weighted_mean <- function(x, w) {
  x <- as.numeric(x)
  w <- as.numeric(w)
  ok <- !is.na(x) & !is.na(w) & is.finite(x) & is.finite(w)
  if (!any(ok)) return(NA_real_)
  stats::weighted.mean(x[ok], w[ok])
}

safe_weighted_sd <- function(x, w) {
  x <- as.numeric(x)
  w <- as.numeric(w)
  ok <- !is.na(x) & !is.na(w) & is.finite(x) & is.finite(w)
  if (!any(ok)) return(NA_real_)
  x <- x[ok]
  w <- w[ok]
  mu <- stats::weighted.mean(x, w)
  sqrt(sum(w * (x - mu)^2) / sum(w))
}

weighted_median <- function(x, w) {
  x <- as.numeric(x)
  w <- as.numeric(w)
  ok <- !is.na(x) & !is.na(w) & is.finite(x) & is.finite(w)
  if (!any(ok)) return(NA_real_)
  x <- x[ok]
  w <- w[ok]
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  cw <- cumsum(w) / sum(w)
  x[which(cw >= 0.5)[1]]
}

safe_tail_n <- function(df, n = 25) {
  if (is.null(df) || nrow(df) == 0) return(df)
  utils::tail(df, n = min(n, nrow(df)))
}

# -----------------------------------------------------------------------------
# Variable helpers
# -----------------------------------------------------------------------------
safe_as_factor <- function(x) {
  if (inherits(x, "haven_labelled") || inherits(x, "labelled")) {
    haven::as_factor(x, levels = "labels")
  } else if (is.factor(x)) {
    x
  } else {
    factor(x)
  }
}

yes_no_unknown_factor <- function(x, yes = 1, no = 0, missing_codes = c(8, 9, 98, 99, 998, 999)) {
  x_num <- suppressWarnings(as.numeric(x))
  out <- dplyr::case_when(
    is.na(x_num) ~ "Unknown",
    x_num %in% missing_codes ~ "Unknown",
    x_num %in% yes ~ "Yes",
    x_num %in% no ~ "No",
    TRUE ~ "Unknown"
  )
  factor(out, levels = c("No", "Yes", "Unknown"))
}

collapse_with_train <- function(train_vec, test_vec, min_n = 50, other_label = "Other/unknown") {
  tr <- as.character(train_vec)
  te <- as.character(test_vec)
  
  tr[is.na(tr) | !nzchar(tr)] <- "Unknown"
  te[is.na(te) | !nzchar(te)] <- "Unknown"
  
  tab <- table(tr)
  rare <- names(tab)[tab < min_n]
  tr[tr %in% rare] <- other_label
  
  allowed <- sort(unique(tr))
  te[!(te %in% allowed)] <- other_label
  
  lvls <- sort(unique(c(tr, te)))
  list(
    train = factor(tr, levels = lvls),
    test  = factor(te, levels = lvls)
  )
}

has_variation <- function(x) {
  x2 <- x[!is.na(x)]
  if (length(x2) == 0) return(FALSE)
  if (is.factor(x2) || is.character(x2)) {
    length(unique(as.character(x2))) > 1
  } else {
    length(unique(x2)) > 1
  }
}

make_psu_split <- function(df, strata_var = "v022", psu_var = "v021", prop_train = 0.80, seed = 20260309) {
  set.seed(seed)
  
  tmp <- df |>
    dplyr::mutate(.strata = .data[[strata_var]], .psu = .data[[psu_var]])
  
  train_psus <- tmp |>
    dplyr::distinct(.strata, .psu) |>
    dplyr::group_by(.strata) |>
    dplyr::summarise(
      .psu = list({
        ids <- unique(.psu)
        n_train <- max(1, floor(length(ids) * prop_train))
        sample(ids, size = n_train, replace = FALSE)
      }),
      .groups = "drop"
    ) |>
    tidyr::unnest(.psu) |>
    dplyr::mutate(split = "Development")
  
  tmp |>
    dplyr::left_join(train_psus, by = c(".strata", ".psu")) |>
    dplyr::mutate(split = dplyr::if_else(is.na(split), "Validation", split)) |>
    dplyr::select(-.strata, -.psu)
}

impute_center_continuous <- function(train, test, raw_var, imp_var, miss_var) {
  med <- weighted_median(train[[raw_var]], train$wt)
  if (is.na(med)) med <- stats::median(train[[raw_var]], na.rm = TRUE)
  if (is.na(med) || !is.finite(med)) med <- 0
  
  train[[miss_var]] <- as.integer(is.na(train[[raw_var]]))
  test[[miss_var]]  <- as.integer(is.na(test[[raw_var]]))
  
  train[[imp_var]] <- ifelse(is.na(train[[raw_var]]), med, train[[raw_var]])
  test[[imp_var]]  <- ifelse(is.na(test[[raw_var]]), med, test[[raw_var]])
  
  train[[imp_var]] <- train[[imp_var]] - med
  test[[imp_var]]  <- test[[imp_var]] - med
  
  list(train = train, test = test, median = med)
}

# -----------------------------------------------------------------------------
# Plot theme / saving
# -----------------------------------------------------------------------------
ndhs_theme <- function() {
  ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 15, colour = "#102A43"),
      plot.subtitle = ggplot2::element_text(size = 10.5, colour = "#334E68"),
      plot.caption = ggplot2::element_text(size = 9, colour = "#486581", hjust = 0),
      axis.title = ggplot2::element_text(face = "bold", colour = "#243B53"),
      axis.text = ggplot2::element_text(colour = "#243B53"),
      strip.text = ggplot2::element_text(face = "bold", colour = "#102A43"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(colour = "#D9E2EC", linewidth = 0.35),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(face = "bold"),
      plot.margin = ggplot2::margin(16, 26, 16, 16)
    )
}

save_png <- function(plot_obj, filename, width = 10, height = 7, dpi = 650) {
  ggplot2::ggsave(
    filename = file.path(analysis_root, "figures", filename),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = dpi,
    units = "in",
    bg = "white",
    limitsize = FALSE
  )
}

save_gt <- function(df, file_stub, title, subtitle = NULL, source_note = NULL) {
  gt_tbl <- df |>
    gt::gt() |>
    gt::tab_header(
      title = gt::md(title),
      subtitle = if (!is.null(subtitle)) gt::md(subtitle) else NULL
    ) |>
    gt::opt_table_font(font = gt::default_fonts()) |>
    gt::tab_options(
      heading.title.font.size = 14,
      heading.subtitle.font.size = 11,
      table.font.size = 11,
      data_row.padding = gt::px(6),
      source_notes.font.size = 10
    )
  
  if (!is.null(source_note)) {
    gt_tbl <- gt_tbl |> gt::tab_source_note(source_note = gt::md(source_note))
  }
  
  gt::gtsave(gt_tbl, filename = file.path(analysis_root, "tables", paste0(file_stub, ".html")))
  readr::write_csv(df, file.path(analysis_root, "tables", paste0(file_stub, ".csv")))
  invisible(gt_tbl)
}

# -----------------------------------------------------------------------------
# Custom weighted ROC / AUC helpers
# -----------------------------------------------------------------------------
weighted_roc_points <- function(y, p, w) {
  y <- as.numeric(y)
  p <- as.numeric(p)
  w <- as.numeric(w)
  
  ok <- !is.na(y) & !is.na(p) & !is.na(w) & is.finite(p) & is.finite(w)
  y <- y[ok]
  p <- p[ok]
  w <- w[ok]
  
  if (length(y) < 2 || length(unique(y)) < 2) {
    return(tibble::tibble(fpr = NA_real_, tpr = NA_real_))
  }
  
  ord <- order(-p, seq_along(p))
  y <- y[ord]
  p <- p[ord]
  w <- w[ord]
  
  pos_w <- sum(w[y == 1], na.rm = TRUE)
  neg_w <- sum(w[y == 0], na.rm = TRUE)
  
  if (pos_w <= 0 || neg_w <= 0) {
    return(tibble::tibble(fpr = NA_real_, tpr = NA_real_))
  }
  
  tp <- cumsum(ifelse(y == 1, w, 0))
  fp <- cumsum(ifelse(y == 0, w, 0))
  
  out <- tibble::tibble(
    score = p,
    fpr = fp / neg_w,
    tpr = tp / pos_w
  )
  
  keep <- c(TRUE, diff(out$score) != 0)
  out <- out[keep, c("fpr", "tpr"), drop = FALSE]
  
  out <- dplyr::bind_rows(
    tibble::tibble(fpr = 0, tpr = 0),
    out,
    tibble::tibble(fpr = 1, tpr = 1)
  )
  
  out
}

weighted_auc <- function(y, p, w) {
  roc <- weighted_roc_points(y, p, w)
  if (nrow(roc) < 2 || anyNA(roc$fpr) || anyNA(roc$tpr)) return(NA_real_)
  sum(diff(roc$fpr) * (head(roc$tpr, -1) + tail(roc$tpr, -1)) / 2)
}

roc_curve_df <- function(y, p, w, dataset_label, model_label) {
  roc <- weighted_roc_points(y, p, w)
  roc$dataset <- dataset_label
  roc$model <- model_label
  roc
}

performance_metrics <- function(y, p, w, dataset_label, model_label) {
  y <- as.numeric(y)
  p <- clip_prob(p)
  w <- as.numeric(w)
  
  auc <- weighted_auc(y, p, w)
  brier <- safe_weighted_mean((y - p)^2, w)
  
  cal_fit <- tryCatch(
    suppressWarnings(stats::glm(y ~ qlogis(p), family = quasibinomial(), weights = w)),
    error = function(e) NULL
  )
  
  int_fit <- tryCatch(
    suppressWarnings(stats::glm(y ~ 1, offset = qlogis(p), family = quasibinomial(), weights = w)),
    error = function(e) NULL
  )
  
  tibble::tibble(
    dataset = dataset_label,
    model = model_label,
    auc = auc,
    brier = brier,
    calibration_intercept = if (!is.null(int_fit)) unname(stats::coef(int_fit)[1]) else NA_real_,
    calibration_slope = if (!is.null(cal_fit)) unname(stats::coef(cal_fit)[2]) else NA_real_,
    mean_observed = safe_weighted_mean(y, w),
    mean_predicted = safe_weighted_mean(p, w)
  )
}

calibration_groups <- function(df, pred_var, y_var = "cs", w_var = "wt", g = 10, dataset_label = "Validation") {
  df |>
    dplyr::filter(!is.na(.data[[pred_var]]), !is.na(.data[[y_var]]), !is.na(.data[[w_var]])) |>
    dplyr::mutate(decile = dplyr::ntile(.data[[pred_var]], g)) |>
    dplyr::group_by(decile) |>
    dplyr::summarise(
      mean_predicted = safe_weighted_mean(.data[[pred_var]], .data[[w_var]]),
      mean_observed  = safe_weighted_mean(.data[[y_var]], .data[[w_var]]),
      weighted_n     = sum(.data[[w_var]], na.rm = TRUE),
      raw_n          = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::mutate(dataset = dataset_label)
}

group_prev_svy <- function(design, var, label) {
  out <- survey::svyby(
    ~cs,
    stats::as.formula(paste0("~", var)),
    design = design,
    FUN = survey::svymean,
    vartype = "se",
    na.rm = TRUE,
    drop.empty.groups = FALSE
  )
  
  se_name <- grep("^(se|se\\.)", names(out), value = TRUE)[1]
  se_vec <- if (!is.na(se_name) && nzchar(se_name)) out[[se_name]] else rep(NA_real_, nrow(out))
  rate <- out$cs
  
  tibble::tibble(
    Group = label,
    Level = as.character(out[[var]]),
    Rate = rate,
    LCL = pmax(0, rate - 1.96 * se_vec),
    UCL = pmin(1, rate + 1.96 * se_vec)
  )
}

# -----------------------------------------------------------------------------
# Table helpers
# -----------------------------------------------------------------------------
table1_continuous <- function(df, var, label) {
  d0 <- df |> dplyr::filter(cs == 0 & !is.na(.data[[var]]))
  d1 <- df |> dplyr::filter(cs == 1 & !is.na(.data[[var]]))
  da <- df |> dplyr::filter(!is.na(.data[[var]]))
  
  tibble::tibble(
    Characteristic = label,
    Level = "",
    Overall = sprintf("%s (%s)", fmt_num(safe_weighted_mean(da[[var]], da$wt), 2), fmt_num(safe_weighted_sd(da[[var]], da$wt), 2)),
    `No caesarean` = sprintf("%s (%s)", fmt_num(safe_weighted_mean(d0[[var]], d0$wt), 2), fmt_num(safe_weighted_sd(d0[[var]], d0$wt), 2)),
    Caesarean = sprintf("%s (%s)", fmt_num(safe_weighted_mean(d1[[var]], d1$wt), 2), fmt_num(safe_weighted_sd(d1[[var]], d1$wt), 2))
  )
}

table1_categorical <- function(df, var, label) {
  lvls <- sort(unique(as.character(df[[var]][!is.na(df[[var]])])))
  
  purrr::map_dfr(lvls, function(lvl) {
    da <- df |> dplyr::filter(!is.na(.data[[var]]))
    d0 <- df |> dplyr::filter(cs == 0 & !is.na(.data[[var]]))
    d1 <- df |> dplyr::filter(cs == 1 & !is.na(.data[[var]]))
    
    n_all <- sum(as.character(da[[var]]) == lvl, na.rm = TRUE)
    n_0   <- sum(as.character(d0[[var]]) == lvl, na.rm = TRUE)
    n_1   <- sum(as.character(d1[[var]]) == lvl, na.rm = TRUE)
    
    p_all <- safe_weighted_mean(as.numeric(as.character(da[[var]]) == lvl), da$wt)
    p_0   <- safe_weighted_mean(as.numeric(as.character(d0[[var]]) == lvl), d0$wt)
    p_1   <- safe_weighted_mean(as.numeric(as.character(d1[[var]]) == lvl), d1$wt)
    
    tibble::tibble(
      Characteristic = label,
      Level = lvl,
      Overall = sprintf("%s (%s)", fmt_n(n_all), fmt_pct(p_all, 1)),
      `No caesarean` = sprintf("%s (%s)", fmt_n(n_0), fmt_pct(p_0, 1)),
      Caesarean = sprintf("%s (%s)", fmt_n(n_1), fmt_pct(p_1, 1))
    )
  })
}

pretty_term <- function(x) {
  out <- x
  out <- gsub("^education_grp", "Education: ", out)
  out <- gsub("^wealth_grp", "Wealth: ", out)
  out <- gsub("^residence_grp", "Residence: ", out)
  out <- gsub("^multiple_gestation", "Multiple gestation: ", out)
  out <- gsub("^prior_cs", "Prior caesarean: ", out)
  out <- gsub("^delivery_sector_grp", "Delivery sector: ", out)
  out <- gsub("^anc_content_cat", "ANC content score: ", out)
  out <- gsub("^early_anc", "Early ANC: ", out)
  out <- gsub("^anc4plus", "ANC 4 or more visits: ", out)
  out <- gsub("^anc8plus", "ANC 8 or more visits: ", out)
  out <- gsub("^iron_use", "Iron in pregnancy: ", out)
  out <- gsub("^sp_malaria", "SP/Fansidar in pregnancy: ", out)
  out <- gsub("^age_imp_c$", "Maternal age (centered)", out)
  out <- gsub("^age_imp_sq$", "Maternal age squared", out)
  out <- gsub("^parity_imp_c$", "Parity (centered)", out)
  out <- gsub("^parity_imp_sq$", "Parity squared", out)
  out <- gsub("^bmi_imp_c$", "BMI (centered)", out)
  out <- gsub("^bmi_imp_sq$", "BMI squared", out)
  out <- gsub("^age_missing$", "Age missing", out)
  out <- gsub("^parity_missing$", "Parity missing", out)
  out <- gsub("^bmi_missing$", "BMI missing", out)
  out <- gsub("`", "", out, fixed = TRUE)
  out
}

# -----------------------------------------------------------------------------
# Read data
# -----------------------------------------------------------------------------
message("Reading BR file: ", br_path)
br <- haven::read_dta(br_path)

required_vars <- c(
  "v005", "v021", "v022", "m17",
  "v012", "v106", "v190", "v025", "v024", "v201", "v445",
  "b0", "v401", "m15", "m42c", "m42d", "m42e", "m42f"
)

missing_vars <- setdiff(required_vars, names(br))
if (length(missing_vars) > 0) {
  stop(sprintf("Missing required BR variables: %s", paste(missing_vars, collapse = ", ")), call. = FALSE)
}

# -----------------------------------------------------------------------------
# Raw derivation from BR
# -----------------------------------------------------------------------------
dat0 <- br |>
  dplyr::mutate(
    row_id = dplyr::row_number(),
    wt = v005 / 1e6,
    
    age_years = suppressWarnings(as.numeric(v012)),
    parity = suppressWarnings(as.numeric(v201)),
    bmi_raw_num = suppressWarnings(as.numeric(v445)),
    
    bmi = dplyr::case_when(
      is.na(bmi_raw_num) ~ NA_real_,
      bmi_raw_num %in% c(9998, 9999) ~ NA_real_,
      bmi_raw_num > 100 ~ bmi_raw_num / 100,
      TRUE ~ bmi_raw_num
    ),
    
    age_years = dplyr::if_else(!is.na(age_years) & (age_years < 10 | age_years > 55), NA_real_, age_years),
    parity    = dplyr::if_else(!is.na(parity)    & (parity < 0 | parity > 20), NA_real_, parity),
    bmi       = dplyr::if_else(!is.na(bmi)       & (bmi < 12 | bmi > 60), NA_real_, bmi),
    
    cs = dplyr::case_when(
      suppressWarnings(as.numeric(m17)) == 1 ~ 1,
      suppressWarnings(as.numeric(m17)) == 0 ~ 0,
      TRUE ~ NA_real_
    ),
    cs_factor = factor(cs, levels = c(0, 1), labels = c("No caesarean", "Caesarean")),
    
    education_raw = safe_as_factor(v106),
    wealth_raw = safe_as_factor(v190),
    residence_raw = safe_as_factor(v025),
    zone_raw = safe_as_factor(v024),
    
    multiple_gestation = factor(
      dplyr::case_when(
        is.na(suppressWarnings(as.numeric(b0))) ~ "Unknown",
        suppressWarnings(as.numeric(b0)) == 0 ~ "Singleton",
        suppressWarnings(as.numeric(b0)) > 0 ~ "Multiple gestation",
        TRUE ~ "Unknown"
      ),
      levels = c("Singleton", "Multiple gestation", "Unknown")
    ),
    
    prior_cs_num = suppressWarnings(as.numeric(v401)),
    prior_cs = factor(
      dplyr::case_when(
        prior_cs_num == 1 ~ "Yes",
        prior_cs_num == 0 ~ "No",
        is.na(prior_cs_num) & !is.na(parity) & parity <= 1 ~ "No",
        TRUE ~ "Unknown"
      ),
      levels = c("No", "Yes", "Unknown")
    ),
    
    anc_bp = yes_no_unknown_factor(m42c),
    anc_urine = yes_no_unknown_factor(m42d),
    anc_blood = yes_no_unknown_factor(m42e),
    anc_heartbeat = yes_no_unknown_factor(m42f),
    
    place_delivery_raw = safe_as_factor(m15),
    delivery_sector_raw = factor(
      dplyr::case_when(
        stringr::str_detect(stringr::str_to_lower(as.character(place_delivery_raw)), "home") ~ "Home",
        stringr::str_detect(stringr::str_to_lower(as.character(place_delivery_raw)), "private|mission|faith|ngo") ~ "Private facility",
        stringr::str_detect(stringr::str_to_lower(as.character(place_delivery_raw)), "public|government|teaching|federal|state|general hospital|health centre|health center|primary health|phc|dispensary|maternity") ~ "Public facility",
        is.na(place_delivery_raw) ~ "Unknown",
        TRUE ~ "Other/unspecified facility"
      ),
      levels = c("Home", "Public facility", "Private facility", "Other/unspecified facility", "Unknown")
    ),
    
    m13_num = if ("m13" %in% names(br)) suppressWarnings(as.numeric(m13)) else NA_real_,
    m14_num = if ("m14" %in% names(br)) suppressWarnings(as.numeric(m14)) else NA_real_,
    m45_num = if ("m45" %in% names(br)) suppressWarnings(as.numeric(m45)) else NA_real_,
    m49a_num = if ("m49a" %in% names(br)) suppressWarnings(as.numeric(m49a)) else NA_real_,
    
    valid_design = !is.na(v021) & !is.na(v022) & !is.na(wt) & wt > 0
  )

# -----------------------------------------------------------------------------
# Clinically reasoned phenotypes
# -----------------------------------------------------------------------------
dat0 <- dat0 |>
  dplyr::mutate(
    education_grp = factor(
      dplyr::case_when(
        stringr::str_detect(stringr::str_to_lower(as.character(education_raw)), "no education|none|primary") ~ "No/primary",
        stringr::str_detect(stringr::str_to_lower(as.character(education_raw)), "secondary") ~ "Secondary",
        stringr::str_detect(stringr::str_to_lower(as.character(education_raw)), "higher|more than secondary") ~ "Higher",
        is.na(education_raw) ~ "Unknown",
        TRUE ~ "Unknown"
      ),
      levels = c("No/primary", "Secondary", "Higher", "Unknown")
    ),
    
    wealth_grp = factor(
      dplyr::case_when(
        stringr::str_detect(stringr::str_to_lower(as.character(wealth_raw)), "poorest|poorer") ~ "Poor",
        stringr::str_detect(stringr::str_to_lower(as.character(wealth_raw)), "middle") ~ "Middle",
        stringr::str_detect(stringr::str_to_lower(as.character(wealth_raw)), "richer|richest") ~ "Rich",
        is.na(wealth_raw) ~ "Unknown",
        TRUE ~ "Unknown"
      ),
      levels = c("Poor", "Middle", "Rich", "Unknown")
    ),
    
    residence_grp = factor(
      dplyr::case_when(
        stringr::str_detect(stringr::str_to_lower(as.character(residence_raw)), "urban") ~ "Urban",
        stringr::str_detect(stringr::str_to_lower(as.character(residence_raw)), "rural") ~ "Rural",
        is.na(residence_raw) ~ "Unknown",
        TRUE ~ "Unknown"
      ),
      levels = c("Rural", "Urban", "Unknown")
    ),
    
    delivery_sector_grp = factor(
      dplyr::case_when(
        as.character(delivery_sector_raw) == "Home" ~ "Home",
        as.character(delivery_sector_raw) == "Public facility" ~ "Public facility",
        as.character(delivery_sector_raw) %in% c("Private facility", "Other/unspecified facility") ~ "Private/other facility",
        TRUE ~ "Unknown"
      ),
      levels = c("Home", "Public facility", "Private/other facility", "Unknown")
    ),
    
    anc_content_yes = rowSums(
      cbind(
        as.integer(anc_bp == "Yes"),
        as.integer(anc_urine == "Yes"),
        as.integer(anc_blood == "Yes"),
        as.integer(anc_heartbeat == "Yes")
      ),
      na.rm = TRUE
    ),
    
    anc_content_cat = factor(
      dplyr::case_when(
        anc_content_yes <= 1 ~ "0-1 components",
        anc_content_yes %in% c(2, 3) ~ "2-3 components",
        anc_content_yes == 4 ~ "4 components",
        TRUE ~ "Unknown"
      ),
      levels = c("0-1 components", "2-3 components", "4 components", "Unknown")
    ),
    
    early_anc = factor(
      dplyr::case_when(
        is.na(m13_num) ~ "Unknown",
        m13_num <= 3 ~ "Yes",
        m13_num > 3 & m13_num < 98 ~ "No",
        TRUE ~ "Unknown"
      ),
      levels = c("No", "Yes", "Unknown")
    ),
    
    anc4plus = factor(
      dplyr::case_when(
        is.na(m14_num) ~ "Unknown",
        m14_num >= 4 & m14_num < 98 ~ "Yes",
        m14_num < 4 ~ "No",
        TRUE ~ "Unknown"
      ),
      levels = c("No", "Yes", "Unknown")
    ),
    
    anc8plus = factor(
      dplyr::case_when(
        is.na(m14_num) ~ "Unknown",
        m14_num >= 8 & m14_num < 98 ~ "Yes",
        m14_num < 8 ~ "No",
        TRUE ~ "Unknown"
      ),
      levels = c("No", "Yes", "Unknown")
    ),
    
    iron_use = factor(
      dplyr::case_when(
        is.na(m45_num) ~ "Unknown",
        m45_num == 1 ~ "Yes",
        m45_num == 0 ~ "No",
        TRUE ~ "Unknown"
      ),
      levels = c("No", "Yes", "Unknown")
    ),
    
    sp_malaria = factor(
      dplyr::case_when(
        is.na(m49a_num) ~ "Unknown",
        m49a_num == 1 ~ "Yes",
        m49a_num == 0 ~ "No",
        TRUE ~ "Unknown"
      ),
      levels = c("No", "Yes", "Unknown")
    )
  )

# -----------------------------------------------------------------------------
# Cohort and missingness
# -----------------------------------------------------------------------------
analysis_vars_raw <- c(
  "age_years", "parity", "bmi",
  "education_grp", "wealth_grp", "residence_grp",
  "multiple_gestation", "prior_cs", "delivery_sector_grp",
  "anc_content_cat", "early_anc", "anc4plus", "anc8plus",
  "iron_use", "sp_malaria"
)

cohort_tbl <- tibble::tibble(
  Section = "Analytic cohort",
  Variable = c(
    "All birth records in BR file",
    "Non-missing outcome (m17)",
    "Valid survey design fields",
    "Final modelling cohort with valid design and outcome"
  ),
  Missing_n = c(NA, NA, NA, NA),
  Missing_pct = c(NA, NA, NA, NA),
  Raw_n = c(
    nrow(dat0),
    sum(!is.na(dat0$cs)),
    sum(!is.na(dat0$cs) & dat0$valid_design),
    sum(!is.na(dat0$cs) & dat0$valid_design)
  )
)

missing_tbl <- purrr::map_dfr(c("cs", analysis_vars_raw), function(v) {
  tibble::tibble(
    Section = "Variable missingness",
    Variable = v,
    Missing_n = sum(is.na(dat0[[v]]) & dat0$valid_design & !is.na(dat0$cs), na.rm = TRUE),
    Missing_pct = mean(is.na(dat0[[v]])[dat0$valid_design & !is.na(dat0$cs)], na.rm = TRUE),
    Raw_n = sum(!is.na(dat0[[v]]) & dat0$valid_design & !is.na(dat0$cs), na.rm = TRUE)
  )
})

table2 <- dplyr::bind_rows(cohort_tbl, missing_tbl) |>
  dplyr::mutate(
    Missing_n = dplyr::if_else(is.na(Missing_n), NA_character_, fmt_n(Missing_n)),
    Missing_pct = dplyr::if_else(is.na(Missing_pct), NA_character_, fmt_pct(Missing_pct, 1)),
    Raw_n = fmt_n(Raw_n)
  )

# -----------------------------------------------------------------------------
# Final analytic cohort
# -----------------------------------------------------------------------------
dat_analytic <- dat0 |>
  dplyr::filter(valid_design, !is.na(cs))

if (nrow(dat_analytic) < 1000) {
  stop("Analytic cohort is unexpectedly small after valid-design/non-missing outcome filtering.", call. = FALSE)
}

# -----------------------------------------------------------------------------
# PSU split
# -----------------------------------------------------------------------------
dat_split <- make_psu_split(dat_analytic, strata_var = "v022", psu_var = "v021", prop_train = 0.80, seed = 20260309)
dev <- dat_split |> dplyr::filter(split == "Development")
val <- dat_split |> dplyr::filter(split == "Validation")

if (nrow(dev) < 1000 || nrow(val) < 500) {
  stop("Development or validation sample is too small after PSU split.", call. = FALSE)
}

# -----------------------------------------------------------------------------
# Train-based continuous imputation + centering + quadratic terms
# -----------------------------------------------------------------------------
tmp <- impute_center_continuous(dev, val, "age_years", "age_imp_c", "age_missing")
dev <- tmp$train
val <- tmp$test

tmp <- impute_center_continuous(dev, val, "parity", "parity_imp_c", "parity_missing")
dev <- tmp$train
val <- tmp$test

tmp <- impute_center_continuous(dev, val, "bmi", "bmi_imp_c", "bmi_missing")
dev <- tmp$train
val <- tmp$test

dev$age_imp_sq <- dev$age_imp_c^2
val$age_imp_sq <- val$age_imp_c^2
dev$parity_imp_sq <- dev$parity_imp_c^2
val$parity_imp_sq <- val$parity_imp_c^2
dev$bmi_imp_sq <- dev$bmi_imp_c^2
val$bmi_imp_sq <- val$bmi_imp_c^2

# -----------------------------------------------------------------------------
# Train-based rare-level collapse / alignment
# -----------------------------------------------------------------------------
factor_model_vars <- c(
  "education_grp", "wealth_grp", "residence_grp",
  "multiple_gestation", "prior_cs", "delivery_sector_grp",
  "anc_content_cat", "early_anc", "anc4plus", "anc8plus",
  "iron_use", "sp_malaria"
)

for (v in factor_model_vars) {
  aligned <- collapse_with_train(dev[[v]], val[[v]], min_n = 50, other_label = "Other/unknown")
  dev[[v]] <- aligned$train
  val[[v]] <- aligned$test
}

# Preserve geography for reporting only
dev$zone_report <- as.character(dev$zone_raw)
val$zone_report <- as.character(val$zone_raw)
dev$zone_report[is.na(dev$zone_report) | !nzchar(dev$zone_report)] <- "Unknown"
val$zone_report[is.na(val$zone_report) | !nzchar(val$zone_report)] <- "Unknown"

# Standardised weights for model fitting
dev$wt_model <- dev$wt / mean(dev$wt, na.rm = TRUE)
val$wt_model <- val$wt / mean(dev$wt, na.rm = TRUE)

dat_model <- dplyr::bind_rows(dev, val) |>
  dplyr::arrange(row_id)

# -----------------------------------------------------------------------------
# Survey design object for descriptive results
# -----------------------------------------------------------------------------
des_all <- survey::svydesign(ids = ~v021, strata = ~v022, weights = ~wt, data = dat_model, nest = TRUE)

# -----------------------------------------------------------------------------
# Table 1 baseline
# -----------------------------------------------------------------------------
table1 <- dplyr::bind_rows(
  table1_continuous(dat_model, "age_years", "Maternal age, years"),
  table1_continuous(dat_model, "parity", "Children ever born (parity)"),
  table1_continuous(dat_model, "bmi", "BMI, kg/m^2"),
  table1_categorical(dat_model, "education_grp", "Education"),
  table1_categorical(dat_model, "wealth_grp", "Wealth"),
  table1_categorical(dat_model, "residence_grp", "Residence"),
  table1_categorical(dat_model, "multiple_gestation", "Multiple gestation"),
  table1_categorical(dat_model, "prior_cs", "Prior caesarean"),
  table1_categorical(dat_model, "delivery_sector_grp", "Delivery sector"),
  table1_categorical(dat_model, "anc_content_cat", "ANC content score"),
  table1_categorical(dat_model, "early_anc", "Early ANC"),
  table1_categorical(dat_model, "anc4plus", "ANC 4 or more visits"),
  table1_categorical(dat_model, "anc8plus", "ANC 8 or more visits"),
  table1_categorical(dat_model, "iron_use", "Iron use in pregnancy"),
  table1_categorical(dat_model, "sp_malaria", "SP/Fansidar in pregnancy")
)

# -----------------------------------------------------------------------------
# Prediction formulas
# -----------------------------------------------------------------------------
base_vars <- c(
  "age_imp_c", "age_imp_sq", "age_missing",
  "parity_imp_c", "parity_imp_sq", "parity_missing",
  "bmi_imp_c", "bmi_imp_sq", "bmi_missing",
  "education_grp", "wealth_grp", "residence_grp"
)

full_vars <- c(
  base_vars,
  "multiple_gestation", "prior_cs", "delivery_sector_grp",
  "anc_content_cat", "early_anc", "anc4plus", "anc8plus",
  "iron_use", "sp_malaria"
)

base_vars <- base_vars[vapply(base_vars, function(v) has_variation(dev[[v]]), logical(1))]
full_vars <- full_vars[vapply(full_vars, function(v) has_variation(dev[[v]]), logical(1))]

if (length(base_vars) == 0) stop("Base predictor list is empty after preprocessing.", call. = FALSE)
if (length(full_vars) == 0) stop("Full predictor list is empty after preprocessing.", call. = FALSE)

base_mm_formula <- stats::as.formula(paste("~", paste(base_vars, collapse = " + ")))
full_mm_formula <- stats::as.formula(paste("~", paste(full_vars, collapse = " + ")))

x_dev_base <- Matrix::sparse.model.matrix(base_mm_formula, data = dev)[, -1, drop = FALSE]
x_val_base <- Matrix::sparse.model.matrix(base_mm_formula, data = val)[, -1, drop = FALSE]

x_dev_full <- Matrix::sparse.model.matrix(full_mm_formula, data = dev)[, -1, drop = FALSE]
x_val_full <- Matrix::sparse.model.matrix(full_mm_formula, data = val)[, -1, drop = FALSE]

y_dev <- dev$cs
y_val <- val$cs

# -----------------------------------------------------------------------------
# Ridge prediction models
# -----------------------------------------------------------------------------
cv_base <- glmnet::cv.glmnet(
  x = x_dev_base,
  y = y_dev,
  family = "binomial",
  alpha = 0,
  weights = dev$wt_model,
  nfolds = 5,
  type.measure = "deviance",
  standardize = TRUE
)

cv_full <- glmnet::cv.glmnet(
  x = x_dev_full,
  y = y_dev,
  family = "binomial",
  alpha = 0,
  weights = dev$wt_model,
  nfolds = 5,
  type.measure = "deviance",
  standardize = TRUE
)

dev$pred_base <- as.numeric(stats::predict(cv_base, newx = x_dev_base, s = "lambda.1se", type = "response"))
val$pred_base <- as.numeric(stats::predict(cv_base, newx = x_val_base, s = "lambda.1se", type = "response"))

dev$pred_full <- as.numeric(stats::predict(cv_full, newx = x_dev_full, s = "lambda.1se", type = "response"))
val$pred_full <- as.numeric(stats::predict(cv_full, newx = x_val_full, s = "lambda.1se", type = "response"))

dev$pred_base <- clip_prob(dev$pred_base)
val$pred_base <- clip_prob(val$pred_base)
dev$pred_full <- clip_prob(dev$pred_full)
val$pred_full <- clip_prob(val$pred_full)

dat_model <- dplyr::bind_rows(dev, val) |>
  dplyr::arrange(row_id)

# -----------------------------------------------------------------------------
# Reduced weighted logistic refit for interpretable OR table / figure
# -----------------------------------------------------------------------------
effect_vars <- c(
  "age_imp_c", "age_imp_sq", "age_missing",
  "parity_imp_c", "parity_imp_sq", "parity_missing",
  "bmi_imp_c", "bmi_imp_sq", "bmi_missing",
  "education_grp", "wealth_grp", "residence_grp",
  "multiple_gestation", "prior_cs", "delivery_sector_grp",
  "anc_content_cat", "early_anc", "anc4plus", "iron_use"
)
effect_vars <- effect_vars[vapply(effect_vars, function(v) has_variation(dev[[v]]), logical(1))]

fit_effect <- NULL
if (length(effect_vars) > 0) {
  effect_formula <- stats::as.formula(paste("cs ~", paste(effect_vars, collapse = " + ")))
  fit_effect <- tryCatch(
    suppressWarnings(
      stats::glm(
        formula = effect_formula,
        data = dev,
        family = quasibinomial(),
        weights = wt_model,
        control = stats::glm.control(maxit = 100, epsilon = 1e-8)
      )
    ),
    error = function(e) NULL
  )
}

# -----------------------------------------------------------------------------
# Table 3 predictor display
# -----------------------------------------------------------------------------
ridge_coef <- as.matrix(stats::coef(cv_full, s = "lambda.1se"))
table3_ridge <- tibble::tibble(
  term = rownames(ridge_coef),
  beta = as.numeric(ridge_coef[, 1])
) |>
  dplyr::filter(term != "(Intercept)") |>
  dplyr::mutate(
    Predictor = pretty_term(term),
    `Penalized OR` = fmt_num(exp(beta), 3),
    `Coefficient` = fmt_num(beta, 4)
  ) |>
  dplyr::arrange(dplyr::desc(abs(beta))) |>
  dplyr::select(Predictor, `Penalized OR`, `Coefficient`)

if (!is.null(fit_effect)) {
  cf <- stats::coef(fit_effect)
  vv <- tryCatch(stats::vcov(fit_effect), error = function(e) NULL)
  
  if (!is.null(vv) && all(is.finite(cf))) {
    se <- sqrt(diag(vv))
    z <- cf / se
    p <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)
    lcl <- cf - 1.96 * se
    ucl <- cf + 1.96 * se
    
    table3 <- tibble::tibble(
      term = names(cf),
      estimate = cf,
      conf.low = lcl,
      conf.high = ucl,
      p.value = p
    ) |>
      dplyr::filter(term != "(Intercept)") |>
      dplyr::mutate(
        Predictor = pretty_term(term),
        `Adjusted OR (95% CI)` = sprintf("%s (%s to %s)", fmt_num(exp(estimate), 2), fmt_num(exp(conf.low), 2), fmt_num(exp(conf.high), 2)),
        `p-value` = fmt_p(p.value)
      ) |>
      dplyr::select(Predictor, `Adjusted OR (95% CI)`, `p-value`)
  } else {
    table3 <- table3_ridge
  }
} else {
  table3 <- table3_ridge
}

# -----------------------------------------------------------------------------
# Table 4 performance
# -----------------------------------------------------------------------------
table4_raw <- dplyr::bind_rows(
  performance_metrics(dev$cs, dev$pred_base, dev$wt, "Development", "Base model"),
  performance_metrics(val$cs, val$pred_base, val$wt, "Validation", "Base model"),
  performance_metrics(dev$cs, dev$pred_full, dev$wt, "Development", "Full model"),
  performance_metrics(val$cs, val$pred_full, val$wt, "Validation", "Full model")
)

table4 <- table4_raw |>
  dplyr::mutate(
    auc = fmt_num(auc, 3),
    brier = fmt_num(brier, 4),
    calibration_intercept = fmt_num(calibration_intercept, 3),
    calibration_slope = fmt_num(calibration_slope, 3),
    mean_observed = fmt_pct(mean_observed, 2),
    mean_predicted = fmt_pct(mean_predicted, 2)
  )

# -----------------------------------------------------------------------------
# Table 5 geography observed vs predicted
# -----------------------------------------------------------------------------
table5 <- dat_model |>
  dplyr::group_by(zone_report) |>
  dplyr::summarise(
    Raw_n = dplyr::n(),
    Weighted_births = sum(wt, na.rm = TRUE),
    Observed_CS = safe_weighted_mean(cs, wt),
    Predicted_CS = safe_weighted_mean(pred_full, wt),
    Expected_CS_per_1000 = 1000 * safe_weighted_mean(pred_full, wt),
    .groups = "drop"
  ) |>
  dplyr::arrange(dplyr::desc(Predicted_CS)) |>
  dplyr::mutate(
    Raw_n = fmt_n(Raw_n),
    Weighted_births = fmt_n(Weighted_births, 0),
    Observed_CS = fmt_pct(Observed_CS, 1),
    Predicted_CS = fmt_pct(Predicted_CS, 1),
    Expected_CS_per_1000 = fmt_num(Expected_CS_per_1000, 1)
  ) |>
  dplyr::rename(`Geographic unit (v024 label)` = zone_report)

# -----------------------------------------------------------------------------
# Figure 1 cohort flow
# -----------------------------------------------------------------------------
fig1_df <- tibble::tibble(
  Stage = factor(
    c("All BR records", "Non-missing m17", "Valid survey design", "Final analytic cohort", "Development set", "Validation set"),
    levels = c("All BR records", "Non-missing m17", "Valid survey design", "Final analytic cohort", "Development set", "Validation set")
  ),
  N = c(
    nrow(dat0),
    sum(!is.na(dat0$cs)),
    sum(!is.na(dat0$cs) & dat0$valid_design),
    nrow(dat_analytic),
    nrow(dev),
    nrow(val)
  )
)

fig1 <- ggplot2::ggplot(fig1_df, ggplot2::aes(x = Stage, y = N)) +
  ggplot2::geom_col(fill = "#2F6B7C", width = 0.72) +
  ggplot2::geom_text(ggplot2::aes(label = fmt_n(N)), vjust = -0.25, size = 3.7, fontface = "bold", colour = "#102A43") +
  ggplot2::scale_y_continuous(labels = scales::comma, expand = ggplot2::expansion(mult = c(0, 0.08))) +
  ggplot2::labs(
    title = "Analytic cohort derivation for NDHS 2024 caesarean-section risk study",
    subtitle = "Counts are raw numbers of birth records retained at each stage",
    x = NULL,
    y = "Number of birth records",
    caption = "Primary analytic file: NGBR8AFL.dta"
  ) +
  ndhs_theme()

save_png(fig1, "fig1_cohort_flow.png", width = 11, height = 7)

# -----------------------------------------------------------------------------
# Figure 2 weighted prevalence by key groups
# -----------------------------------------------------------------------------
fig2_df <- dplyr::bind_rows(
  group_prev_svy(des_all, "residence_grp", "Residence"),
  group_prev_svy(des_all, "wealth_grp", "Wealth"),
  group_prev_svy(des_all, "delivery_sector_grp", "Delivery sector"),
  group_prev_svy(des_all, "prior_cs", "Prior caesarean")
) |>
  dplyr::mutate(
    Level = forcats::fct_inorder(Level),
    Rate_pct = 100 * Rate,
    LCL_pct = 100 * LCL,
    UCL_pct = 100 * UCL
  )

fig2 <- ggplot2::ggplot(fig2_df, ggplot2::aes(x = forcats::fct_reorder(Level, Rate_pct), y = Rate_pct)) +
  ggplot2::geom_col(fill = "#3E8E7E", width = 0.72) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = LCL_pct, ymax = UCL_pct), width = 0.18, colour = "#1F2933") +
  ggplot2::coord_flip() +
  ggplot2::facet_wrap(~ Group, scales = "free_y", ncol = 2) +
  ggplot2::scale_y_continuous(labels = function(x) paste0(x, "%"), expand = ggplot2::expansion(mult = c(0, 0.08))) +
  ggplot2::labs(
    title = "Weighted caesarean-section prevalence across key maternal and service strata",
    subtitle = "Bars show survey-weighted prevalence; whiskers show approximate 95% confidence intervals",
    x = NULL,
    y = "Weighted prevalence of caesarean section",
    caption = "Descriptive estimates use the DHS complex survey design."
  ) +
  ndhs_theme()

save_png(fig2, "fig2_weighted_prevalence_key_groups.png", width = 12, height = 9)

# -----------------------------------------------------------------------------
# Figure 3 predictor plot
# -----------------------------------------------------------------------------
fig3 <- NULL

if (!is.null(fit_effect)) {
  cf <- stats::coef(fit_effect)
  vv <- tryCatch(stats::vcov(fit_effect), error = function(e) NULL)
  
  if (!is.null(vv) && all(is.finite(cf))) {
    se <- sqrt(diag(vv))
    lcl <- cf - 1.96 * se
    ucl <- cf + 1.96 * se
    
    plot_df <- tibble::tibble(
      term = names(cf),
      estimate = cf,
      conf.low = lcl,
      conf.high = ucl
    ) |>
      dplyr::filter(term != "(Intercept)") |>
      dplyr::mutate(
        Predictor = pretty_term(term),
        OR = exp(estimate),
        LCL = exp(conf.low),
        UCL = exp(conf.high)
      ) |>
      dplyr::filter(is.finite(OR), is.finite(LCL), is.finite(UCL)) |>
      dplyr::arrange(OR)
    
    plot_df <- safe_tail_n(plot_df, 25)
    
    if (nrow(plot_df) > 0) {
      plot_df <- plot_df |>
        dplyr::mutate(Predictor = factor(Predictor, levels = Predictor))
      
      fig3 <- ggplot2::ggplot(plot_df, ggplot2::aes(x = Predictor, y = OR, ymin = LCL, ymax = UCL)) +
        ggplot2::geom_hline(yintercept = 1, linetype = 2, colour = "#7B8794") +
        ggplot2::geom_pointrange(colour = "#7C3AED", linewidth = 0.4) +
        ggplot2::coord_flip() +
        ggplot2::scale_y_log10() +
        ggplot2::labs(
          title = "Adjusted predictor profile from reduced weighted logistic refit",
          subtitle = "Top displayed terms by magnitude; prediction engine remains weighted ridge logistic regression",
          x = NULL,
          y = "Adjusted odds ratio (log scale)",
          caption = "Interpretation plot for manuscript presentation."
        ) +
        ndhs_theme()
    }
  }
}

if (is.null(fig3)) {
  ridge_plot_df <- tibble::tibble(
    term = rownames(ridge_coef),
    beta = as.numeric(ridge_coef[, 1])
  ) |>
    dplyr::filter(term != "(Intercept)") |>
    dplyr::mutate(
      Predictor = pretty_term(term),
      Penalized_OR = exp(beta)
    ) |>
    dplyr::filter(is.finite(Penalized_OR)) |>
    dplyr::arrange(abs(beta))
  
  ridge_plot_df <- safe_tail_n(ridge_plot_df, 25)
  
  if (nrow(ridge_plot_df) > 0) {
    ridge_plot_df <- ridge_plot_df |>
      dplyr::mutate(Predictor = factor(Predictor, levels = Predictor))
    
    fig3 <- ggplot2::ggplot(ridge_plot_df, ggplot2::aes(x = Predictor, y = Penalized_OR)) +
      ggplot2::geom_hline(yintercept = 1, linetype = 2, colour = "#7B8794") +
      ggplot2::geom_point(colour = "#7C3AED", size = 2.4) +
      ggplot2::coord_flip() +
      ggplot2::scale_y_log10() +
      ggplot2::labs(
        title = "Penalized predictor profile from the full ridge model",
        subtitle = "Top displayed terms by absolute coefficient magnitude",
        x = NULL,
        y = "Penalized odds ratio (log scale)",
        caption = "These are shrunken model coefficients from the prediction model, not causal effects."
      ) +
      ndhs_theme()
  } else {
    fig3 <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 1, y = 1, label = "No predictor plot available") +
      ggplot2::theme_void()
  }
}

save_png(fig3, "fig3_predictor_plot.png", width = 12, height = 10)

# -----------------------------------------------------------------------------
# Figure 4 ROC base vs full on validation
# -----------------------------------------------------------------------------
roc_df <- dplyr::bind_rows(
  roc_curve_df(y_val, val$pred_base, val$wt, "Validation", "Base model"),
  roc_curve_df(y_val, val$pred_full, val$wt, "Validation", "Full model")
)

perf_val <- dplyr::bind_rows(
  performance_metrics(y_val, val$pred_base, val$wt, "Validation", "Base model"),
  performance_metrics(y_val, val$pred_full, val$wt, "Validation", "Full model")
)

sub_lab <- paste0(
  "Validation AUCs — Base: ", fmt_num(perf_val$auc[perf_val$model == "Base model"], 3),
  " | Full: ", fmt_num(perf_val$auc[perf_val$model == "Full model"], 3)
)

fig4 <- ggplot2::ggplot(roc_df, ggplot2::aes(x = fpr, y = tpr, colour = model)) +
  ggplot2::geom_abline(intercept = 0, slope = 1, linetype = 2, colour = "#9FB3C8") +
  ggplot2::geom_line(linewidth = 1.1) +
  ggplot2::scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  ggplot2::scale_colour_manual(values = c("Base model" = "#1D4ED8", "Full model" = "#DC2626")) +
  ggplot2::labs(
    title = "Validation discrimination: base versus full caesarean-risk models",
    subtitle = sub_lab,
    x = "False-positive rate",
    y = "True-positive rate",
    colour = NULL,
    caption = "Prediction uses weighted ridge logistic regression fitted in the development subset."
  ) +
  ndhs_theme()

save_png(fig4, "fig4_validation_roc_base_vs_full.png", width = 10, height = 7)

# -----------------------------------------------------------------------------
# Figure 5 calibration full model
# -----------------------------------------------------------------------------
cal_df <- calibration_groups(val, pred_var = "pred_full", y_var = "cs", w_var = "wt", g = 10, dataset_label = "Validation")
max_cal <- suppressWarnings(max(c(cal_df$mean_predicted, cal_df$mean_observed), na.rm = TRUE))
if (!is.finite(max_cal) || max_cal <= 0) max_cal <- 0.1

fig5 <- ggplot2::ggplot(cal_df, ggplot2::aes(x = mean_predicted, y = mean_observed)) +
  ggplot2::geom_abline(intercept = 0, slope = 1, linetype = 2, colour = "#9FB3C8") +
  ggplot2::geom_point(size = 2.8, colour = "#DC2626") +
  ggplot2::geom_line(linewidth = 0.8, colour = "#DC2626") +
  ggplot2::scale_x_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, max_cal * 1.05)
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, max_cal * 1.05)
  ) +
  ggplot2::labs(
    title = "Validation calibration of the full caesarean-risk model",
    subtitle = "Observed versus mean predicted risk across deciles of predicted probability",
    x = "Mean predicted risk",
    y = "Observed weighted risk",
    caption = "Perfect calibration lies on the 45-degree line."
  ) +
  ndhs_theme()

save_png(fig5, "fig5_validation_calibration_full_model.png", width = 10, height = 7)

# -----------------------------------------------------------------------------
# Save tables
# -----------------------------------------------------------------------------
save_gt(
  table1,
  "table1_baseline_characteristics",
  title = "**Table 1. Baseline maternal, antenatal, and delivery characteristics by caesarean-section status**",
  subtitle = "Categorical variables are raw n with weighted percentage; continuous variables are weighted mean (SD).",
  source_note = "Descriptive summaries are based on the final analytic cohort with valid outcome and survey-design variables."
)

save_gt(
  table2,
  "table2_missingness_and_cohort",
  title = "**Table 2. Cohort derivation and variable missingness**",
  subtitle = "Missingness is evaluated among records with valid survey design and non-missing outcome.",
  source_note = "Prediction modeling uses train-based imputation for continuous variables and explicit Unknown categories for categorical variables."
)

save_gt(
  table3,
  "table3_model_predictors",
  title = "**Table 3. Predictor profile for the final model**",
  subtitle = "Reduced weighted logistic refit is shown when stable; otherwise penalized coefficients from the ridge model are displayed.",
  source_note = "The primary prediction engine is weighted ridge logistic regression."
)

save_gt(
  table4,
  "table4_model_performance",
  title = "**Table 4. Model discrimination and calibration for base and full models**",
  subtitle = "Performance is shown in development and validation subsets.",
  source_note = "Base model = sociodemographic and maternal background features; full model = base model plus obstetric and antenatal phenotypes."
)

save_gt(
  table5,
  "table5_observed_vs_predicted_by_geography",
  title = "**Table 5. Observed versus predicted caesarean-section burden by geographic unit**",
  subtitle = "Predicted risk comes from the full weighted ridge logistic model; geography is used for reporting only.",
  source_note = "Geographic variable corresponds to the v024 label in the BR file."
)

# -----------------------------------------------------------------------------
# Save models and supplementary outputs
# -----------------------------------------------------------------------------
saveRDS(cv_base, file.path(analysis_root, "models", "cv_base_ridge_glmnet.rds"))
saveRDS(cv_full, file.path(analysis_root, "models", "cv_full_ridge_glmnet.rds"))
if (!is.null(fit_effect)) saveRDS(fit_effect, file.path(analysis_root, "models", "reduced_weighted_logistic_refit.rds"))

readr::write_csv(dev, file.path(analysis_root, "supplementary", "development_predictions.csv"))
readr::write_csv(val, file.path(analysis_root, "supplementary", "validation_predictions.csv"))

ridge_coef_tbl <- tibble::tibble(
  term = rownames(ridge_coef),
  beta = as.numeric(ridge_coef[, 1])
)
readr::write_csv(ridge_coef_tbl, file.path(analysis_root, "supplementary", "full_ridge_coefficients.csv"))

phenotype_dictionary <- tibble::tribble(
  ~phenotype, ~source_variables, ~definition,
  "cs", "m17", "Delivery by caesarean section; binary outcome",
  "multiple_gestation", "b0", "Singleton / multiple gestation / unknown",
  "prior_cs", "v401, v201", "Previous caesarean where known; parity<=1 with missing v401 classified as No",
  "delivery_sector_grp", "m15", "Home / public facility / private-other facility / unknown",
  "anc_content_cat", "m42c,m42d,m42e,m42f", "Count of four ANC content components: 0-1 / 2-3 / 4 / unknown",
  "early_anc", "m13", "First ANC by 3 months versus later versus unknown",
  "anc4plus", "m14", "Four or more ANC visits versus fewer versus unknown",
  "anc8plus", "m14", "Eight or more ANC visits versus fewer versus unknown",
  "iron_use", "m45", "Iron tablets/syrup in pregnancy yes/no/unknown",
  "sp_malaria", "m49a", "SP/Fansidar in pregnancy yes/no/unknown",
  "education_grp", "v106", "No-primary / secondary / higher / unknown",
  "wealth_grp", "v190", "Poor / middle / rich / unknown",
  "residence_grp", "v025", "Rural / urban / unknown",
  "age_imp_c", "v012", "Train-median-imputed and centered maternal age",
  "parity_imp_c", "v201", "Train-median-imputed and centered parity",
  "bmi_imp_c", "v445", "Train-median-imputed and centered BMI"
)
readr::write_csv(phenotype_dictionary, file.path(analysis_root, "supplementary", "phenotype_dictionary.csv"))

# -----------------------------------------------------------------------------
# Session log
# -----------------------------------------------------------------------------
sink(file.path(analysis_root, "logs", "sessionInfo.txt"))
print(sessionInfo())
sink()

message("Done. Outputs written to: ", analysis_root)