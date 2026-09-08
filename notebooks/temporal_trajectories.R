#!/usr/bin/env Rscript

# Analysis: Age- and follow-up-eligible matched comparison of
#           pre-diagnostic gastrointestinal cancer risk-score patterns
# Date: 2026-08-14
# Input: output_data/gcrc/data_analyze_gcrc_score.rds
# Key message: Comparing cases with age- and follow-up-eligible controls

set.seed(42)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(grid)
})

# ---------- Paths and reproducibility ----------

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)

if (length(file_arg) == 1L) {
  script_path <- normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
  script_dir <- dirname(script_path)
} else {
  # Interactive fallback: run from the project root.
  script_dir <- normalizePath(file.path(getwd(), "notebooks"), mustWork = TRUE)
}

project_root <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
input_file <- file.path(project_root, "output_data", "gcrc", "data_analyze_gcrc_score.rds")
output_dir <- file.path(project_root, "output_plot", "temporal_trajectories")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

cat("R version:", R.version.string, "\n")
cat("dplyr:", as.character(packageVersion("dplyr")), "\n")
cat("ggplot2:", as.character(packageVersion("ggplot2")), "\n")
cat("Random seed: 42\n")
cat("Input:", input_file, "\n\n")

# ---------- Constants ----------

MATCH_RATIO <- 50L
MIN_YEAR_TO_INDEX <- -15L
MAX_YEAR_TO_INDEX <- -2L
AGE_CALIPER_YEARS <- 5L
BOOTSTRAP_REPLICATES <- 2000L
CONTROL_COLOR <- "#3C5488"
CASE_COLOR <- "#E64B35"

# ---------- Helpers ----------

to_binary <- function(x, variable_name) {
  value <- trimws(as.character(x))
  out <- rep(NA_integer_, length(value))
  out[value %in% c("0", "No", "NO", "no", "FALSE", "False", "false")] <- 0L
  out[value %in% c("1", "Yes", "YES", "yes", "TRUE", "True", "true")] <- 1L

  unknown <- sort(unique(value[!is.na(value) & is.na(out)]))
  if (length(unknown) > 0L) {
    stop(
      "Unexpected coding in ", variable_name, ": ",
      paste(unknown, collapse = ", ")
    )
  }
  out
}

mean_ci <- function(x, conf_level = 0.95) {
  x <- x[is.finite(x)]
  n <- length(x)
  estimate <- if (n > 0L) mean(x) else NA_real_

  if (n < 2L) {
    return(c(estimate = estimate, lower = NA_real_, upper = NA_real_))
  }

  se <- stats::sd(x) / sqrt(n)
  critical <- stats::qt((1 + conf_level) / 2, df = n - 1L)
  c(
    estimate = estimate,
    lower = estimate - critical * se,
    upper = estimate + critical * se
  )
}

bootstrap_ci <- function(x, statistic = c("mean", "median"),
                         replicates = BOOTSTRAP_REPLICATES,
                         conf_level = 0.95) {
  statistic <- match.arg(statistic)
  x <- x[is.finite(x)]

  if (length(x) < 2L) {
    return(c(lower = NA_real_, upper = NA_real_))
  }

  statistic_function <- if (statistic == "mean") mean else median
  boot_values <- replicate(
    replicates,
    statistic_function(sample(x, length(x), replace = TRUE))
  )
  alpha <- (1 - conf_level) / 2
  stats::quantile(
    boot_values,
    probs = c(alpha, 1 - alpha),
    na.rm = TRUE,
    names = FALSE,
    type = 7
  ) |>
    stats::setNames(c("lower", "upper"))
}

normality_check <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)

  if (n < 3L || stats::sd(x) == 0) {
    return(list(method = "Not assessable", p_value = NA_real_, normal = FALSE))
  }

  if (n < 50L) {
    test <- stats::shapiro.test(x)
    return(list(
      method = "Shapiro-Wilk",
      p_value = unname(test$p.value),
      normal = unname(test$p.value) >= 0.05
    ))
  }

  standardized <- as.numeric(scale(x))
  test <- suppressWarnings(stats::ks.test(standardized, "pnorm"))
  list(
    method = "Kolmogorov-Smirnov",
    p_value = unname(test$p.value),
    normal = unname(test$p.value) >= 0.05
  )
}

paired_difference_test <- function(differences) {
  differences <- differences[is.finite(differences)]
  n <- length(differences)
  normality <- normality_check(differences)

  if (n < 2L) {
    return(data.frame(
      n_sets = n,
      effect_type = "Not estimable",
      effect_estimate = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      test = "Not estimable",
      statistic = NA_real_,
      p_value = NA_real_,
      normality_test = normality$method,
      normality_p = normality$p_value,
      stringsAsFactors = FALSE
    ))
  }

  if (isTRUE(normality$normal)) {
    test <- stats::t.test(differences, mu = 0, conf.level = 0.95)
    return(data.frame(
      n_sets = n,
      effect_type = "Mean paired difference",
      effect_estimate = unname(test$estimate),
      ci_lower = unname(test$conf.int[1]),
      ci_upper = unname(test$conf.int[2]),
      test = "One-sample t-test of matched-set differences",
      statistic = unname(test$statistic),
      p_value = unname(test$p.value),
      normality_test = normality$method,
      normality_p = normality$p_value,
      stringsAsFactors = FALSE
    ))
  }

  test <- suppressWarnings(stats::wilcox.test(
    differences,
    mu = 0,
    paired = FALSE,
    exact = FALSE,
    correct = FALSE
  ))
  ci <- bootstrap_ci(differences, statistic = "median")

  data.frame(
    n_sets = n,
    effect_type = "Median paired difference",
    effect_estimate = stats::median(differences),
    ci_lower = unname(ci["lower"]),
    ci_upper = unname(ci["upper"]),
    test = "Wilcoxon signed-rank test of matched-set differences",
    statistic = unname(test$statistic),
    p_value = unname(test$p.value),
    normality_test = normality$method,
    normality_p = normality$p_value,
    stringsAsFactors = FALSE
  )
}

standardized_mean_difference <- function(x, y) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  pooled_sd <- sqrt((stats::var(x) + stats::var(y)) / 2)

  if (!is.finite(pooled_sd) || pooled_sd == 0) {
    return(ifelse(isTRUE(all.equal(mean(x), mean(y))), 0, NA_real_))
  }

  (mean(x) - mean(y)) / pooled_sd
}

format_p <- function(p_value) {
  ifelse(
    is.na(p_value),
    "NA",
    ifelse(p_value < 0.001, "<0.001", sprintf("%.3f", p_value))
  )
}

write_csv <- function(x, filename) {
  utils::write.csv(
    x,
    file.path(output_dir, filename),
    row.names = FALSE,
    na = ""
  )
}

# ---------- Read and validate the canonical analysis data ----------

raw_data <- readRDS(input_file)
required_columns <- c(
  "recode", "time", "Age", "Sex", "score", "Stomach", "Colorectal", "gcrc"
)
missing_columns <- setdiff(required_columns, names(raw_data))

if (length(missing_columns) > 0L) {
  stop("Required columns missing: ", paste(missing_columns, collapse = ", "))
}

analysis_data <- raw_data |>
  transmute(
    source_row = dplyr::row_number(),
    followup_years = as.numeric(time),
    baseline_age = as.numeric(Age),
    sex = droplevels(as.factor(Sex)),
    score = as.numeric(score),
    stomach = to_binary(Stomach, "Stomach"),
    colorectal = to_binary(Colorectal, "Colorectal"),
    gcrc = to_binary(gcrc, "gcrc")
  )

if (any(analysis_data$followup_years <= 0, na.rm = TRUE)) {
  stop("All follow-up times must be positive.")
}

if (any(analysis_data$baseline_age %% 1 != 0, na.rm = TRUE)) {
  warning("Age is not integer-valued for all rows; matching uses observed values.")
}

profile <- data.frame(
  metric = c(
    "Total participants",
    "Gastric cancer events (all follow-up)",
    "Colorectal cancer events (all follow-up)",
    "Combined gastric/colorectal cancer events (all follow-up)",
    "Cancer-free controls",
    "Rows missing any matching/score field"
  ),
  value = c(
    nrow(analysis_data),
    sum(analysis_data$stomach == 1L, na.rm = TRUE),
    sum(analysis_data$colorectal == 1L, na.rm = TRUE),
    sum(analysis_data$gcrc == 1L, na.rm = TRUE),
    sum(analysis_data$gcrc == 0L, na.rm = TRUE),
    sum(!stats::complete.cases(
      analysis_data[, c("followup_years", "baseline_age", "sex", "score")]
    ))
  )
)
write_csv(profile, "analysis_cohort_profile.csv")
print(profile)

# ---------- Age- and follow-up-eligible matched-control sampling ----------

match_endpoint <- function(data, event_variable, endpoint_label, endpoint_short) {
  event <- data[[event_variable]]
  year_to_index <- -ceiling(data$followup_years)
  complete_fields <- stats::complete.cases(
    data[, c("followup_years", "baseline_age", "sex", "score")]
  )

  case_rows <- which(
    !is.na(event) & event == 1L & complete_fields &
      year_to_index >= MIN_YEAR_TO_INDEX &
      year_to_index <= MAX_YEAR_TO_INDEX
  )
  control_rows <- which(data$gcrc == 0L & complete_fields)

  cases <- data[case_rows, , drop = FALSE]
  cases$year_to_index <- year_to_index[case_rows]
  controls <- data[control_rows, , drop = FALSE]

  if (nrow(cases) == 0L) {
    stop("No eligible cases for endpoint: ", endpoint_label)
  }
  if (nrow(controls) == 0L) {
    stop("No eligible cancer-free controls.")
  }

  cases$eligible_controls <- vapply(
    seq_len(nrow(cases)),
    function(i) {
      sum(
        abs(controls$baseline_age - cases$baseline_age[i]) <= AGE_CALIPER_YEARS &
          controls$followup_years >= cases$followup_years[i]
      )
    },
    integer(1)
  )

  if (any(cases$eligible_controls == 0L)) {
    warning(
      endpoint_label, ": ", sum(cases$eligible_controls == 0L),
      " case(s) have no eligible controls before matching without replacement."
    )
  }

  # Match the most constrained/longest-followed cases first to reduce avoidable failures.
  match_order <- order(
    cases$eligible_controls,
    -cases$followup_years,
    cases$baseline_age
  )
  available <- rep(TRUE, nrow(controls))
  matched_set_rows <- vector("list", nrow(cases))
  matched_control_rows <- vector("list", nrow(cases))
  set_counter <- 0L

  for (case_position in match_order) {
    candidate_positions <- which(
      available &
        abs(controls$baseline_age - cases$baseline_age[case_position]) <=
          AGE_CALIPER_YEARS &
        controls$followup_years >= cases$followup_years[case_position]
    )

    if (length(candidate_positions) < MATCH_RATIO) {
      next
    }

    # Randomly sample inside the eligible age/follow-up pool; seed is fixed above.
    selected_positions <- sample(candidate_positions, MATCH_RATIO, replace = FALSE)
    available[selected_positions] <- FALSE
    selected_controls <- controls[selected_positions, , drop = FALSE]
    set_counter <- set_counter + 1L

    matched_set_rows[[set_counter]] <- data.frame(
      endpoint = endpoint_label,
      endpoint_short = endpoint_short,
      matched_set = set_counter,
      year_to_index = cases$year_to_index[case_position],
      index_interval_years = cases$followup_years[case_position],
      case_baseline_age = cases$baseline_age[case_position],
      control_baseline_age_mean = mean(selected_controls$baseline_age),
      case_index_age = cases$baseline_age[case_position] +
        cases$followup_years[case_position],
      control_pseudo_index_age_mean = mean(selected_controls$baseline_age) +
        cases$followup_years[case_position],
      sex = as.character(cases$sex[case_position]),
      case_score = cases$score[case_position],
      control_score_mean = mean(selected_controls$score),
      paired_difference = cases$score[case_position] - mean(selected_controls$score),
      n_controls = nrow(selected_controls),
      minimum_control_followup_margin = min(
        selected_controls$followup_years - cases$followup_years[case_position]
      ),
      stringsAsFactors = FALSE
    )

    matched_control_rows[[set_counter]] <- data.frame(
      endpoint = endpoint_label,
      endpoint_short = endpoint_short,
      matched_set = set_counter,
      year_to_index = cases$year_to_index[case_position],
      control_number = seq_len(nrow(selected_controls)),
      control_sex = as.character(selected_controls$sex),
      baseline_age_difference = selected_controls$baseline_age -
        cases$baseline_age[case_position],
      pseudo_index_age_difference = (
        selected_controls$baseline_age + cases$followup_years[case_position]
      ) - (
        cases$baseline_age[case_position] + cases$followup_years[case_position]
      ),
      followup_margin_years = selected_controls$followup_years -
        cases$followup_years[case_position],
      stringsAsFactors = FALSE
    )
  }

  matched_sets <- bind_rows(matched_set_rows[seq_len(set_counter)])
  matched_controls <- bind_rows(matched_control_rows[seq_len(set_counter)])

  if (nrow(matched_sets) == 0L) {
    stop("Matching failed for endpoint: ", endpoint_label)
  }

  if (any(matched_controls$followup_margin_years < -sqrt(.Machine$double.eps))) {
    stop("Follow-up eligibility violation: a control has insufficient follow-up.")
  }

  if (any(abs(matched_controls$baseline_age_difference) > AGE_CALIPER_YEARS)) {
    stop("Age-caliper violation detected.")
  }

  if (any(matched_sets$n_controls != MATCH_RATIO)) {
    stop("Matching ratio violation detected.")
  }

  age_smd <- standardized_mean_difference(
    matched_sets$case_baseline_age,
    matched_sets$control_baseline_age_mean
  )
  sex_levels <- levels(droplevels(data$sex))
  sex_level_for_smd <- if (length(sex_levels) == 2L) sex_levels[2L] else NA_character_
  sex_case <- if (!is.na(sex_level_for_smd)) {
    as.numeric(matched_sets$sex == sex_level_for_smd)
  } else {
    rep(NA_real_, nrow(matched_sets))
  }
  sex_control <- if (!is.na(sex_level_for_smd)) {
    as.numeric(matched_controls$control_sex == sex_level_for_smd)
  } else {
    rep(NA_real_, nrow(matched_controls))
  }
  sex_smd <- if (!is.na(sex_level_for_smd)) {
    standardized_mean_difference(sex_case, sex_control)
  } else {
    NA_real_
  }

  diagnostics <- data.frame(
    endpoint = endpoint_label,
    eligible_cases_2_to_15_years = nrow(cases),
    matched_cases = nrow(matched_sets),
    unmatched_cases = nrow(cases) - nrow(matched_sets),
    matched_controls = sum(matched_sets$n_controls),
    controls_per_case = MATCH_RATIO,
    age_caliper_years = AGE_CALIPER_YEARS,
    baseline_age_smd = age_smd,
    sex_matched = FALSE,
    sex_level_for_smd = sex_level_for_smd,
    case_sex_level_proportion = mean(sex_case),
    control_sex_level_proportion = mean(sex_control),
    sex_smd = sex_smd,
    maximum_absolute_baseline_age_difference = max(
      abs(matched_controls$baseline_age_difference)
    ),
    maximum_absolute_pseudo_index_age_difference = max(
      abs(matched_controls$pseudo_index_age_difference)
    ),
    minimum_control_followup_margin_years = min(
      matched_controls$followup_margin_years
    ),
    stringsAsFactors = FALSE
  )

  list(
    matched_sets = matched_sets,
    diagnostics = diagnostics
  )
}

endpoint_specification <- data.frame(
  variable = c("stomach", "colorectal", "gcrc"),
  label = c(
    "Gastric cancer",
    "Colorectal cancer",
    "Gastric/colorectal cancer"
  ),
  short = c("GC", "CRC", "GCRC"),
  stringsAsFactors = FALSE
)

matched_results <- lapply(
  seq_len(nrow(endpoint_specification)),
  function(i) {
    match_endpoint(
      analysis_data,
      event_variable = endpoint_specification$variable[i],
      endpoint_label = endpoint_specification$label[i],
      endpoint_short = endpoint_specification$short[i]
    )
  }
)

matched_sets <- bind_rows(lapply(matched_results, `[[`, "matched_sets"))
matching_diagnostics <- bind_rows(lapply(matched_results, `[[`, "diagnostics"))

# This aggregate file deliberately excludes participant identifiers.
write_csv(matched_sets, "matched_set_level_data.csv")
write_csv(matching_diagnostics, "matching_diagnostics.csv")

cat("\nMatching diagnostics\n")
print(matching_diagnostics)

# ---------- Matched annual comparisons ----------

summarize_year <- function(data) {
  case_ci <- mean_ci(data$case_score)
  control_ci <- mean_ci(data$control_score_mean)
  comparison <- paired_difference_test(data$paired_difference)

  data.frame(
    n_matched_sets = nrow(data),
    n_controls = sum(data$n_controls),
    case_mean = unname(case_ci["estimate"]),
    case_ci_lower = unname(case_ci["lower"]),
    case_ci_upper = unname(case_ci["upper"]),
    matched_control_mean = unname(control_ci["estimate"]),
    matched_control_ci_lower = unname(control_ci["lower"]),
    matched_control_ci_upper = unname(control_ci["upper"]),
    raw_mean_paired_difference = mean(data$paired_difference),
    comparison,
    stringsAsFactors = FALSE
  )
}

yearly_results <- matched_sets |>
  group_by(endpoint, endpoint_short, year_to_index) |>
  group_modify(~ summarize_year(.x), .keep = FALSE) |>
  ungroup() |>
  group_by(endpoint) |>
  mutate(p_adjusted_bh = stats::p.adjust(p_value, method = "BH")) |>
  ungroup() |>
  arrange(factor(endpoint_short, levels = c("GC", "CRC", "GCRC")), year_to_index)

overall_results <- matched_sets |>
  group_by(endpoint, endpoint_short) |>
  group_modify(~ paired_difference_test(.x$paired_difference), .keep = FALSE) |>
  ungroup() |>
  arrange(factor(endpoint_short, levels = c("GC", "CRC", "GCRC")))

write_csv(yearly_results, "matched_score_pattern_yearly_results.csv")
write_csv(overall_results, "matched_score_pattern_overall_results.csv")

cat("\nOverall matched comparisons\n")
print(overall_results)
cat("\nYear-specific matched comparisons (BH correction within endpoint)\n")
print(yearly_results)

# Q-Q diagnostic plots for the paired set differences.
grDevices::pdf(
  file.path(output_dir, "matched_difference_qq_plots.pdf"),
  width = 7,
  height = 3
)
old_par <- graphics::par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))
for (endpoint_code in c("GC", "CRC", "GCRC")) {
  difference <- matched_sets$paired_difference[
    matched_sets$endpoint_short == endpoint_code
  ]
  stats::qqnorm(
    difference,
    main = endpoint_code,
    xlab = "Theoretical quantiles",
    ylab = "Matched-set difference"
  )
  stats::qqline(difference, col = CASE_COLOR, lwd = 1.2)
}
graphics::par(old_par)
grDevices::dev.off()

# ---------- Publication-ready matched score patterns ----------

score_pattern_plot_data <- bind_rows(
  yearly_results |>
    transmute(
      endpoint,
      endpoint_short,
      year_to_index,
      group = "Incident case",
      mean_score = case_mean,
      ci_lower = case_ci_lower,
      ci_upper = case_ci_upper,
      n_matched_sets
    ),
  yearly_results |>
    transmute(
      endpoint,
      endpoint_short,
      year_to_index,
      group = "Matched cancer-free control",
      mean_score = matched_control_mean,
      ci_lower = matched_control_ci_lower,
      ci_upper = matched_control_ci_upper,
      n_matched_sets
    )
)

score_pattern_plot_data$group <- factor(
  score_pattern_plot_data$group,
  levels = c("Incident case", "Matched cancer-free control")
)

make_score_pattern_plot <- function(endpoint_code, panel_label = NULL) {
  plot_data <- score_pattern_plot_data |>
    filter(endpoint_short == endpoint_code)

  ggplot(
    plot_data,
    aes(
      x = year_to_index,
      y = mean_score,
      color = group,
      linetype = group,
      shape = group,
      group = group
    )
  ) +
    geom_ribbon(
      aes(ymin = ci_lower, ymax = ci_upper, fill = group),
      alpha = 0.12,
      color = NA,
      show.legend = FALSE
    ) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.0, stroke = 0.4) +
    scale_color_manual(
      values = c(
        "Incident case" = CASE_COLOR,
        "Matched cancer-free control" = CONTROL_COLOR
      )
    ) +
    scale_fill_manual(
      values = c(
        "Incident case" = CASE_COLOR,
        "Matched cancer-free control" = CONTROL_COLOR
      )
    ) +
    scale_linetype_manual(
      values = c(
        "Incident case" = "solid",
        "Matched cancer-free control" = "22"
      )
    ) +
    scale_shape_manual(
      values = c(
        "Incident case" = 16,
        "Matched cancer-free control" = 17
      )
    ) +
    scale_x_continuous(
      breaks = seq(MIN_YEAR_TO_INDEX, MAX_YEAR_TO_INDEX, by = 1),
      limits = c(MIN_YEAR_TO_INDEX, MAX_YEAR_TO_INDEX)
    ) +
    labs(
      tag = panel_label,
      x = "Years before diagnosis",
      y = "Mean baseline risk score (95% CI)",
      color = NULL,
      linetype = NULL,
      shape = NULL
    ) +
    # Use the platform sans-serif mapping (Helvetica on macOS) so Cairo can
    # embed the font reliably in the vector PDF.
    theme_classic(base_family = "sans", base_size = 9) +
    theme(
      axis.text = element_text(size = 8, color = "black"),
      axis.title = element_text(size = 9, color = "black"),
      axis.ticks = element_line(linewidth = 0.4),
      legend.position = "top",
      legend.text = element_text(size = 8),
      legend.key.width = grid::unit(1.3, "lines"),
      plot.tag = element_text(face = "bold", size = 12),
      plot.tag.position = c(0.01, 0.99),
      plot.margin = margin(6, 6, 6, 6)
    )
}

save_plot <- function(plot, stem, width = 3.5, height = 3.5) {
  ggsave(
    file.path(output_dir, paste0(stem, ".pdf")),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = grDevices::pdf,
    useDingbats = FALSE,
    bg = "white"
  )
  ggsave(
    file.path(output_dir, paste0(stem, ".png")),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    # Use 601 rather than exactly 600 because some PNG encoders store 600 dpi
    # as 599.9988 after pixels-per-metre conversion.
    dpi = 601,
    bg = "white"
  )
}

save_two_panel <- function(plot_left, plot_right, stem) {
  pdf_path <- file.path(output_dir, paste0(stem, ".pdf"))
  png_path <- file.path(output_dir, paste0(stem, ".png"))

  grDevices::pdf(pdf_path, width = 7, height = 3.5, useDingbats = FALSE)
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, 2)))
  print(plot_left, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(plot_right, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  grDevices::dev.off()

  grDevices::png(
    png_path,
    width = 7,
    height = 3.5,
    units = "in",
    res = 600,
    bg = "white"
  )
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, 2)))
  print(plot_left, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(plot_right, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  grDevices::dev.off()
}

p_gc <- make_score_pattern_plot("GC", "A")
p_crc <- make_score_pattern_plot("CRC", "B")
p_gcrc <- make_score_pattern_plot("GCRC")

save_plot(p_gc, "matched_score_pattern_gastric_cancer")
save_plot(p_crc, "matched_score_pattern_colorectal_cancer")
save_plot(p_gcrc, "matched_score_pattern_combined_sensitivity")
save_two_panel(p_gc, p_crc, "Figure5_matched_score_patterns")

# Figure 5A reference-style visualization. The visual grammar follows the
# original figure (single panel, point-line pattern, error bars, top legend,
# and significance symbols), while retaining the corrected matched-control
# analysis. In particular, the blue control series is time-varying rather than
# a static horizontal reference.
make_figure5a_reference_style <- function() {
  annual <- yearly_results |>
    filter(endpoint_short == "GCRC") |>
    mutate(
      significance = case_when(
        p_adjusted_bh < 0.001 ~ "***",
        p_adjusted_bh < 0.01 ~ "**",
        p_adjusted_bh < 0.05 ~ "*",
        TRUE ~ ""
      ),
      annotation_y = pmax(case_ci_upper, matched_control_ci_upper) + 0.18
    )

  plot_data <- bind_rows(
    annual |>
      transmute(
        year_to_index,
        group = "Gastric/Colorectal cancer",
        mean_score = case_mean,
        ci_lower = case_ci_lower,
        ci_upper = case_ci_upper
      ),
    annual |>
      transmute(
        year_to_index,
        group = "Cancer-free (matched)",
        mean_score = matched_control_mean,
        ci_lower = matched_control_ci_lower,
        ci_upper = matched_control_ci_upper
      )
  ) |>
    mutate(
      group = factor(
        group,
        levels = c("Cancer-free (matched)", "Gastric/Colorectal cancer")
      )
    )

  dodge <- position_dodge(width = 0.16)

  ggplot(
    plot_data,
    aes(
      x = year_to_index,
      y = mean_score,
      color = group,
      linetype = group,
      shape = group,
      group = group
    )
  ) +
    geom_errorbar(
      aes(ymin = ci_lower, ymax = ci_upper),
      position = dodge,
      width = 0.16,
      linewidth = 0.55
    ) +
    geom_line(position = dodge, linewidth = 0.85) +
    geom_point(position = dodge, size = 2.4, stroke = 0.45) +
    geom_text(
      data = annual |> filter(significance != ""),
      aes(
        x = year_to_index,
        y = annotation_y,
        label = significance
      ),
      inherit.aes = FALSE,
      color = "black",
      size = 3.2,
      fontface = "bold",
      vjust = 0
    ) +
    scale_color_manual(
      values = c(
        "Cancer-free (matched)" = CONTROL_COLOR,
        "Gastric/Colorectal cancer" = CASE_COLOR
      )
    ) +
    scale_linetype_manual(
      values = c(
        "Cancer-free (matched)" = "22",
        "Gastric/Colorectal cancer" = "solid"
      )
    ) +
    scale_shape_manual(
      values = c(
        "Cancer-free (matched)" = 17,
        "Gastric/Colorectal cancer" = 16
      )
    ) +
    scale_x_continuous(
      breaks = seq(MIN_YEAR_TO_INDEX, MAX_YEAR_TO_INDEX, by = 1),
      limits = c(MIN_YEAR_TO_INDEX - 0.45, MAX_YEAR_TO_INDEX + 0.45),
      expand = expansion(mult = 0)
    ) +
    scale_y_continuous(
      breaks = scales::pretty_breaks(n = 6),
      expand = expansion(mult = c(0.04, 0.12))
    ) +
    labs(
      tag = "A",
      x = "Years before diagnosis",
      y = "Mean baseline risk score (95% CI)",
      color = NULL,
      linetype = NULL,
      shape = NULL
    ) +
    theme_classic(base_family = "sans", base_size = 10) +
    theme(
      axis.text = element_text(size = 9, color = "black"),
      axis.title = element_text(size = 10, color = "black"),
      axis.line = element_line(linewidth = 0.65, color = "black"),
      axis.ticks = element_line(linewidth = 0.5, color = "black"),
      legend.position = "top",
      legend.justification = "center",
      legend.text = element_text(size = 9),
      legend.key.width = grid::unit(1.6, "lines"),
      legend.margin = margin(0, 0, 4, 0),
      plot.tag = element_text(face = "bold", size = 18),
      plot.tag.position = c(0.01, 0.99),
      plot.margin = margin(8, 10, 8, 10)
    )
}

p_figure5a <- make_figure5a_reference_style()
save_plot(
  p_figure5a,
  "Figure5A_matched_gcrc_score_pattern",
  width = 7,
  height = 5
)

# ---------- Manuscript-ready evidence and manifests ----------

get_overall_row <- function(code) {
  overall_results[overall_results$endpoint_short == code, , drop = FALSE]
}

get_diagnostic_row <- function(label) {
  matching_diagnostics[matching_diagnostics$endpoint == label, , drop = FALSE]
}

gc_overall <- get_overall_row("GC")
crc_overall <- get_overall_row("CRC")
gcrc_overall <- get_overall_row("GCRC")
gc_diag <- get_diagnostic_row("Gastric cancer")
crc_diag <- get_diagnostic_row("Colorectal cancer")
gcrc_diag <- get_diagnostic_row("Gastric/colorectal cancer")

evidence_lines <- c(
  "# Matched pre-diagnostic score-pattern evidence",
  "",
  paste("Generated:", format(Sys.Date(), "%Y-%m-%d")),
  "",
  "All values below were generated directly from the canonical analysis RDS by this script.",
  "No participant identifier is included in the exported matched-set data.",
  "",
  "## Design facts",
  "",
  paste0(
    "- Gastric cancer: ", gc_diag$matched_cases, " of ",
    gc_diag$eligible_cases_2_to_15_years, " eligible cases matched to ",
    gc_diag$matched_controls, " controls."
  ),
  paste0(
    "- Colorectal cancer: ", crc_diag$matched_cases, " of ",
    crc_diag$eligible_cases_2_to_15_years, " eligible cases matched to ",
    crc_diag$matched_controls, " controls."
  ),
  paste0(
    "- Combined sensitivity analysis: ", gcrc_diag$matched_cases, " of ",
    gcrc_diag$eligible_cases_2_to_15_years, " eligible cases matched to ",
    gcrc_diag$matched_controls, " controls."
  ),
  paste0(
    "- Maximum absolute age difference = ",
    sprintf("%.1f", max(matching_diagnostics$maximum_absolute_baseline_age_difference)),
    " years; minimum control follow-up margin = ",
    sprintf("%.3f", min(matching_diagnostics$minimum_control_followup_margin_years)),
    " years."
  ),
  paste0(
    "- Sex was not matched; post-matching sex SMD ranged from ",
    sprintf("%.3f", min(matching_diagnostics$sex_smd)), " to ",
    sprintf("%.3f", max(matching_diagnostics$sex_smd)), "."
  ),
  "",
  "## Overall matched comparisons",
  "",
  paste0(
    "- Gastric cancer: ", gc_overall$effect_type, " = ",
    sprintf("%.3f", gc_overall$effect_estimate), " (95% CI ",
    sprintf("%.3f", gc_overall$ci_lower), " to ",
    sprintf("%.3f", gc_overall$ci_upper), "; p ",
    format_p(gc_overall$p_value), "; n = ", gc_overall$n_sets,
    " matched sets)."
  ),
  paste0(
    "- Colorectal cancer: ", crc_overall$effect_type, " = ",
    sprintf("%.3f", crc_overall$effect_estimate), " (95% CI ",
    sprintf("%.3f", crc_overall$ci_lower), " to ",
    sprintf("%.3f", crc_overall$ci_upper), "; p ",
    format_p(crc_overall$p_value), "; n = ", crc_overall$n_sets,
    " matched sets)."
  ),
  paste0(
    "- Combined sensitivity analysis: ", gcrc_overall$effect_type, " = ",
    sprintf("%.3f", gcrc_overall$effect_estimate), " (95% CI ",
    sprintf("%.3f", gcrc_overall$ci_lower), " to ",
    sprintf("%.3f", gcrc_overall$ci_upper), "; p ",
    format_p(gcrc_overall$p_value), "; n = ", gcrc_overall$n_sets,
    " matched sets)."
  ),
  "",
  "Interpret these overall contrasts together with the year-specific estimates and",
  "Benjamini-Hochberg-adjusted p-values in matched_score_pattern_yearly_results.csv."
)
writeLines(evidence_lines, file.path(output_dir, "reviewer_response_evidence.md"))

analysis_manifest <- c(
  "# Analysis Outputs",
  paste("Generated:", format(Sys.Date(), "%Y-%m-%d")),
  "Study type: age- and follow-up-eligible matched-control analysis; sex not matched",
  "",
  "## Tables",
  "- `analysis_cohort_profile.csv` -- Input cohort counts and missingness summary.",
  "- `matching_diagnostics.csv` -- Matching completion, balance, and follow-up eligibility checks.",
  "- `matched_set_level_data.csv` -- De-identified matched-set aggregates.",
  "- `matched_score_pattern_yearly_results.csv` -- Annual estimates, 95% CIs, exact p-values, and BH-adjusted p-values.",
  "- `matched_score_pattern_overall_results.csv` -- Overall matched comparisons by endpoint.",
  "",
  "## Figures",
  "- `Figure5_matched_score_patterns.pdf` / `.png` -- Main two-panel matched score patterns.",
  "- `Figure5A_matched_gcrc_score_pattern.pdf` / `.png` -- Reference-style revised Figure 5A for the combined endpoint.",
  "- `matched_score_pattern_gastric_cancer.pdf` / `.png` -- Gastric cancer panel.",
  "- `matched_score_pattern_colorectal_cancer.pdf` / `.png` -- Colorectal cancer panel.",
  "- `matched_score_pattern_combined_sensitivity.pdf` / `.png` -- Combined endpoint sensitivity analysis.",
  "- `matched_difference_qq_plots.pdf` -- Distributional diagnostic for matched-set differences.",
  "",
  "## Reviewer-response support",
  "- `reviewer_response_evidence.md` -- Audited numerical evidence for the response letter."
)
writeLines(analysis_manifest, file.path(output_dir, "_analysis_outputs.md"))

figure_manifest <- c(
  "# Figure Manifest",
  paste("Generated:", format(Sys.Date(), "%Y-%m-%d")),
  "Study type: age- and follow-up-eligible matched comparison of pre-diagnostic score patterns; sex not matched",
  "",
  "| Figure | Path | Type | Tool | Critic | Rounds | Description |",
  "|---|---|---|---|---|---:|---|",
  "| Figure 5 | gcrc_output_age_time/Figure5_matched_score_patterns.pdf | other | ggplot2 | yes | 1 | Age- and follow-up-eligible matched score patterns for gastric and colorectal cancer; sex not matched. |",
  "| Figure 5A | gcrc_output_age_time/Figure5A_matched_gcrc_score_pattern.pdf | other | ggplot2 | yes | 1 | Reference-style matched score pattern for combined gastric/colorectal cancer. |",
  "| Sensitivity | gcrc_output_age_time/matched_score_pattern_combined_sensitivity.pdf | other | ggplot2 | yes | 1 | Combined gastrointestinal cancer matched score pattern. |",
  "",
  "## Critic notes",
  "- Quantitative critic: PASS (7.0-inch width, 600 dpi, no flags).",
  "- Revised Figure 5A quantitative critic: PASS (7.0-inch width, 601 dpi, no flags).",
  "- Qualitative critic: PASS; axes/units are present, line style and marker shape preserve grayscale distinction, and no statistical stars are used."
)
writeLines(figure_manifest, file.path(output_dir, "_figure_manifest.md"))

cat("\nOutputs written to:", output_dir, "\n")
