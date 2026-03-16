# 01_main_logistic_analysis.R
# Primary multivariable model for internal vs external outflow procedures.

source("R/00_helpers.R")

input_file <- "data/glaucoma_surgery.xlsx"   # replace with your local authorized dataset
output_dir <- "results/main"
make_output_dir(output_dir)

cohort <- build_analytic_cohort(input_file)
df_analysis <- cohort$df_analysis
df_model_cc <- cohort$df_model_cc

fit_main <- glm(
  Surgery_bin ~ Age10 + Sex + GlaucomaType + Year_c,
  data = df_model_cc,
  family = binomial
)

or_table <- broom::tidy(fit_main, conf.int = TRUE, exponentiate = TRUE) %>%
  dplyr::mutate(
    term_pretty = dplyr::case_when(
      term == "Age10" ~ "Age (per 10 years)",
      term == "SexFemale" ~ "Female (vs Male)",
      term == "GlaucomaTypePOAG" ~ "POAG (vs PACG)",
      term == "GlaucomaTypeSG" ~ "SG (vs PACG)",
      term == "GlaucomaTypePG" ~ "PG (vs PACG)",
      term == "GlaucomaTypeOthers" ~ "Others (vs PACG)",
      term == "Year_c" ~ "Calendar year (per 1 year)",
      TRUE ~ term
    )
  )

# Discrimination and calibration
y <- df_model_cc$Surgery_bin
pred_prob <- stats::predict(fit_main, type = "response")
auc_obj <- pROC::roc(response = y, predictor = pred_prob)
auc_ci <- pROC::ci.auc(auc_obj)
brier <- mean((pred_prob - y)^2)

hl <- ResourceSelection::hoslem.test(x = y, y = pred_prob, g = 10)

p_clip <- pmin(pmax(pred_prob, 1e-6), 1 - 1e-6)
logit_pred <- qlogis(p_clip)

cal_intercept_fit <- glm(y ~ 1, family = binomial, offset = logit_pred)
cal_slope_fit <- glm(y ~ logit_pred, family = binomial)

calibration_summary <- data.frame(
  metric = c("AUC", "AUC_low", "AUC_high", "Brier", "HL_chisq", "HL_df", "HL_p",
             "Calibration_intercept", "Calibration_slope"),
  value = c(as.numeric(auc_obj$auc), auc_ci[1], auc_ci[3], brier,
            as.numeric(hl$statistic), as.numeric(hl$parameter), as.numeric(hl$p.value),
            unname(coef(cal_intercept_fit)[1]), unname(coef(cal_slope_fit)["logit_pred"]))
)

# Forest plot
forest_data <- or_table %>%
  dplyr::filter(term != "(Intercept)") %>%
  dplyr::mutate(
    sig = ifelse(p.value < 0.05, "Significant", "Not significant"),
    term_pretty = factor(
      term_pretty,
      levels = rev(c(
        "Age (per 10 years)",
        "Female (vs Male)",
        "POAG (vs PACG)",
        "SG (vs PACG)",
        "PG (vs PACG)",
        "Others (vs PACG)",
        "Calendar year (per 1 year)"
      ))
    )
  )

p_forest <- ggplot2::ggplot(
  forest_data,
  ggplot2::aes(y = term_pretty, x = estimate, xmin = conf.low, xmax = conf.high, shape = sig)
) +
  ggplot2::geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.6, colour = "grey40") +
  ggplot2::geom_errorbarh(height = 0.15, size = 0.6, colour = "black") +
  ggplot2::geom_point(size = 2.7, stroke = 0.8, colour = "black") +
  ggplot2::scale_x_log10(breaks = c(0.5, 1, 2, 4), expand = c(0.02, 0.02)) +
  ggplot2::scale_shape_manual(values = c("Significant" = 16, "Not significant" = 1)) +
  ggplot2::labs(
    x = "Adjusted odds ratio (log scale)",
    y = "",
    title = "Factors associated with choosing internal vs external outflow surgery"
  ) +
  ggplot2::theme_classic(base_size = 13)

# Calibration plot
cal_df <- data.frame(y = y, pred = p_clip) %>%
  dplyr::mutate(decile = dplyr::ntile(pred, 10))

cal_summary <- cal_df %>%
  dplyr::group_by(decile) %>%
  dplyr::summarise(
    n = dplyr::n(),
    pred_mean = mean(pred),
    obs_rate = mean(y),
    events = sum(y),
    .groups = "drop"
  )

ci <- binom::binom.wilson(cal_summary$events, cal_summary$n, conf.level = 0.95)
cal_summary$obs_low <- ci$lower
cal_summary$obs_high <- ci$upper

p_cal <- ggplot2::ggplot(cal_summary, ggplot2::aes(x = pred_mean, y = obs_rate)) +
  ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = obs_low, ymax = obs_high), width = 0.00) +
  ggplot2::geom_point(size = 2) +
  ggplot2::scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  ggplot2::labs(
    x = "Predicted probability (internal outflow)",
    y = "Observed proportion (internal outflow)",
    title = "Calibration plot for the primary multivariable model"
  ) +
  ggplot2::theme_classic(base_size = 12)

writexl::write_xlsx(
  list(
    cohort_counts = data.frame(
      step = c("Eligible procedures", "Internal/external subset", "Complete-case sample"),
      n = c(nrow(df_analysis), nrow(cohort$df_model), nrow(df_model_cc))
    ),
    odds_ratios = or_table,
    calibration_summary = calibration_summary,
    calibration_deciles = cal_summary
  ),
  path = file.path(output_dir, "main_logistic_outputs.xlsx")
)

ggplot2::ggsave(file.path(output_dir, "forest_plot.pdf"), p_forest, width = 6, height = 5, units = "in")
ggplot2::ggsave(file.path(output_dir, "forest_plot.tiff"), p_forest, width = 6, height = 5, units = "in",
                dpi = 600, compression = "lzw")
ggplot2::ggsave(file.path(output_dir, "calibration_plot.pdf"), p_cal, width = 6, height = 5, units = "in")
ggplot2::ggsave(file.path(output_dir, "calibration_plot.tiff"), p_cal, width = 6, height = 5, units = "in",
                dpi = 600, compression = "lzw")

message("Main analysis complete. Outputs saved in: ", output_dir)
