# 03_sensitivity_excluding_lens_only.R
# Trend analysis excluding lens-only procedures (P and P+I).

source("R/00_helpers.R")

input_file <- "data/glaucoma_surgery.xlsx"   # replace with your local authorized dataset
output_dir <- "results/sensitivity_excluding_lens_only"
make_output_dir(output_dir)

cohort <- build_analytic_cohort(input_file)
df_model_cc <- cohort$df_model_cc

lens_only_codes <- c("P", "P+I")

df_trend_base <- df_model_cc %>%
  dplyr::filter(SurgeryClass %in% c("Internal", "External")) %>%
  dplyr::filter(!is.na(Year_model), Year_model >= YEAR_MIN, Year_model <= YEAR_MAX)

trend_main <- df_trend_base %>%
  dplyr::group_by(Year_model) %>%
  dplyr::summarise(
    n = dplyr::n(),
    internal_n = sum(SurgeryClass == "Internal"),
    internal_pct = internal_n / n,
    .groups = "drop"
  ) %>%
  dplyr::arrange(Year_model)

trend_sensitivity <- df_trend_base %>%
  dplyr::filter(!(SurgeryType_std %in% lens_only_codes)) %>%
  dplyr::group_by(Year_model) %>%
  dplyr::summarise(
    n = dplyr::n(),
    internal_n = sum(SurgeryClass == "Internal"),
    internal_pct = internal_n / n,
    .groups = "drop"
  ) %>%
  dplyr::arrange(Year_model)

trend_test <- function(dat, label) {
  dat <- dat %>% dplyr::mutate(Surgery_bin = dplyr::if_else(SurgeryClass == "Internal", 1L, 0L))
  fit_u <- glm(Surgery_bin ~ Year_c, data = dat, family = binomial)
  fit_a <- glm(Surgery_bin ~ Age10 + Sex + GlaucomaType + Year_c, data = dat, family = binomial)

  dplyr::bind_rows(
    broom::tidy(fit_u, conf.int = TRUE, exponentiate = TRUE) %>%
      dplyr::filter(term == "Year_c") %>%
      dplyr::transmute(analysis = label, model = "Unadjusted", OR = estimate,
                       CI_low = conf.low, CI_high = conf.high, p = p.value),
    broom::tidy(fit_a, conf.int = TRUE, exponentiate = TRUE) %>%
      dplyr::filter(term == "Year_c") %>%
      dplyr::transmute(analysis = label, model = "Adjusted", OR = estimate,
                       CI_low = conf.low, CI_high = conf.high, p = p.value)
  )
}

trend_test_table <- dplyr::bind_rows(
  trend_test(df_trend_base, "Main (including lens-only)"),
  trend_test(df_trend_base %>% dplyr::filter(!(SurgeryType_std %in% lens_only_codes)),
             "Sensitivity (excluding lens-only)")
)

trend_compare <- dplyr::bind_rows(
  trend_main %>% dplyr::mutate(analysis = "Main (including P/P+I)"),
  trend_sensitivity %>% dplyr::mutate(analysis = "Sensitivity (excluding P/P+I)")
)

p_trend <- ggplot2::ggplot(trend_compare, ggplot2::aes(x = Year_model, y = internal_pct, linetype = analysis)) +
  ggplot2::geom_line(linewidth = 0.9) +
  ggplot2::geom_point(size = 2) +
  ggplot2::scale_x_continuous(breaks = YEAR_MIN:YEAR_MAX) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  ggplot2::labs(
    x = "Year",
    y = "Proportion of internal approach",
    linetype = ""
  ) +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(legend.position = "top")

writexl::write_xlsx(
  list(
    trend_main = trend_main,
    trend_sensitivity = trend_sensitivity,
    year_effect_models = trend_test_table
  ),
  path = file.path(output_dir, "sensitivity_excluding_lens_only.xlsx")
)

ggplot2::ggsave(file.path(output_dir, "trend_overlay.pdf"), p_trend, width = 7, height = 4, units = "in")
ggplot2::ggsave(file.path(output_dir, "trend_overlay.tiff"), p_trend, width = 7, height = 4, units = "in",
                dpi = 600, compression = "lzw")

message("Lens-only sensitivity analysis complete. Outputs saved in: ", output_dir)
