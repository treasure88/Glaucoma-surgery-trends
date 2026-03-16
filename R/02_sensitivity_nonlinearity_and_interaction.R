# 02_sensitivity_nonlinearity_and_interaction.R
# Sensitivity analyses for linearity assumptions and subtype-specific temporal trends.

source("R/00_helpers.R")

input_file <- "data/glaucoma_surgery.xlsx"   # replace with your local authorized dataset
output_dir <- "results/sensitivity_nonlinearity"
make_output_dir(output_dir)

cohort <- build_analytic_cohort(input_file)
df_model_cc <- cohort$df_model_cc

# Linearity checks
fit_linear <- glm(
  Surgery_bin ~ Age10 + Sex + GlaucomaType + Year_c,
  data = df_model_cc,
  family = binomial
)

fit_age_spline <- glm(
  Surgery_bin ~ splines::ns(Age10, df = 4) + Sex + GlaucomaType + Year_c,
  data = df_model_cc,
  family = binomial
)

fit_year_spline <- glm(
  Surgery_bin ~ Age10 + Sex + GlaucomaType + splines::ns(Year_c, df = 4),
  data = df_model_cc,
  family = binomial
)

fit_piecewise_year <- glm(
  Surgery_bin ~ Age10 + Sex + GlaucomaType + Year_c + pmax(Year_model - 2020, 0),
  data = df_model_cc,
  family = binomial
)

fit_year_factor <- glm(
  Surgery_bin ~ Age10 + Sex + GlaucomaType + Year_f,
  data = df_model_cc,
  family = binomial
)

lr_table <- data.frame(
  comparison = c(
    "Age linear vs age spline",
    "Year linear vs year spline",
    "Year linear vs piecewise year",
    "Year linear vs year factor"
  ),
  p_value = c(
    anova(fit_linear, fit_age_spline, test = "LRT")[2, "Pr(>Chi)"],
    anova(fit_linear, fit_year_spline, test = "LRT")[2, "Pr(>Chi)"],
    anova(fit_linear, fit_piecewise_year, test = "LRT")[2, "Pr(>Chi)"],
    anova(fit_linear, fit_year_factor, test = "LRT")[2, "Pr(>Chi)"]
  )
)

# Subtype by year interaction
fit_interaction <- glm(
  Surgery_bin ~ Age10 + Sex + GlaucomaType * Year_c,
  data = df_model_cc,
  family = binomial
)

coef_table <- broom::tidy(fit_interaction, conf.int = TRUE, exponentiate = TRUE)

extract_subtype_annual_or <- function(model) {
  coefs <- stats::coef(model)
  vc <- stats::vcov(model)
  types <- levels(df_model_cc$GlaucomaType)

  rows <- lapply(types, function(gt) {
    if (gt == "PACG") {
      beta <- coefs["Year_c"]
      se <- sqrt(vc["Year_c", "Year_c"])
    } else {
      interaction_term <- paste0("GlaucomaType", gt, ":Year_c")
      if (!interaction_term %in% names(coefs)) {
        interaction_term <- paste0("Year_c:GlaucomaType", gt)
      }
      beta <- coefs["Year_c"] + coefs[interaction_term]
      se <- sqrt(
        vc["Year_c", "Year_c"] +
          vc[interaction_term, interaction_term] +
          2 * vc["Year_c", interaction_term]
      )
    }
    data.frame(
      GlaucomaType = gt,
      annual_or = exp(beta),
      lcl = exp(beta - 1.96 * se),
      ucl = exp(beta + 1.96 * se)
    )
  })

  do.call(rbind, rows)
}

interaction_summary <- extract_subtype_annual_or(fit_interaction)

# Predicted probabilities by subtype and calendar year
pred_grid <- tidyr::expand_grid(
  Year_model = YEAR_MIN:YEAR_MAX,
  GlaucomaType = levels(df_model_cc$GlaucomaType)
) %>%
  dplyr::mutate(
    Year_c = Year_model - YEAR_MIN,
    Age10 = median(df_model_cc$Age10, na.rm = TRUE),
    Sex = factor("Male", levels = levels(df_model_cc$Sex)),
    GlaucomaType = factor(GlaucomaType, levels = levels(df_model_cc$GlaucomaType))
  )

pr <- predict(fit_interaction, newdata = pred_grid, type = "link", se.fit = TRUE)
pred_grid <- pred_grid %>%
  dplyr::mutate(
    fit_link = pr$fit,
    se_link = pr$se.fit,
    p = plogis(fit_link),
    lo = plogis(fit_link - 1.96 * se_link),
    hi = plogis(fit_link + 1.96 * se_link)
  )

p_traj <- ggplot2::ggplot(pred_grid, ggplot2::aes(x = Year_model, y = p, group = GlaucomaType)) +
  ggplot2::geom_line(linewidth = 0.9) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  ggplot2::facet_wrap(~ GlaucomaType, ncol = 3) +
  ggplot2::scale_x_continuous(breaks = YEAR_MIN:YEAR_MAX) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  ggplot2::labs(
    x = "Calendar year",
    y = "Adjusted probability of internal outflow surgery",
    title = "Subtype-specific adjusted temporal trajectories"
  ) +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

writexl::write_xlsx(
  list(
    likelihood_ratio_tests = lr_table,
    interaction_coefficients = coef_table,
    subtype_specific_annual_or = interaction_summary,
    predicted_probabilities = pred_grid
  ),
  path = file.path(output_dir, "sensitivity_nonlinearity_and_interaction.xlsx")
)

ggplot2::ggsave(file.path(output_dir, "subtype_trajectories.pdf"), p_traj, width = 12, height = 7, units = "in")
ggplot2::ggsave(file.path(output_dir, "subtype_trajectories.tiff"), p_traj, width = 12, height = 7, units = "in",
                dpi = 600, compression = "lzw")

message("Sensitivity analyses complete. Outputs saved in: ", output_dir)
