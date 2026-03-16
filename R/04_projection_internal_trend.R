# 04_projection_internal_trend.R
# Projection of adjusted internal outflow probability through 2030.
# This script uses patient-level bootstrap sampling with replacement.

source("R/00_helpers.R")

input_file <- "data/glaucoma_surgery.xlsx"   # replace with your local authorized dataset
output_dir <- "results/projection"
make_output_dir(output_dir)

cohort <- build_analytic_cohort(input_file)
dat <- cohort$df_model_cc

if (!("patient_id" %in% names(dat))) {
  if ("患者唯一号" %in% names(dat)) {
    dat <- dat %>% dplyr::rename(patient_id = 患者唯一号)
  } else {
    stop("A patient ID column is required.")
  }
}

dat <- dat %>%
  dplyr::mutate(
    patient_id = as.character(patient_id),
    Surgery_bin = as.integer(Surgery_bin),
    Age10 = suppressWarnings(as.numeric(Age10)),
    Sex = factor(as.character(Sex), levels = c("Male", "Female")),
    GlaucomaType = factor(as.character(GlaucomaType), levels = c("PACG", "POAG", "SG", "PG", "Others")),
    Year_c = as.integer(Year_c)
  ) %>%
  dplyr::filter(
    !is.na(patient_id), !is.na(Surgery_bin), !is.na(Age10),
    !is.na(Sex), !is.na(GlaucomaType), !is.na(Year_c)
  )

year_all <- 2015:2030
B <- 1000

predict_standardized <- function(fit, dat, years) {
  res <- lapply(levels(dat$GlaucomaType), function(gt) {
    dat_gt <- dat %>% dplyr::filter(GlaucomaType == gt)
    sapply(years, function(y) {
      dat_tmp <- dat_gt
      dat_tmp$Year_c <- as.integer(y - YEAR_MIN)
      mean(stats::predict(fit, newdata = dat_tmp, type = "response"), na.rm = TRUE)
    })
  })
  res <- do.call(rbind, res)
  colnames(res) <- years

  tidyr::expand_grid(
    GlaucomaType = levels(dat$GlaucomaType),
    Year = years
  ) %>%
    dplyr::mutate(prob = as.vector(t(res)))
}

fit0 <- glm(
  Surgery_bin ~ Age10 + Sex + GlaucomaType + Year_c,
  data = dat,
  family = binomial
)

point <- predict_standardized(fit0, dat, year_all) %>%
  dplyr::rename(prob_hat = prob)

ids <- unique(dat$patient_id)
n_ids <- length(ids)

boot_mat <- matrix(NA_real_, nrow = B, ncol = nrow(point))

for (b in seq_len(B)) {
  sampled_ids <- sample(ids, n_ids, replace = TRUE)
  dat_b <- do.call(
    rbind,
    lapply(sampled_ids, function(id) dat[dat$patient_id == id, , drop = FALSE])
  )

  fit_b <- glm(
    Surgery_bin ~ Age10 + Sex + GlaucomaType + Year_c,
    data = dat_b,
    family = binomial
  )

  boot_mat[b, ] <- predict_standardized(fit_b, dat_b, year_all)$prob
}

ci_df <- t(apply(boot_mat, 2, function(x) stats::quantile(x, probs = c(0.025, 0.975), na.rm = TRUE)))
colnames(ci_df) <- c("lcl", "ucl")

pred_df <- cbind(point, as.data.frame(ci_df)) %>%
  dplyr::mutate(
    prob = prob_hat,
    Period = dplyr::if_else(Year <= YEAR_MAX, "Observed", "Forecast"),
    GlaucomaType = factor(GlaucomaType, levels = c("PACG", "POAG", "SG", "PG", "Others"))
  )

palette_okabe <- c(
  "PACG" = "#0072B2",
  "POAG" = "#009E73",
  "SG" = "#E69F00",
  "PG" = "#CC79A7",
  "Others" = "#7F7F7F"
)

label_df <- pred_df %>%
  dplyr::filter(Year == 2030) %>%
  dplyr::arrange(prob)

p_trend <- ggplot2::ggplot(
  pred_df,
  ggplot2::aes(x = Year, y = prob, colour = GlaucomaType, fill = GlaucomaType)
) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = lcl, ymax = ucl, group = GlaucomaType), alpha = 0.12, colour = NA) +
  ggplot2::geom_line(
    data = pred_df %>% dplyr::filter(Year <= YEAR_MAX),
    ggplot2::aes(group = GlaucomaType),
    linewidth = 1.0
  ) +
  ggplot2::geom_line(
    data = pred_df %>% dplyr::filter(Year > YEAR_MAX),
    ggplot2::aes(group = GlaucomaType),
    linewidth = 1.0,
    linetype = "dashed"
  ) +
  ggplot2::geom_vline(xintercept = YEAR_MAX + 0.5, linetype = "dotted", linewidth = 0.7, colour = "grey40") +
  ggplot2::geom_text(
    data = label_df,
    ggplot2::aes(label = as.character(GlaucomaType), x = 2030.2, y = prob),
    inherit.aes = FALSE,
    size = 3.2,
    colour = "grey15",
    hjust = 0
  ) +
  ggplot2::scale_colour_manual(values = palette_okabe, drop = FALSE) +
  ggplot2::scale_fill_manual(values = palette_okabe, drop = FALSE) +
  ggplot2::scale_x_continuous(breaks = seq(2015, 2030, by = 2), expand = ggplot2::expansion(mult = c(0.01, 0.08))) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  ggplot2::labs(
    title = "Observed and projected trends in internal outflow probability",
    subtitle = "Shaded bands indicate 95% confidence intervals from patient-level bootstrap sampling.",
    x = "Year",
    y = "Adjusted probability of internal outflow"
  ) +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(
    legend.position = "bottom",
    plot.margin = ggplot2::margin(8, 24, 8, 8)
  ) +
  ggplot2::coord_cartesian(xlim = c(2015, 2031), clip = "off") +
  ggplot2::guides(fill = "none")

write.csv(pred_df, file.path(output_dir, "projected_internal_trend.csv"), row.names = FALSE)
ggplot2::ggsave(file.path(output_dir, "projected_internal_trend.pdf"), p_trend, width = 9, height = 6, units = "in")
ggplot2::ggsave(file.path(output_dir, "projected_internal_trend.tiff"), p_trend, width = 9, height = 6, units = "in",
                dpi = 600, compression = "lzw")

message("Projection analysis complete. Outputs saved in: ", output_dir)
