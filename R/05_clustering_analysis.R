# 05_clustering_analysis.R
# Patient-level clinical profile clustering using Gower distance and PAM.

source("R/00_helpers.R")

input_file <- "data/glaucoma_surgery.xlsx"   # replace with your local authorized dataset
output_dir <- "results/clustering"
make_output_dir(output_dir)

cohort <- build_analytic_cohort(input_file)
dat_proc <- cohort$df_model_cc %>% dplyr::mutate(.row_id = dplyr::row_number())

if (!("patient_id" %in% names(dat_proc))) {
  if ("患者唯一号" %in% names(dat_proc)) {
    dat_proc <- dat_proc %>% dplyr::rename(patient_id = 患者唯一号)
  } else {
    stop("A patient ID column is required.")
  }
}

dat_proc <- dat_proc %>%
  dplyr::mutate(patient_id = as.character(patient_id))

if (!("Age" %in% names(dat_proc))) {
  dat_proc <- dat_proc %>% dplyr::mutate(Age = Age10 * 10)
}

has_reop <- "Reop" %in% names(dat_proc)

if (has_reop) {
  dat_proc <- dat_proc %>%
    dplyr::mutate(
      Reop_num = dplyr::case_when(
        is.na(Reop) ~ NA_real_,
        is.numeric(Reop) ~ as.numeric(Reop),
        is.integer(Reop) ~ as.numeric(Reop),
        is.factor(Reop) ~ as.numeric(as.character(Reop)),
        is.character(Reop) & Reop %in% c("是", "Yes", "YES", "1") ~ 1,
        is.character(Reop) & Reop %in% c("否", "No", "NO", "0") ~ 0,
        TRUE ~ suppressWarnings(as.numeric(Reop))
      )
    )
} else {
  dat_proc$Reop_num <- NA_real_
}

dat_pat <- dat_proc %>%
  dplyr::filter(!is.na(patient_id)) %>%
  dplyr::group_by(patient_id) %>%
  dplyr::arrange(Year_model, .by_group = TRUE) %>%
  dplyr::summarise(
    age_index = dplyr::first(Age),
    glaucoma = dplyr::first(GlaucomaType),
    reop_any = if (has_reop) suppressWarnings(max(Reop_num, na.rm = TRUE)) else NA_real_,
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    reop_any = ifelse(is.infinite(reop_any), NA, reop_any),
    AgeGroup = dplyr::case_when(
      is.na(age_index) ~ NA_character_,
      age_index < 40 ~ "<40",
      age_index < 60 ~ "40-59",
      age_index < 80 ~ "60-79",
      TRUE ~ ">=80"
    ),
    Reop_any = dplyr::case_when(
      !has_reop ~ NA_character_,
      is.na(reop_any) ~ NA_character_,
      reop_any >= 1 ~ "Yes",
      TRUE ~ "No"
    )
  ) %>%
  dplyr::mutate(
    AgeGroup = factor(AgeGroup, levels = c("<40", "40-59", "60-79", ">=80")),
    glaucoma = factor(as.character(glaucoma), levels = c("PACG", "POAG", "SG", "PG", "Others")),
    Reop_any = if (has_reop) factor(Reop_any, levels = c("No", "Yes")) else NULL
  )

cluster_df <- dat_pat %>%
  dplyr::select(patient_id, AgeGroup, glaucoma, dplyr::any_of("Reop_any")) %>%
  dplyr::filter(!is.na(AgeGroup), !is.na(glaucoma)) %>%
  { if (has_reop) dplyr::filter(., !is.na(Reop_any)) else . }

vars_clust <- c("AgeGroup", "glaucoma")
if (has_reop) vars_clust <- c(vars_clust, "Reop_any")

gower_dist <- cluster::daisy(cluster_df[, vars_clust, drop = FALSE], metric = "gower")

ks <- 2:6
avg_sil <- sapply(ks, function(k) {
  fit_k <- cluster::pam(gower_dist, k = k, diss = TRUE)
  fit_k$silinfo$avg.width
})

k_final <- ks[which.max(avg_sil)]
pam_fit <- cluster::pam(gower_dist, k = k_final, diss = TRUE)
cluster_df$Cluster <- factor(pam_fit$clustering, labels = paste0("C", seq_len(k_final)))

cluster_size <- cluster_df %>%
  dplyr::count(Cluster, name = "n_patients") %>%
  dplyr::mutate(pct = round(n_patients / sum(n_patients) * 100, 1))

profile_age <- cluster_df %>%
  dplyr::count(Cluster, AgeGroup) %>%
  dplyr::group_by(Cluster) %>%
  dplyr::mutate(pct = round(n / sum(n) * 100, 1)) %>%
  dplyr::ungroup()

profile_glaucoma <- cluster_df %>%
  dplyr::count(Cluster, glaucoma) %>%
  dplyr::group_by(Cluster) %>%
  dplyr::mutate(pct = round(n / sum(n) * 100, 1)) %>%
  dplyr::ungroup()

if (has_reop) {
  profile_reop <- cluster_df %>%
    dplyr::count(Cluster, Reop_any) %>%
    dplyr::group_by(Cluster) %>%
    dplyr::mutate(pct = round(n / sum(n) * 100, 1)) %>%
    dplyr::ungroup()
}

dat_reg <- dat_proc %>%
  dplyr::left_join(cluster_df %>% dplyr::select(patient_id, Cluster), by = "patient_id") %>%
  dplyr::filter(!is.na(Cluster), !is.na(Surgery_bin), !is.na(Sex), !is.na(Year_c))

ext_valid_internal <- dat_reg %>%
  dplyr::group_by(Cluster) %>%
  dplyr::summarise(
    n_procedures = dplyr::n(),
    internal_n = sum(Surgery_bin == 1, na.rm = TRUE),
    internal_pct = round(internal_n / n_procedures * 100, 1),
    .groups = "drop"
  ) %>%
  dplyr::arrange(internal_pct)

# Use the cluster with the lowest internal procedure proportion as the reference,
# but do not hard-code a specific cluster label because labels are data-dependent.
dat_reg <- dat_reg %>%
  dplyr::mutate(
    Cluster = factor(Cluster, levels = ext_valid_internal$Cluster),
    Sex = factor(Sex, levels = c("Male", "Female"))
  )

fit_cluster <- glm(
  Surgery_bin ~ Cluster + Sex + Year_c,
  data = dat_reg,
  family = binomial
)

fit_cluster_gee <- geepack::geeglm(
  Surgery_bin ~ Cluster + Sex + Year_c,
  data = dat_reg,
  id = patient_id,
  family = binomial(link = "logit"),
  corstr = "exchangeable"
)

beta <- stats::coef(fit_cluster_gee)
vc_robust <- fit_cluster_gee$geese$vbeta
se_robust <- sqrt(diag(vc_robust))

gee_or <- data.frame(
  term = names(beta),
  OR = exp(beta),
  LCL = exp(beta - 1.96 * se_robust),
  UCL = exp(beta + 1.96 * se_robust)
)

p_sil <- ggplot2::ggplot(data.frame(k = ks, avg_silhouette = avg_sil),
                         ggplot2::aes(x = k, y = avg_silhouette)) +
  ggplot2::geom_line() +
  ggplot2::geom_point(size = 2) +
  ggplot2::scale_x_continuous(breaks = ks) +
  ggplot2::labs(x = "k", y = "Average silhouette width",
                title = "Silhouette analysis (PAM on Gower distance)") +
  ggplot2::theme_classic(base_size = 12)

ggplot2::ggsave(file.path(output_dir, "silhouette_curve.pdf"), p_sil, width = 6, height = 4, units = "in")
ggplot2::ggsave(file.path(output_dir, "silhouette_curve.tiff"), p_sil, width = 6, height = 4, units = "in",
                dpi = 600, compression = "lzw")

writexl::write_xlsx(
  list(
    cluster_size = cluster_size,
    profile_age = profile_age,
    profile_glaucoma = profile_glaucoma,
    profile_reop = if (exists("profile_reop")) profile_reop else data.frame(),
    internal_pct_by_cluster = ext_valid_internal,
    logistic_cluster = broom::tidy(fit_cluster, conf.int = TRUE, exponentiate = TRUE),
    gee_cluster = gee_or
  ),
  path = file.path(output_dir, "clustering_outputs.xlsx")
)

message("Clustering analysis complete. Outputs saved in: ", output_dir)
