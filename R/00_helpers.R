# 00_helpers.R
# Shared utilities for the glaucoma-surgery-trends repository.
# This file centralizes package loading, date parsing, procedure-code
# standardization, and analytic cohort construction.

load_required_packages <- function(pkgs) {
  missing_pkgs <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (length(missing_pkgs) > 0) install.packages(missing_pkgs)
  invisible(lapply(pkgs, library, character.only = TRUE))
}

required_pkgs <- c(
  "readxl", "dplyr", "tidyr", "stringr", "ggplot2", "broom",
  "writexl", "pROC", "ResourceSelection", "binom", "scales",
  "splines", "sandwich", "lmtest", "cluster", "geepack"
)
load_required_packages(required_pkgs)

YEAR_MIN <- 2015L
YEAR_MAX <- 2024L

INTERNAL_CODES <- c(
  "GSL", "GSL+GT", "P+I+CLASS", "P+I+GATT", "3T", "ABIC", "CLASS",
  "GATT", "GT", "MAT", "P+I", "P+I+GSL", "P+I+GSL+GT", "P+I+GT",
  "PCP", "P", "SPI"
)

EXTERNAL_CODES <- c(
  "GDIS", "P+I+GDIS", "P+I+TRAB", "TRAB", "TRAB+TRABECULOTOMY"
)

OTHER_CODES <- c("光凝", "冷冻", "UCP", "滤过泡修补")

parse_any_date <- function(x) {
  if (inherits(x, "Date") || inherits(x, "POSIXt")) return(as.Date(x))
  if (is.numeric(x)) return(as.Date(x, origin = "1899-12-30"))

  x_chr <- as.character(x)
  x_chr <- stringr::str_trim(x_chr)
  x_chr[x_chr == "" | is.na(x_chr)] <- NA_character_

  is_num_str <- !is.na(x_chr) & stringr::str_detect(x_chr, "^\\d{4,6}$")
  out <- rep(as.Date(NA), length(x_chr))

  if (any(is_num_str)) {
    out[is_num_str] <- as.Date(as.numeric(x_chr[is_num_str]), origin = "1899-12-30")
  }

  idx <- which(!is_num_str & !is.na(x_chr))
  if (length(idx) > 0) {
    s <- x_chr[idx]
    date_part <- stringr::str_extract(
      s,
      "(19\\d{2}|20\\d{2})[年\\-/\\.](\\d{1,2})[月\\-/\\.](\\d{1,2})"
    )
    date_part2 <- stringr::str_extract(
      s,
      "(19\\d{2}|20\\d{2})年(\\d{1,2})月(\\d{1,2})日"
    )
    date_part[is.na(date_part)] <- date_part2[is.na(date_part)]

    date_norm <- date_part
    date_norm <- gsub("年", "-", date_norm)
    date_norm <- gsub("月", "-", date_norm)
    date_norm <- gsub("日", "", date_norm)
    date_norm <- gsub("[/\\.]", "-", date_norm)
    out[idx] <- as.Date(date_norm, format = "%Y-%m-%d")
  }

  out
}

extract_year_fallback <- function(x) {
  y <- stringr::str_extract(as.character(x), "\\b(19\\d{2}|20\\d{2})\\b")
  suppressWarnings(as.integer(y))
}

standardize_surgery_code <- function(x) {
  s <- stringr::str_trim(as.character(x))
  s[s == "" | is.na(s)] <- NA_character_
  s <- stringr::str_replace_all(s, "–|—", "-")
  s <- stringr::str_replace_all(s, "\\s+", "")

  s <- dplyr::if_else(
    !is.na(s) & (toupper(s) == "UCP" | stringr::str_detect(s, "超声.*睫状体.*成形|睫状体成形")),
    "UCP", s
  )

  s <- dplyr::if_else(
    !is.na(s) & stringr::str_detect(s, "冷凝|冷冻|睫状体冷凝|经巩膜.*冷"),
    "冷冻", s
  )

  s <- dplyr::if_else(
    !is.na(s) & (
      stringr::str_detect(toupper(s), "TSCPC|CPC") |
        stringr::str_detect(s, "光凝|激光|微脉冲|经巩膜.*光")
    ),
    "光凝", s
  )

  s_up <- toupper(s)
  s <- dplyr::if_else(s %in% c("冷冻", "光凝"), s, s_up)
  s <- dplyr::if_else(!is.na(s) & s == "PHACO", "P", s)
  s
}

classify_glaucoma_type <- function(dx) {
  dx_upper <- toupper(as.character(dx))
  dplyr::case_when(
    is.na(dx) ~ "Others",
    dx_upper == "PACG" ~ "PACG",
    dx_upper == "POAG" ~ "POAG",
    dx_upper %in% c("CG", "PG", "PCG", "JOAG") ~ "PG",
    dx_upper == "SG" ~ "SG",
    TRUE ~ "Others"
  )
}

build_analytic_cohort <- function(input_file, year_min = YEAR_MIN, year_max = YEAR_MAX) {
  df_raw <- readxl::read_xlsx(input_file)

  df_clean <- df_raw %>%
    dplyr::mutate(
      SurgeryType_raw = 手术,
      SurgeryType_std = standardize_surgery_code(手术),
      GlaucomaType = classify_glaucoma_type(主要诊断),
      Sex = dplyr::recode(性别, "男" = "Male", "女" = "Female", .default = NA_character_),
      Age = suppressWarnings(as.numeric(年龄)),
      Age10 = Age / 10,
      Reop = dplyr::recode(是否再次手术, "是" = 1L, "否" = 0L, .default = NA_integer_),
      SurgeryDate1 = parse_any_date(手术及操作日期1),
      SurgeryYear1 = suppressWarnings(as.integer(format(SurgeryDate1, "%Y"))),
      SurgeryYear1 = dplyr::if_else(
        is.na(SurgeryYear1),
        extract_year_fallback(手术及操作日期1),
        SurgeryYear1
      ),
      AdmDate1 = parse_any_date(入院时间),
      AdmYear = suppressWarnings(as.integer(format(AdmDate1, "%Y"))),
      AdmYear = dplyr::if_else(is.na(AdmYear), extract_year_fallback(入院时间), AdmYear),
      Year_model = dplyr::if_else(!is.na(SurgeryYear1), SurgeryYear1, AdmYear),
      Year_c = Year_model - year_min
    ) %>%
    dplyr::filter(!is.na(SurgeryType_std)) %>%
    dplyr::mutate(
      SurgeryClass = dplyr::case_when(
        SurgeryType_std %in% INTERNAL_CODES ~ "Internal",
        SurgeryType_std %in% EXTERNAL_CODES ~ "External",
        SurgeryType_std %in% OTHER_CODES ~ "Other",
        TRUE ~ "Other"
      )
    ) %>%
    dplyr::filter(!is.na(Year_model), Year_model >= year_min, Year_model <= year_max)

  df_analysis <- df_clean %>%
    dplyr::select(
      患者唯一号, 病案号, 登记号, 就诊ID,
      Sex, Age, Age10, GlaucomaType,
      SurgeryType_raw, SurgeryType_std, SurgeryClass,
      手术及操作日期1, SurgeryDate1, SurgeryYear1,
      入院时间, AdmYear, Year_model, Year_c,
      Reop, 眼别, 出院时间
    )

  df_model <- df_analysis %>%
    dplyr::filter(SurgeryClass %in% c("Internal", "External")) %>%
    dplyr::mutate(
      Surgery_bin = dplyr::if_else(SurgeryClass == "Internal", 1L, 0L),
      Sex = factor(Sex, levels = c("Male", "Female")),
      GlaucomaType = factor(GlaucomaType, levels = c("PACG", "POAG", "SG", "PG", "Others")),
      Year_f = factor(Year_model)
    )

  df_model_cc <- df_model %>%
    dplyr::filter(
      !is.na(Surgery_bin),
      !is.na(Age10),
      !is.na(Sex),
      !is.na(GlaucomaType),
      !is.na(Year_c)
    )

  list(
    df_analysis = df_analysis,
    df_model = df_model,
    df_model_cc = df_model_cc
  )
}

make_output_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}
