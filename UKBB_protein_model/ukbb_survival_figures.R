library(survminer)
library(survival)
library(data.table)
library(dplyr)
library(tidyr)
library(gridExtra)
library(forestmodel)
library(ggpubr)
library(scales)

`%!in%` <- function(x, y) !(`%in%`(x, y))

cli_defaults <- list(
  ukbb_metadata = "UKBB_FINNGEN_test_condition_excl_proteom.tsv",
  excluded_samples = "w88907_20250818.csv",
  predictions = "pred_pheb_combat_1_lasso_stlseakclrah_v3.txt",
  participant_file = "participant_1.tsv.gz",
  output_dir = "FIGURES_STLSEAKCLRAH_PHEB_NONEGATIVE_UKBB_excl"
)

print_usage <- function(defaults) {
  arg_names <- names(defaults)
  arg_lines <- paste0(
    "  --", gsub("_", "-", arg_names), " <value>",
    "    default: ", unlist(defaults, use.names = FALSE)
  )

  cat(paste(c(
    "Usage:",
    "  Rscript ukbb_survival_figures_public.R [arguments]",
    "",
    "Optional arguments:",
    arg_lines,
    "",
    "Arguments can be passed as --name value or --name=value. Hyphens and underscores are equivalent in argument names."
  ), collapse = "\n"), "\n")
}

parse_cli_args <- function(defaults) {
  raw_args <- commandArgs(trailingOnly = TRUE)

  if (any(raw_args %in% c("-h", "--help"))) {
    print_usage(defaults)
    quit(status = 0)
  }

  parsed <- list()
  i <- 1
  while (i <= length(raw_args)) {
    token <- raw_args[[i]]
    if (!startsWith(token, "--")) {
      stop(paste0("Unexpected argument format: ", token), call. = FALSE)
    }

    token <- sub("^--", "", token)
    if (grepl("=", token, fixed = TRUE)) {
      pieces <- strsplit(token, "=", fixed = TRUE)[[1]]
      key <- pieces[[1]]
      value <- paste(pieces[-1], collapse = "=")
    } else {
      key <- token
      if (i == length(raw_args) || startsWith(raw_args[[i + 1]], "--")) {
        stop(paste0("Missing value for --", key), call. = FALSE)
      }
      value <- raw_args[[i + 1]]
      i <- i + 1
    }

    key <- gsub("-", "_", key)
    if (!key %in% names(defaults)) {
      stop(paste0("Unknown argument: --", gsub("_", "-", key)), call. = FALSE)
    }

    parsed[[key]] <- value
    i <- i + 1
  }

  values <- defaults
  for (key in names(parsed)) {
    values[[key]] <- parsed[[key]]
  }

  values
}

cli <- parse_cli_args(cli_defaults)

ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
}

save_plot_pdf <- function(plot, path, width, height) {
  pdf(path, width = width, height = height)
  print(plot)
  dev.off()
}

save_plot_png <- function(plot, path, width, height, res = 300) {
  png(path, width = width, height = height, units = "in", res = res)
  print(plot)
  dev.off()
}

data4 <- as.data.frame(fread(cli$ukbb_metadata, header = TRUE))
EXCLUDED <- read.table(cli$excluded_samples, sep = "\t", header = FALSE)
colnames(EXCLUDED)[1] <- "eid"
data4 <- data4 %>% filter(eid %!in% EXCLUDED$eid)

pheb1000 <- as.data.frame(fread(cli$predictions, header = TRUE, sep = "\t"))
colnames(pheb1000)[2] <- "pheb1000"

p1 <- as.data.frame(fread(cli$participant_file, header = TRUE, sep = "\t"))
p1 <- p1 %>% select(eid, p53_i0, p191, p52, p34)
p1$p34 <- paste(p1$p34, p1$p52, "15", sep = "-")
p1 <- p1 %>% mutate(
  AGE = as.numeric(
    difftime(
      as.POSIXct(as.Date(p1$p53_i0, "%Y-%m-%d")),
      as.POSIXct(as.Date(p1$p34, "%Y-%m-%d")),
      units = "weeks"
    )
  ) / 52.25
)

data4 <- as.data.frame(merge(data4, p1, by = "eid"))
data4 <- as.data.frame(merge(data4, pheb1000, by = "eid"))
data4 <- data4 %>% mutate(pheb1000 = case_when(pheb1000 < 0 ~ 0, TRUE ~ pheb1000))

data1 <- data4 %>% select(
  eid, pheb1000, SEX, BL_AGE,
  E4_ENDONUTRMET, I9_HYPERTENSION, I9_CVD, E4_DM2, M13_POLYARTHROPATHIES,
  J10_COPD, AUTOIMMUNE, M13_SLE, K11_LIVER, J10_ASTHMA, I9_MI,
  E4_ENDONUTRMET_AGE, I9_HYPERTENSION_AGE, I9_CVD_AGE, E4_DM2_AGE,
  M13_POLYARTHROPATHIES_AGE, J10_COPD_AGE, AUTOIMMUNE_AGE, M13_SLE_AGE,
  K11_LIVER_AGE, J10_ASTHMA_AGE, I9_MI_AGE
)

KM25_age <- function(data1, protein, pheno, SEX1, SEX2, phenostring = "protein_level", pos1, mod, mod1) {
  data1$protein <- data1[, protein]
  data1$PHENO <- data1[, paste(pheno, sep = "")]
  d <- data1 %>% select(eid, pheb1000, BL_AGE, DEATH, DEATH_AGE)

  if (mod == "normal") {
    data1$PHENO_AGE <- data1[, paste(pheno, "_AGE", sep = "")]
  } else {
    data1$PHENO_AGE <- data1$AGE
  }

  data1$DEATH_AGE_DIFF <- data1[, paste(pheno, "_AGE", sep = "")] - data1$BL_AGE
  data1 <- data1 %>% filter(DEATH_AGE_DIFF > 0 & DEATH_AGE_DIFF < 17)
  data1 <- data1 %>% filter(!is.na(protein))
  data1 <- data1 %>% filter(!is.na(PHENO))
  data1 <- data1 %>% filter(SEX %in% c(SEX1, SEX2))

  if (mod1 == "paper") {
    quantiles <- data1 %>%
      group_by(BL_AGE) %>%
      dplyr::summarise(
        q10 = quantile(protein, 0.10, na.rm = TRUE),
        q90 = quantile(protein, 0.90, na.rm = TRUE),
        q45 = quantile(protein, 0.45, na.rm = TRUE),
        q55 = quantile(protein, 0.55, na.rm = TRUE)
      )
    data1 <- data1 %>%
      left_join(quantiles, by = "BL_AGE") %>%
      mutate(
        protein_level_25 = case_when(
          protein < q10 ~ "0-10%",
          protein > q45 & protein < q55 ~ "45-55%",
          protein > q90 ~ "90-100%"
        )
      )
  } else {
    quantiles <- data1 %>%
      group_by(BL_AGE) %>%
      dplyr::summarise(
        q25 = quantile(protein, 0.25, na.rm = TRUE),
        q75 = quantile(protein, 0.75, na.rm = TRUE)
      )

    data1 <- data1 %>%
      left_join(quantiles, by = "BL_AGE") %>%
      mutate(
        protein_level_25 = case_when(
          protein < q25 ~ "0-25%",
          protein >= q25 & protein < q75 ~ "25-75%",
          protein >= q75 ~ "75-100%"
        )
      )
  }

  data1 <- data1[!is.na(data1$protein_level_25), ]
  fit <- survival::survfit(Surv(data1$PHENO_AGE, data1$PHENO) ~ protein_level_25, data = data1)

  data1 <<- data1

  pval_result <- survminer::surv_pvalue(fit, data = data1)
  pval_text <- paste(pval_result$pval.txt)
  summ <- summary(fit, times = pos1)
  event_at_65 <- 1 - summ$surv
  pos2 <- min(1, max(event_at_65, na.rm = TRUE))
  pos2 <- pos2 + pos2 * 0.5
  pos <- pos2

  pheno_names <- c(
    "AB1_SEPSIS" = "Septicaemia",
    "AB1_SEPSIS_CONDITION" = "Condition for Implicit Sepsis (Organ dysfunction codes)",
    "AUTOIMMUNE" = "Autoimmune diseases",
    "C3_CANCER" = "Malignant neoplasm",
    "D3_ANAEMIA" = "Anaemias",
    "DEATH" = "Death",
    "E4_DM2" = "Type II diabetes",
    "F5_DEMENTIA" = "Dementia",
    "FG_HYPERTENSION" = "Hypertensive diseases (excluding secondary)",
    "G6_ALZHEIMER" = "Alzheimer disease",
    "G6_PARKINSON" = "Parkinson's disease",
    "J10_PNEUMONIA" = "All pneumoniae",
    "K11_LIVER" = "Diseases of liver",
    "M13_GOUT" = "Gout",
    "N14_RENFAIL" = "Renal failure",
    "H7_MACULADEGEN" = "Degeneration of macula and posterior pole",
    "M13_OSTEOPOROSIS" = "Osteoporosis",
    "C3_LUNG_NONSMALL_EXALLC" = "Non-small cell lung cancer",
    "M13_SOFTOVERUSE" = "Soft tissue disorders related to use, overuse and pressure",
    "F5_PERSOTH" = "Other and unspecified disorders of adult personality and behaviour",
    "E4_DIABINSIPIDUS" = "Diabetes insipidus"
  )

  title_text <- ifelse(pheno %in% names(pheno_names), pheno_names[pheno], pheno)

  if (mod1 == "paper" & mod == "normal") {
    g <- ggsurvplot(
      fit,
      data = data1,
      conf.int = TRUE,
      risk.table = TRUE,
      xlim = c(40, pos1),
      ylim = c(0, pos2),
      risk.table.col = "strata",
      fun = "event",
      theme = theme_bw(),
      legend = "none",
      palette = c("#4F7942", "#999999", "#CC5500"),
      legend.labs = c(
        paste(phenostring, "=0-10%", sep = ""),
        paste(phenostring, "=45-55%", sep = ""),
        paste(phenostring, "=90-100%", sep = "")
      ),
      ylab = paste("Probability", sep = ""),
      xlab = "AGE OF EVENT",
      title = title_text
    )
    g$plot <- g$plot + theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))
    g$table <- g$table + theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))
  } else if (mod1 == "paper" & mod != "normal") {
    g <- ggsurvplot(
      fit,
      data = data1,
      conf.int = TRUE,
      risk.table = TRUE,
      xlim = c(40, pos1),
      ylim = c(0, pos2),
      risk.table.col = "strata",
      fun = "event",
      theme = theme_bw(),
      legend = "none",
      palette = c("#4F7942", "#999999", "#CC5500"),
      legend.labs = c(
        paste(phenostring, "=0-10%", sep = ""),
        paste(phenostring, "=45-55%", sep = ""),
        paste(phenostring, "=90-100%", sep = "")
      ),
      ylab = paste("Probability", sep = ""),
      xlab = "AGE OF ENROLLMENT",
      title = title_text
    )
    g$plot <- g$plot + theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))
    g$table <- g$table + theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))
  } else if (mod1 != "paper" & mod == "normal") {
    g <- ggsurvplot(
      fit,
      data = data1,
      conf.int = TRUE,
      risk.table = TRUE,
      xlim = c(40, pos1),
      ylim = c(0, pos2),
      risk.table.col = "strata",
      fun = "event",
      theme = theme_bw(),
      legend = "none",
      palette = c("#4F7942", "#999999", "#CC5500"),
      legend.labs = c(
        paste(phenostring, "=0-25%", sep = ""),
        paste(phenostring, "=25-75%", sep = ""),
        paste(phenostring, "=75-100%", sep = "")
      ),
      ylab = paste("Probability", sep = ""),
      xlab = "AGE OF EVENT",
      title = title_text
    )
    g$plot <- g$plot + theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))
    g$table <- g$table + theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))
  } else if (mod1 != "paper" & mod != "normal") {
    g <- ggsurvplot(
      fit,
      data = data1,
      conf.int = TRUE,
      risk.table = TRUE,
      xlim = c(40, pos1),
      ylim = c(0, pos2),
      risk.table.col = "strata",
      fun = "event",
      theme = theme_bw(),
      legend = "none",
      palette = c("#4F7942", "#999999", "#CC5500"),
      legend.labs = c(
        paste(phenostring, "=0-25%", sep = ""),
        paste(phenostring, "=25-75%", sep = ""),
        paste(phenostring, "=75-100%", sep = "")
      ),
      ylab = paste("Probability", sep = ""),
      xlab = "AGE OF ENROLLMENT",
      title = title_text
    )
    g$plot <- g$plot + theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))
    g$table <- g$table + theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))
  }

  g$plot <- g$plot + annotate("text", x = 50, y = pos, label = pval_text, size = 5, hjust = 0)
  return(g$plot)
}

KM25_group <- function(data1, protein, pheno, SEX1, SEX2, phenostring = "protein_level", AGE_) {
  data1$protein <- data1[, protein]
  data1$PHENO <- data1[, paste(pheno, sep = "")]
  data1$PHENO_AGE <- data1[, paste(pheno, "_AGE", sep = "")]
  data1$DEATH_AGE_DIFF <- data1[, paste(pheno, "_AGE", sep = "")] - data1$BL_AGE

  data1 <- data1 %>% filter(DEATH_AGE_DIFF > 0 & DEATH_AGE_DIFF < 17)
  data1 <- data1 %>% filter(!is.na(protein))
  data1 <- data1 %>% filter(!is.na(PHENO))
  data1 <- data1 %>% filter(SEX %in% c(SEX1, SEX2))
  data1 <- data1 %>% mutate(AGE_cat = case_when(
    BL_AGE >= 40 & BL_AGE < 50 ~ "40-50",
    BL_AGE >= 50 & BL_AGE < 60 ~ "50-60",
    BL_AGE >= 60 & BL_AGE < 70 ~ "60-70"
  ))

  data1 <- data1 %>% filter(AGE_cat == AGE_)

  quantiles <- data1 %>%
    group_by(BL_AGE) %>%
    dplyr::summarise(
      q25 = quantile(protein, 0.25, na.rm = TRUE),
      q75 = quantile(protein, 0.75, na.rm = TRUE)
    )

  data1 <- data1 %>%
    left_join(quantiles, by = "BL_AGE") %>%
    mutate(
      protein_level_25 = case_when(
        protein < q25 ~ "0-25%",
        protein >= q25 & protein < q75 ~ "25-75%",
        protein >= q75 ~ "75-100%"
      )
    )

  fit <- survival::survfit(Surv(data1$DEATH_AGE_DIFF, data1$PHENO) ~ protein_level_25, data = data1)

  data1 <<- data1

  summ <- summary(fit, times = 16)
  event_at_65 <- 1 - summ$surv
  pos <- min(1, max(event_at_65, na.rm = TRUE))
  pos <- pos + pos * 0.5

  pheno_names <- c(
    "AB1_SEPSIS" = "Septicaemia",
    "AB1_SEPSIS_CONDITION" = "Condition for Implicit Sepsis (Organ dysfunction codes)",
    "AUTOIMMUNE" = "Autoimmune diseases",
    "C3_CANCER" = "Malignant neoplasm",
    "D3_ANAEMIA" = "Anaemias",
    "DEATH" = "Death",
    "E4_DM2" = "Type II diabetes",
    "F5_DEMENTIA" = "Dementia",
    "FG_HYPERTENSION" = "Hypertensive diseases (excluding secondary)",
    "G6_ALZHEIMER" = "Alzheimer disease",
    "G6_PARKINSON" = "Parkinson's disease",
    "J10_PNEUMONIA" = "All pneumoniae",
    "K11_LIVER" = "Diseases of liver",
    "M13_GOUT" = "Gout",
    "N14_RENFAIL" = "Renal failure",
    "H7_MACULADEGEN" = "Degeneration of macula and posterior pole",
    "M13_OSTEOPOROSIS" = "Osteoporosis",
    "C3_LUNG_NONSMALL_EXALLC" = "Non-small cell lung cancer",
    "M13_SOFTOVERUSE" = "Soft tissue disorders related to use, overuse and pressure",
    "F5_PERSOTH" = "Other and unspecified disorders of adult personality and behaviour",
    "E4_DIABINSIPIDUS" = "Diabetes insipidus"
  )

  title_text <- ifelse(pheno %in% names(pheno_names), pheno_names[pheno], pheno)

  g1 <- ggsurvplot(
    fit,
    data = data1,
    pval = TRUE,
    pval.coord = c(4, pos),
    conf.int = TRUE,
    risk.table = TRUE,
    risk.table.col = "strata",
    fun = "event",
    theme = theme_bw(),
    legend = "none",
    palette = c("#4F7942", "#999999", "#CC5500"),
    legend.labs = c(
      paste(phenostring, "=0-25%", " AGE_cat=", AGE_, sep = ""),
      paste(phenostring, "=25-75%", " AGE_cat=", AGE_, sep = ""),
      paste(phenostring, "=75-100%", " AGE_cat=", AGE_, sep = "")
    ),
    ylab = paste("Probability", sep = ""),
    xlab = "Follow up age",
    title = title_text
  )

  g1$plot <- g1$plot + theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))
  g1$table <- g1$table + theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))

  return(g1$plot)
}

BOX <- function(data1, protein, pheno, SEX1, SEX2) {
  data1$protein <- data1[, protein]
  data1$PHENO <- data1[, paste(pheno, sep = "")]
  data1$DEATH_AGE_DIFF <- data1[, paste(pheno, "_AGE", sep = "")] - data1$BL_AGE
  data1 <- data1 %>% filter(!is.na(protein))
  data1 <- data1 %>% filter(SEX %in% c(SEX1, SEX2))

  quantiles <- data1 %>%
    group_by(BL_AGE) %>%
    dplyr::summarise(
      q25 = quantile(protein, 0.25, na.rm = TRUE),
      q75 = quantile(protein, 0.75, na.rm = TRUE)
    )

  data1 <- data1 %>%
    left_join(quantiles, by = "BL_AGE") %>%
    mutate(
      protein_level_25 = case_when(
        protein < q25 ~ "0-25%",
        protein >= q25 & protein < q75 ~ "25-75%",
        protein >= q75 ~ "75-100%"
      )
    )

  data1 <- data1 %>% mutate(AGE_cat = case_when(
    BL_AGE >= 40 & BL_AGE < 50 ~ "40-50",
    BL_AGE >= 50 & BL_AGE < 60 ~ "50-60",
    BL_AGE >= 60 & BL_AGE < 70 ~ "60-70"
  ))

  data1 <- data1 %>% mutate(PHENO1 = case_when(DEATH_AGE_DIFF > 0 & PHENO == 1 ~ 0, TRUE ~ PHENO))

  cols <- c("#009E73", "#E69F00", "#999999", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

  data1 %>%
    mutate(PHENO1 = case_when(PHENO1 == 1 ~ "yes", PHENO1 == 0 ~ "no")) %>%
    filter(!is.na(AGE_cat) & !is.na(PHENO1)) %>%
    ggplot(aes(AGE_cat, protein, fill = PHENO1)) +
    geom_boxplot(outlier.size = 0.1, size = 0.1) +
    theme_classic() +
    ylab(paste(protein)) +
    xlab(paste("AGE group")) +
    scale_fill_manual(values = cols) +
    labs(fill = paste(pheno)) +
    stat_compare_means(
      aes(group = PHENO1),
      method = "t.test",
      label = "p.signif",
      hide.ns = TRUE
    ) +
    stat_compare_means(
      comparisons = list(c("40-50", "50-60"), c("50-60", "60-70")),
      method = "t.test",
      label = "p.signif",
      hide.ns = TRUE,
      label.y = c(max(data1$protein), max(data1$protein) + 1, max(data1$protein) + 2)
    )
}

BOX1 <- function(data1, protein, SEX1, SEX2) {
  data1$protein <- data1[, protein]
  data1 <- data1 %>% filter(!is.na(protein))
  data1 <- data1 %>% filter(SEX %in% c(SEX1, SEX2))

  quantiles <- data1 %>%
    group_by(BL_AGE) %>%
    dplyr::summarise(
      q25 = quantile(protein, 0.25, na.rm = TRUE),
      q75 = quantile(protein, 0.75, na.rm = TRUE)
    )

  data1 <- data1 %>%
    left_join(quantiles, by = "BL_AGE") %>%
    mutate(
      protein_level_25 = case_when(
        protein < q25 ~ "0-25%",
        protein >= q25 & protein < q75 ~ "25-75%",
        protein >= q75 ~ "75-100%"
      )
    )

  data1 <- data1 %>% mutate(AGE_cat = case_when(
    BL_AGE >= 40 & BL_AGE < 45 ~ "40-45",
    BL_AGE >= 45 & BL_AGE < 50 ~ "45-50",
    BL_AGE >= 50 & BL_AGE < 55 ~ "50-55",
    BL_AGE >= 55 & BL_AGE < 60 ~ "55-60",
    BL_AGE >= 60 & BL_AGE < 65 ~ "60-65",
    BL_AGE >= 65 & BL_AGE < 70 ~ "65-70"
  ))

  m <- cor.test(data1$protein, data1$BL_AGE, method = "pearson")
  b_p <- m$estimate
  p_p <- m$p.value

  data1 %>%
    filter(!is.na(AGE_cat)) %>%
    ggplot(aes(AGE_cat, protein)) +
    geom_boxplot(outlier.size = 0.1, size = 0.1) +
    theme_classic() +
    ylab(paste(protein)) +
    xlab(paste("AGE group")) +
    stat_compare_means(
      comparisons = list(c("40-45", "45-50"), c("45-50", "50-55"), c("50-55", "55-60"), c("55-60", "60-65"), c("60-65", "65-70")),
      method = "t.test",
      label = "p.signif",
      hide.ns = TRUE,
      label.y = c(15, 17, 19, 21, 23)
    ) +
    annotate("text", x = 0.5, y = 23, label = paste("R2: ", b_p, sep = ""), size = 4, hjust = 0) +
    annotate("text", x = 0.5, y = 21, label = paste("P: ", p_p, sep = ""), size = 4, hjust = 0)
}

COX <- function(data4, protein, pheno, age, SEX1, SEX2, var = "protein_level") {
  pheno_names <- c(
    "AB1_SEPSIS" = "Septicaemia",
    "AB1_SEPSIS_CONDITION" = "Condition for Implicit Sepsis (Organ dysfunction codes)",
    "AUTOIMMUNE" = "Autoimmune diseases",
    "C3_CANCER" = "Malignant neoplasm",
    "D3_ANAEMIA" = "Anaemias",
    "DEATH" = "Death",
    "E4_DM2" = "Type II diabetes",
    "F5_DEMENTIA" = "Dementia",
    "FG_HYPERTENSION" = "Hypertensive diseases (excluding secondary)",
    "G6_ALZHEIMER" = "Alzheimer disease",
    "G6_PARKINSON" = "Parkinson's disease",
    "J10_PNEUMONIA" = "All pneumoniae",
    "K11_LIVER" = "Diseases of liver",
    "M13_GOUT" = "Gout",
    "N14_RENFAIL" = "Renal failure",
    "H7_MACULADEGEN" = "Degeneration of macula and posterior pole",
    "M13_OSTEOPOROSIS" = "Osteoporosis",
    "C3_LUNG_NONSMALL_EXALLC" = "Non-small cell lung cancer",
    "M13_SOFTOVERUSE" = "Soft tissue disorders related to use, overuse and pressure",
    "F5_PERSOTH" = "Other and unspecified disorders of adult personality and behaviour",
    "E4_DIABINSIPIDUS" = "Diabetes insipidus"
  )
  title_text <- ifelse(pheno %in% names(pheno_names), pheno_names[pheno], pheno)

  data4$DEATH_AGE_DIFF <- data4[, paste(pheno, "_AGE", sep = "")] - data4$BL_AGE
  data4$protein <- data4[, protein]

  data1 <- data4 %>% filter(DEATH_AGE_DIFF > 0 & DEATH_AGE_DIFF < 18)
  data1 <- data1 %>% filter(!is.na(protein))

  su <- Surv(data1$DEATH_AGE_DIFF, data1[, paste(pheno, sep = "")])

  data1 <- data1 %>% mutate(AGE_cat = case_when(
    BL_AGE >= 40 & BL_AGE < 50 ~ "40-50",
    BL_AGE >= 50 & BL_AGE < 60 ~ "50-60",
    BL_AGE >= 60 & BL_AGE < 70 ~ "60-70"
  ))

  data1 <- data1 %>% filter(AGE_cat == age & SEX %in% c(SEX1, SEX2))

  quantiles <- data1 %>%
    group_by(BL_AGE) %>%
    dplyr::summarise(
      q25 = quantile(protein, 0.25, na.rm = TRUE),
      q75 = quantile(protein, 0.75, na.rm = TRUE)
    )

  data1 <- data1 %>%
    left_join(quantiles, by = "BL_AGE") %>%
    mutate(
      protein_level_25 = case_when(
        protein < q25 ~ "0-25%",
        protein >= q25 & protein < q75 ~ "25-75%",
        protein >= q75 ~ "75-100%"
      )
    )

  names(data1)[names(data1) == "protein_level_25"] <- var
  data1$SEX <- ifelse(data1$SEX == 1, "males", "females")

  names(data1)[names(data1) == "BL_AGE"] <- "AGE"
  names(data1)[names(data1) == var] <- "%TemB Percentile"

  su <- Surv(data1$DEATH_AGE_DIFF, data1[, paste(pheno, sep = "")])
  v_bt <- "`%TemB Percentile`"

  if (SEX1 == 0 | SEX2 == 0) {
    model <- as.formula(paste("su", paste(v_bt, "AGE", sep = "+"), sep = "~"))
    m <- coxph(model, data = data1)

    g1 <- forest_model(
      m,
      covariates = c("`%TemB Percentile`", "AGE"),
      factor_separate_line = TRUE,
      recalculate_width = TRUE,
      recalculate_height = FALSE,
      merge_models = TRUE
    )
  } else {
    model <- as.formula(paste("su", paste(v_bt, "AGE", "SEX", sep = "+"), sep = "~"))
    m <- coxph(model, data = data1)

    panels <- default_forest_panels(model = m, factor_separate_line = TRUE)
    panels <- append(panels, list(forest_panel(width = 0.08)), after = 8)
    panels[[length(panels)]]$width <- 0.01

    g1 <- forest_model(
      m,
      covariates = c("`%TemB Percentile`", "AGE", "SEX"),
      panels = panels,
      factor_separate_line = TRUE,
      recalculate_width = TRUE,
      recalculate_height = FALSE,
      merge_models = TRUE
    ) +
      theme(
        plot.margin = margin(t = 8, r = 20, b = 8, l = 8),
        panel.spacing.x = unit(1.5, "lines")
      )
  }

  g1 <- g1 + ggtitle(title_text) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  hr_vals <- exp(coef(m))
  hr_min <- min(hr_vals, na.rm = TRUE)
  hr_max <- max(hr_vals, na.rm = TRUE)
  hr_min <- hr_min * 0.95
  hr_max <- hr_max * 1.05
  hr_min <- round(hr_min, 1)
  hr_max <- round(hr_max, 1)
  hr_breaks <- seq(hr_min, hr_max, length.out = 4)
  hr_breaks <- round(hr_breaks, 1)

  g1 <- g1 + scale_x_continuous(breaks = log(hr_breaks), labels = hr_breaks)

  return(g1)
}

DIS <- c(
  "DEATH", "H7_MACULADEGEN", "E4_DM2", "AUTOIMMUNE", "M13_GOUT", "AB1_SEPSIS",
  "N14_RENFAIL", "K11_LIVER", "J10_PNEUMONIA", "FG_HYPERTENSION", "D3_ANAEMIA",
  "G6_ALZHEIMER", "C3_CANCER", "F5_DEMENTIA", "G6_PARKINSON", "M13_OSTEOPOROSIS",
  "C3_LUNG_NONSMALL_EXALLC", "F5_PERSOTH", "E4_DIABINSIPIDUS", "CD2_BENIGN"
)

ensure_dir(cli$output_dir)

for (name in DIS) {
  disease_dir <- file.path(cli$output_dir, name)
  if (dir.exists(disease_dir)) {
    unlink(disease_dir, recursive = TRUE)
  }
  dir.create(disease_dir, recursive = TRUE, showWarnings = FALSE)
}

for (name in DIS) {
  disease_dir <- file.path(cli$output_dir, name)

  p <- KM25_age(data4, "pheb1000", name, 1, 2, "pheb_level", 82, mod = "normal", mod1 = "pap")
  save_plot_pdf(p, file.path(disease_dir, "pheb_40_82_from_0.pdf"), width = 3, height = 3)

  p <- KM25_group(data4, "pheb1000", name, 1, 2, "pheb_level", "40-50")
  save_plot_pdf(p, file.path(disease_dir, "pheb_40_50_KM.pdf"), width = 3, height = 3)

  p <- KM25_group(data4, "pheb1000", name, 1, 2, "pheb_level", "50-60")
  save_plot_pdf(p, file.path(disease_dir, "pheb_50_60_KM.pdf"), width = 3, height = 3)

  p <- KM25_group(data4, "pheb1000", name, 1, 2, "pheb_level", "60-70")
  save_plot_pdf(p, file.path(disease_dir, "pheb_60_70_KM.pdf"), width = 3, height = 3)

  p <- COX(data4, "pheb1000", name, "40-50", 1, 2, "pheb_level_25")
  save_plot_pdf(p, file.path(disease_dir, "pheb_40_50_COX.pdf"), width = 7, height = 4)

  p <- COX(data4, "pheb1000", name, "50-60", 1, 2, "pheb_level_25")
  save_plot_pdf(p, file.path(disease_dir, "pheb_50_60_COX.pdf"), width = 7, height = 4)

  p <- COX(data4, "pheb1000", name, "60-70", 1, 2, "pheb_level_25")
  save_plot_pdf(p, file.path(disease_dir, "pheb_60_70_COX.pdf"), width = 7, height = 4)
}

for (name in DIS) {
  disease_dir <- file.path(cli$output_dir, name)

  p <- KM25_age(data4, "pheb1000", name, 1, 2, "pheb_level", 82, mod = "normal", mod1 = "pap")
  save_plot_png(p, file.path(disease_dir, "pheb_40_82_from_0.png"), width = 3, height = 3)

  p <- KM25_group(data4, "pheb1000", name, 1, 2, "pheb_level", "40-50")
  save_plot_png(p, file.path(disease_dir, "pheb_40_50_KM.png"), width = 3, height = 3)

  p <- KM25_group(data4, "pheb1000", name, 1, 2, "pheb_level", "50-60")
  save_plot_png(p, file.path(disease_dir, "pheb_50_60_KM.png"), width = 3, height = 3)

  p <- KM25_group(data4, "pheb1000", name, 1, 2, "pheb_level", "60-70")
  save_plot_png(p, file.path(disease_dir, "pheb_60_70_KM.png"), width = 3, height = 3)

  p <- COX(data4, "pheb1000", name, "40-50", 1, 2, "pheb_level_25")
  save_plot_png(p, file.path(disease_dir, "pheb_40_50_COX.png"), width = 7, height = 4)

  p <- COX(data4, "pheb1000", name, "50-60", 1, 2, "pheb_level_25")
  save_plot_png(p, file.path(disease_dir, "pheb_50_60_COX.png"), width = 7, height = 4)

  p <- COX(data4, "pheb1000", name, "60-70", 1, 2, "pheb_level_25")
  save_plot_png(p, file.path(disease_dir, "pheb_60_70_COX.png"), width = 7, height = 4)
}
