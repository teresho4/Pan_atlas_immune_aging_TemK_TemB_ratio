args <- commandArgs(trailingOnly = TRUE)

library(data.table)
library(survival)
library(dplyr)
library(tidyr)
library(readxl)
library(openxlsx)
library(ROCR)

`%!in%` <- function(x, y) !(`%in%`(x, y))

print_usage <- function() {
  cat(paste(c(
    "Usage:",
    "  Rscript phenotype_association_clean.R <dataset> <min_age> <adjustment> <predictor> <model> <sex> <label>",
    "",
    "Arguments:",
    "  phenotypes set      finngen or ukbb",
    "  min_age      minimum baseline age to include",
    "  adjustment   no or adj",
    "  predictor    predictor column name, also used to read pred_<predictor>.txt",
    "  model        reg, cox, or coxr where applicable",
    "  sex          all, male, or female",
    "  label        free-text label used in the output file name",
    "",
    "Required files in the current working directory:",
    "  UKBB_FINNGEN_test_condition_excl_proteom.tsv",
    "  pred_<predictor>.txt",
    "  w88907_20250818.csv",
    "  EXCLUDED_CONDITIONS.tsv",
    "",
    "Additional files for dataset = finngen:",
    "  finngen_R11_endpoint_core_noncore_1.0_ordered.tsv",
    "",
    "Additional files for dataset = ukbb:",
    "  participant_1.tsv.gz ... participant_59.tsv.gz",
    "  or participant_1_parsed.tsv.gz ... participant_59_parsed.tsv.gz",
    "  pheno.tsv"
  ), collapse = "\n"), "\n")
}

if (length(args) == 0 || any(args %in% c("-h", "--help"))) {
  print_usage()
  quit(status = 0)
}

if (length(args) < 7) {
  print_usage()
  stop("Expected 7 positional arguments.", call. = FALSE)
}

metadata_file <- "UKBB_FINNGEN_test_condition_excl_proteom.tsv"
prediction_file <- paste0("pred_", args[4], ".txt")
excluded_file <- "w88907_20250818.csv"
excluded_conditions_file <- "EXCLUDED_CONDITIONS.tsv"
finngen_dictionary_file <- "finngen_R11_endpoint_core_noncore_1.0_ordered.tsv"
ukbb_dictionary_file <- "pheno.tsv"

make_output_file <- function() {
  if (args[3] == "no") {
    paste0(args[4], "VS", args[1], "_", args[3], "_", args[2], ".txt")
  } else {
    paste0(args[4], "VS", args[1], "_", args[3], "_", args[2], "_", args[5], "_", args[6], "_", args[7], ".txt")
  }
}

phen_fast_cox_rev <- function(pheno_data, i, X) {
  print(i)
  name <- colnames(pheno_data)[i]
  print(name)

  control <- NA
  case <- NA

  pheno_data[, name] <- as.numeric(pheno_data[, name])

  if ((!grepl("_AGE|_YEAR|_NEVT", name) | name == "BL_AGE") &&
      paste0(name, "_AGE") %in% colnames(pheno_data) &&
      length(table(pheno_data[, name])) > 1 &&
      colnames(pheno_data)[i] %!in% c("eid", "BL_AGE", "SEX")) {
    pheno_data1 <- pheno_data %>%
      dplyr::select(name, paste0(name, "_AGE"), X, "BL_AGE", "SEX") %>%
      drop_na()

    pheno_data1$DEATH_AGE_DIFF <- pheno_data1[, paste0(name, "_AGE")] - pheno_data1$BL_AGE
    pheno_data1 <- pheno_data1 %>% filter(DEATH_AGE_DIFF > 0)

    if (length(table(pheno_data1[, name])) == 2) {
      control <- as.numeric(table(pheno_data1[, name])[1])
      case <- as.numeric(table(pheno_data1[, name])[2])

      su <- Surv(as.numeric(pheno_data1$BL_AGE), as.numeric(pheno_data1[, name]))
      model <- as.formula(paste("su", paste(X, paste0(name, "_AGE"), "SEX", sep = "+"), sep = "~"))

      try({
        m <- coxph(model, data = pheno_data1, na.action = na.exclude)
        p <- summary(m)$coefficients[1, 5]
        print(p)
        b <- summary(m)$coefficients[1, 2]
        se <- summary(m)$coefficients[1, 3]
        ci <- concordance(m)$concordance
      })
    }
  }

  if (!exists("b")) b <- NA
  if (!exists("p")) p <- NA
  if (!exists("ci")) ci <- NA
  if (!exists("se")) se <- NA

  c(names = name, HR = b, pp = p, se = se, ci = ci, control = control, case = case)
}

phen_fast_ukbb <- function(pheno_data, i, X) {
  print(i)
  name <- colnames(pheno_data)[i]
  print(name)

  control <- NA
  case <- NA

  pheno_data[, name] <- as.numeric(pheno_data[, name])

  if (length(table(pheno_data[, name])) > 1 &&
      colnames(pheno_data)[i] %!in% c("eid", "BL_AGE", "SEX") &&
      sum(is.na(pheno_data[, name])) < 0.9 * nrow(pheno_data)) {
    pheno_data1 <- pheno_data %>% dplyr::select(name, X, "BL_AGE", "SEX")

    if (length(table(pheno_data1[, name])) == 2) {
      control <- as.numeric(table(pheno_data1[, name])[1])
      case <- as.numeric(table(pheno_data1[, name])[2])
      model <- as.formula(paste(name, paste(X, "BL_AGE", "SEX", sep = "+"), sep = "~"))

      try({
        m <- glm(model, data = pheno_data1, family = "binomial", na.action = na.exclude)
        p <- summary(m)$coefficients[2, 4]
        b <- summary(m)$coefficients[2, 1]
      })
    } else if (length(table(pheno_data1[, name])) > 2) {
      control <- nrow(pheno_data1 %>% select(name) %>% drop_na())
      case <- NA

      model <- as.formula(paste(name, paste(X, "BL_AGE", "SEX", sep = "+"), sep = "~"))

      try({
        m <- cor.test(pheno_data1[, name], pheno_data1[, X], method = "pearson")
        b_p <- m$estimate
        p_p <- m$p.value

        m <- cor.test(pheno_data1[, name], pheno_data1[, X], method = "spearman")
        b_s <- m$estimate
        p_s <- m$p.value

        m <- lm(model, data = pheno_data1, na.action = na.exclude)
        p <- summary(m)$coefficients[2, 4]
        b <- summary(m)$coefficients[2, 1]
      })
    }
  }

  if (!exists("b")) b <- NA
  if (!exists("p")) p <- NA
  if (!exists("b_p")) b_p <- NA
  if (!exists("p_p")) p_p <- NA
  if (!exists("b_s")) b_s <- NA
  if (!exists("p_s")) p_s <- NA

  c(
    names = name,
    bb = b,
    pp = p,
    r2_pearson = b_p,
    p_pearson = p_p,
    r2_spearman = b_s,
    p_spearman = p_s,
    control = control,
    case = case
  )
}

phen_fast_cox_fish <- function(pheno_data, i, X) {
  print(i)
  name <- colnames(pheno_data)[i]
  print(name)

  control <- NA
  case <- NA
  d00 <- NA
  d01 <- NA
  d10 <- NA
  d11 <- NA

  pheno_data[, name] <- as.numeric(pheno_data[, name])

  if ((!grepl("_AGE|_YEAR|_NEVT", name) | name == "BL_AGE") &&
      paste0(name, "_AGE") %in% colnames(pheno_data) &&
      length(table(pheno_data[, name])) > 1 &&
      colnames(pheno_data)[i] %!in% c("eid", "BL_AGE", "SEX")) {
    pheno_data1 <- pheno_data %>%
      dplyr::select(name, paste0(name, "_AGE"), X, "BL_AGE", "SEX") %>%
      drop_na()

    pheno_data1$DEATH_AGE_DIFF <- pheno_data1[, paste0(name, "_AGE")] - pheno_data1$BL_AGE
    pheno_data1 <- pheno_data1 %>% filter(DEATH_AGE_DIFF > 0)

    if (length(table(pheno_data1[, name])) == 2) {
      control <- as.numeric(table(pheno_data1[, name])[1])
      case <- as.numeric(table(pheno_data1[, name])[2])

      group_columns <- c(name, X)
      d <- pheno_data1 %>%
        group_by(across(all_of(group_columns))) %>%
        dplyr::summarise(N = n(), .groups = "drop")

      if (nrow(d %>% filter(!!sym(name) == 0 & !!sym(X) == 0)) > 0) d00 <- (d %>% filter(!!sym(name) == 0 & !!sym(X) == 0))$N
      if (nrow(d %>% filter(!!sym(name) == 0 & !!sym(X) == 1)) > 0) d01 <- (d %>% filter(!!sym(name) == 0 & !!sym(X) == 1))$N
      if (nrow(d %>% filter(!!sym(name) == 1 & !!sym(X) == 0)) > 0) d10 <- (d %>% filter(!!sym(name) == 1 & !!sym(X) == 0))$N
      if (nrow(d %>% filter(!!sym(name) == 1 & !!sym(X) == 1)) > 0) d11 <- (d %>% filter(!!sym(name) == 1 & !!sym(X) == 1))$N

      su <- Surv(as.numeric(pheno_data1$DEATH_AGE_DIFF), as.numeric(pheno_data1[, name]))
      model <- as.formula(paste("su", paste(X, "BL_AGE", "SEX", sep = "+"), sep = "~"))

      try({
        m <- coxph(model, data = pheno_data1, na.action = na.exclude)
        p <- summary(m)$coefficients[1, 5]
        print(p)
        b <- summary(m)$coefficients[1, 2]
        se <- summary(m)$coefficients[1, 3]
        ci <- concordance(m)$concordance
      })
    }
  }

  if (!exists("b")) b <- NA
  if (!exists("p")) p <- NA
  if (!exists("ci")) ci <- NA
  if (!exists("se")) se <- NA

  c(names = name, HR = b, pp = p, se = se, ci = ci, control = control, case = case, d00 = d00, d01 = d01, d10 = d10, d11 = d11)
}

phen_fast_fish <- function(pheno_data, i, X) {
  print(i)
  name <- colnames(pheno_data)[i]
  print(name)

  control <- NA
  case <- NA
  d00 <- NA
  d01 <- NA
  d10 <- NA
  d11 <- NA

  pheno_data[, name] <- as.numeric(pheno_data[, name])

  if ((!grepl("_AGE|_YEAR|_NEVT", name) | name == "BL_AGE") &&
      paste0(name, "_AGE") %in% colnames(pheno_data) &&
      length(table(pheno_data[, name])) > 1 &&
      colnames(pheno_data)[i] %!in% c("eid", "BL_AGE", "SEX")) {
    pheno_data1 <- pheno_data %>% dplyr::select(name, X, "BL_AGE", "SEX", paste0(name, "_AGE"))
    pheno_data1$DEATH_AGE_DIFF <- pheno_data1[, paste0(name, "_AGE")] - pheno_data1$BL_AGE
    pheno_data1$PHENO <- pheno_data1[, name]

    pheno_data1 <- pheno_data1 %>% filter((DEATH_AGE_DIFF >= 0 & PHENO == 0) | (DEATH_AGE_DIFF <= 0 & PHENO == 1))
    pheno_data1 <- pheno_data1 %>% drop_na()

    pheno_data1$DEATH_AGE_DIFF <- NULL
    pheno_data1[, paste0(name, "_AGE")] <- NULL
    pheno_data1$PHENO <- NULL

    if (length(table(pheno_data1[, name])) == 2) {
      control <- as.numeric(table(pheno_data1[, name])[1])
      case <- as.numeric(table(pheno_data1[, name])[2])

      group_columns <- c(name, X)
      d <- pheno_data1 %>%
        group_by(across(all_of(group_columns))) %>%
        dplyr::summarise(N = n(), .groups = "drop")

      if (nrow(d %>% filter(!!sym(name) == 0 & !!sym(X) == 0)) > 0) d00 <- (d %>% filter(!!sym(name) == 0 & !!sym(X) == 0))$N
      if (nrow(d %>% filter(!!sym(name) == 0 & !!sym(X) == 1)) > 0) d01 <- (d %>% filter(!!sym(name) == 0 & !!sym(X) == 1))$N
      if (nrow(d %>% filter(!!sym(name) == 1 & !!sym(X) == 0)) > 0) d10 <- (d %>% filter(!!sym(name) == 1 & !!sym(X) == 0))$N
      if (nrow(d %>% filter(!!sym(name) == 1 & !!sym(X) == 1)) > 0) d11 <- (d %>% filter(!!sym(name) == 1 & !!sym(X) == 1))$N

      model <- as.formula(paste(name, paste(X, "BL_AGE", "SEX", sep = "+"), sep = "~"))

      try({
        m <- glm(model, data = pheno_data1, family = "binomial", na.action = na.exclude)
        p <- summary(m)$coefficients[2, 4]
        b <- summary(m)$coefficients[2, 1]
      })
    } else if (length(table(pheno_data1[, name])) > 2) {
      control <- nrow(pheno_data1 %>% select(name) %>% drop_na())
      case <- NA

      ranks <- rank(pheno_data1[, name], ties.method = "average")
      pheno_data1[, name] <- qnorm((ranks - 0.5) / nrow(pheno_data1))

      model <- as.formula(paste(name, paste(X, "BL_AGE", "SEX", sep = "+"), sep = "~"))

      try({
        m <- cor.test(pheno_data1[, name], pheno_data1[, X], method = "pearson")
        b_p <- m$estimate
        p_p <- m$p.value

        m <- cor.test(pheno_data1[, name], pheno_data1[, X], method = "spearman")
        b_s <- m$estimate
        p_s <- m$p.value

        m <- lm(model, data = pheno_data1, na.action = na.exclude)
        p <- summary(m)$coefficients[2, 4]
        b <- summary(m)$coefficients[2, 1]
      })
    }
  }

  if (!exists("b")) b <- NA
  if (!exists("p")) p <- NA
  if (!exists("b_p")) b_p <- NA
  if (!exists("p_p")) p_p <- NA
  if (!exists("b_s")) b_s <- NA
  if (!exists("p_s")) p_s <- NA

  c(
    names = name,
    bb = b,
    pp = p,
    r2_pearson = b_p,
    p_pearson = p_p,
    r2_spearman = b_s,
    p_spearman = p_s,
    control = control,
    case = case,
    d00 = d00,
    d01 = d01,
    d10 = d10,
    d11 = d11
  )
}

phen_fast_cox <- function(pheno_data, i, X) {
  print(i)
  name <- colnames(pheno_data)[i]
  print(name)

  control <- NA
  case <- NA

  pheno_data[, name] <- as.numeric(pheno_data[, name])

  if ((!grepl("_AGE|_YEAR|_NEVT", name) | name == "BL_AGE") &&
      paste0(name, "_AGE") %in% colnames(pheno_data) &&
      length(table(pheno_data[, name])) > 1 &&
      colnames(pheno_data)[i] %!in% c("eid", "BL_AGE", "SEX")) {
    pheno_data1 <- pheno_data %>%
      dplyr::select(name, paste0(name, "_AGE"), X, "BL_AGE", "SEX") %>%
      drop_na()

    pheno_data1$DEATH_AGE_DIFF <- pheno_data1[, paste0(name, "_AGE")] - pheno_data1$BL_AGE
    pheno_data1 <- pheno_data1 %>% filter(DEATH_AGE_DIFF > 0)

    if (length(table(pheno_data1[, name])) == 2) {
      control <- as.numeric(table(pheno_data1[, name])[1])
      case <- as.numeric(table(pheno_data1[, name])[2])

      su <- Surv(as.numeric(pheno_data1$DEATH_AGE_DIFF), as.numeric(pheno_data1[, name]))
      model <- as.formula(paste("su", paste(X, "BL_AGE", "SEX", sep = "+"), sep = "~"))

      try({
        m <- coxph(model, data = pheno_data1, na.action = na.exclude)
        p <- summary(m)$coefficients[1, 5]
        print(p)
        b <- summary(m)$coefficients[1, 2]
        se <- summary(m)$coefficients[1, 3]
        ci <- concordance(m)$concordance
      })
    }
  }

  if (!exists("b")) b <- NA
  if (!exists("p")) p <- NA
  if (!exists("ci")) ci <- NA
  if (!exists("se")) se <- NA

  c(names = name, HR = b, pp = p, se = se, ci = ci, control = control, case = case)
}

phen_fast_uni <- function(pheno_data, i, X) {
  print(i)
  name <- colnames(pheno_data)[i]
  print(name)

  control <- NA
  case <- NA

  pheno_data[, name] <- as.numeric(pheno_data[, name])

  if ((!grepl("_AGE|_YEAR|_NEVT", name) | name == "BL_AGE") &&
      paste0(name, "_AGE") %in% colnames(pheno_data) &&
      length(table(pheno_data[, name])) > 1 &&
      colnames(pheno_data)[i] %!in% c("eid")) {
    pheno_data1 <- pheno_data %>% dplyr::select(name, X, paste0(name, "_AGE"))
    pheno_data1$DEATH_AGE_DIFF <- pheno_data1[, paste0(name, "_AGE")] - pheno_data1$BL_AGE
    pheno_data1$PHENO <- pheno_data1[, name]

    pheno_data1 <- pheno_data1 %>% filter((DEATH_AGE_DIFF >= 0 & PHENO == 0) | (DEATH_AGE_DIFF <= 0 & PHENO == 1))
    pheno_data1 <- pheno_data1 %>% drop_na()
    pheno_data1$DEATH_AGE_DIFF <- NULL
    pheno_data1[, paste0(name, "_AGE")] <- NULL
    pheno_data1$PHENO <- NULL

    if (length(table(pheno_data1[, name])) == 2) {
      control <- as.numeric(table(pheno_data1[, name])[1])
      case <- as.numeric(table(pheno_data1[, name])[2])
      model <- as.formula(paste(name, X, sep = "~"))

      try({
        m <- glm(model, data = pheno_data1, family = "binomial", na.action = na.exclude)
        p <- summary(m)$coefficients[2, 4]
        b <- summary(m)$coefficients[2, 1]
      })
    } else if (length(table(pheno_data1[, name])) > 2) {
      control <- nrow(pheno_data1 %>% select(name) %>% drop_na())
      case <- NA

      ranks <- rank(pheno_data1[, name], ties.method = "average")
      pheno_data1[, name] <- qnorm((ranks - 0.5) / nrow(pheno_data1))

      model <- as.formula(paste(name, X, sep = "~"))

      try({
        m <- cor.test(pheno_data1[, name], pheno_data1[, X], method = "pearson")
        b_p <- m$estimate
        p_p <- m$p.value

        m <- cor.test(pheno_data1[, name], pheno_data1[, X], method = "spearman")
        b_s <- m$estimate
        p_s <- m$p.value

        m <- lm(model, data = pheno_data1, na.action = na.exclude)
        p <- summary(m)$coefficients[2, 4]
        b <- summary(m)$coefficients[2, 1]
      })
    }
  }

  if (!exists("b")) b <- NA
  if (!exists("p")) p <- NA
  if (!exists("b_p")) b_p <- NA
  if (!exists("p_p")) p_p <- NA
  if (!exists("b_s")) b_s <- NA
  if (!exists("p_s")) p_s <- NA

  c(
    names = name,
    bb = b,
    pp = p,
    r2_pearson = b_p,
    p_pearson = p_p,
    r2_spearman = b_s,
    p_spearman = p_s,
    control = control,
    case = case
  )
}

phen_fast <- function(pheno_data, i, X) {
  print(i)
  name <- colnames(pheno_data)[i]
  print(name)

  control <- NA
  case <- NA

  pheno_data[, name] <- as.numeric(pheno_data[, name])

  if ((!grepl("_AGE|_YEAR|_NEVT", name) | name == "BL_AGE") &&
      paste0(name, "_AGE") %in% colnames(pheno_data) &&
      length(table(pheno_data[, name])) > 1 &&
      colnames(pheno_data)[i] %!in% c("eid", "BL_AGE", "SEX")) {
    pheno_data1 <- pheno_data %>% dplyr::select(name, X, "BL_AGE", "SEX", paste0(name, "_AGE"))
    pheno_data1$DEATH_AGE_DIFF <- pheno_data1[, paste0(name, "_AGE")] - pheno_data1$BL_AGE
    pheno_data1$PHENO <- pheno_data1[, name]

    pheno_data1 <- pheno_data1 %>% filter((DEATH_AGE_DIFF >= 0 & PHENO == 0) | (DEATH_AGE_DIFF <= 0 & PHENO == 1))
    pheno_data1 <- pheno_data1 %>% drop_na()

    pheno_data1$DEATH_AGE_DIFF <- NULL
    pheno_data1[, paste0(name, "_AGE")] <- NULL
    pheno_data1$PHENO <- NULL

    if (length(table(pheno_data1[, name])) == 2) {
      control <- as.numeric(table(pheno_data1[, name])[1])
      case <- as.numeric(table(pheno_data1[, name])[2])
      model <- as.formula(paste(name, paste(X, "BL_AGE", "SEX", sep = "+"), sep = "~"))

      try({
        m <- glm(model, data = pheno_data1, family = "binomial", na.action = na.exclude)
        p <- summary(m)$coefficients[2, 4]
        b <- summary(m)$coefficients[2, 1]
      })
    } else if (length(table(pheno_data1[, name])) > 2) {
      control <- nrow(pheno_data1 %>% select(name) %>% drop_na())
      case <- NA

      ranks <- rank(pheno_data1[, name], ties.method = "average")
      pheno_data1[, name] <- qnorm((ranks - 0.5) / nrow(pheno_data1))

      model <- as.formula(paste(name, paste(X, "BL_AGE", "SEX", sep = "+"), sep = "~"))

      try({
        m <- cor.test(pheno_data1[, name], pheno_data1[, X], method = "pearson")
        b_p <- m$estimate
        p_p <- m$p.value

        m <- cor.test(pheno_data1[, name], pheno_data1[, X], method = "spearman")
        b_s <- m$estimate
        p_s <- m$p.value

        m <- lm(model, data = pheno_data1, na.action = na.exclude)
        p <- summary(m)$coefficients[2, 4]
        b <- summary(m)$coefficients[2, 1]
      })
    }
  }

  if (!exists("b")) b <- NA
  if (!exists("p")) p <- NA
  if (!exists("b_p")) b_p <- NA
  if (!exists("p_p")) p_p <- NA
  if (!exists("b_s")) b_s <- NA
  if (!exists("p_s")) p_s <- NA

  c(
    names = name,
    bb = b,
    pp = p,
    r2_pearson = b_p,
    p_pearson = p_p,
    r2_spearman = b_s,
    p_spearman = p_s,
    control = control,
    case = case
  )
}

if (args[1] == "finngen") {
  data3 <- as.data.frame(fread(metadata_file, header = TRUE))
  pheb1 <- as.data.frame(fread(prediction_file, header = TRUE, sep = "\t"))
  data3 <- as.data.frame(merge(pheb1, data3, by = "eid"))

  EXCLUDED <- read.table(excluded_file, sep = "\t", header = FALSE)
  colnames(EXCLUDED)[1] <- "eid"
  data3 <- data3 %>% filter(eid %!in% EXCLUDED$eid)

  EXCLUDED_CONDITIONS <- read.table(excluded_conditions_file, sep = "\t", header = TRUE)
  data3 <- data3[, !(names(data3) %in% EXCLUDED_CONDITIONS$x)]

  data3 <- data3 %>% filter(BL_AGE >= as.numeric(args[2]))

  if (args[6] == "male") {
    data3 <- data3 %>% filter(SEX == 1)
  }

  if (args[6] == "female") {
    data3 <- data3 %>% filter(SEX == 2)
  }

  if (args[3] == "no") {
    columns <- 1:ncol(data3)
    res <- sapply(columns, function(column) phen_fast_uni(data3, column, args[4]))
    res <- as.data.frame(t(res))
    res$pp <- as.numeric(res$pp)
    res$control <- as.numeric(res$control)
    res$case <- as.numeric(res$case)
    res <- res %>% dplyr::arrange(pp) %>% filter(!is.na(pp))

    dict <- as.data.frame(fread(finngen_dictionary_file, header = TRUE))
    dict <- dict %>% select(NAME, LONGNAME)
    colnames(dict) <- c("names", "FEATURE")
    res <- as.data.frame(merge(res, dict, by = "names", all.x = TRUE))
    res <- res %>% dplyr::arrange(pp)

    write.table(res, make_output_file(), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  }

  if (args[3] == "adj") {
    columns <- 1:ncol(data3)

    if (args[5] == "reg" && args[4] != "p30000_cat_i0") {
      res <- sapply(columns, function(column) phen_fast(data3, column, args[4]))
    } else if (args[5] == "reg" && args[4] == "p30000_cat_i0") {
      res <- sapply(columns, function(column) phen_fast_fish(data3, column, args[4]))
    } else if (args[5] == "cox" && args[4] != "p30000_cat_i0") {
      res <- sapply(columns, function(column) phen_fast_cox(data3, column, args[4]))
    } else if (args[5] == "coxr" && args[4] != "p30000_cat_i0") {
      res <- sapply(columns, function(column) phen_fast_cox_rev(data3, column, args[4]))
    } else if (args[5] == "cox" && args[4] == "p30000_cat_i0") {
      res <- sapply(columns, function(column) phen_fast_cox_fish(data3, column, args[4]))
    }

    res <- as.data.frame(t(res))
    res$pp <- as.numeric(res$pp)
    res$control <- as.numeric(res$control)
    res$case <- as.numeric(res$case)
    res <- res %>% dplyr::arrange(pp) %>% filter(!is.na(pp))

    dict <- as.data.frame(fread(finngen_dictionary_file, header = TRUE))
    dict <- dict %>% select(NAME, LONGNAME)
    colnames(dict) <- c("names", "FEATURE")
    res <- as.data.frame(merge(res, dict, by = "names", all.x = TRUE))
    res <- res %>% dplyr::arrange(pp)

    write.table(res, make_output_file(), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  }
}

if (args[1] == "ukbb") {
  data3 <- as.data.frame(fread(metadata_file, header = TRUE, select = c("eid", "SEX", "BL_AGE")))
  pheb1 <- as.data.frame(fread(prediction_file, header = TRUE, sep = "\t"))
  data3 <- as.data.frame(merge(pheb1, data3, by = "eid"))

  data3 <- data3 %>% mutate(!!args[4] := case_when(!!sym(args[4]) < 0 ~ 0, TRUE ~ !!sym(args[4])))

  EXCLUDED <- read.table(excluded_file, sep = "\t", header = TRUE)
  colnames(EXCLUDED)[1] <- "eid"
  data3 <- data3 %>% filter(eid %!in% EXCLUDED$eid)

  EXCLUDED_CONDITIONS <- read.table(excluded_conditions_file, sep = "\t", header = TRUE)
  data3 <- data3[, !(names(data3) %in% EXCLUDED_CONDITIONS$x)]

  data3 <- data3 %>% filter(BL_AGE >= as.numeric(args[2]))

  if (args[6] == "male") {
    data3 <- data3 %>% filter(SEX == 1)
  }

  if (args[6] == "female") {
    data3 <- data3 %>% filter(SEX == 2)
  }

  RES <- data.frame(matrix(ncol = 5, nrow = 0))
  colnames(RES) <- c("names", "bb", "sese", "pp", "cici")

  for (i in 1:59) {
    parsed_file <- paste0("participant_", i, "_parsed.tsv.gz")
    raw_file <- paste0("participant_", i, ".tsv.gz")

    if (file.exists(parsed_file)) {
      pheno_data <- fread(parsed_file, header = TRUE, sep = "\t")
    } else {
      pheno_data <- fread(raw_file, header = TRUE, sep = "\t")
    }

    pheno_data <- as.data.frame(pheno_data)
    pheno_data <- pheno_data[, !colSums(is.na(pheno_data)) == nrow(pheno_data)]

    if ("p21003_i0" %in% colnames(pheno_data)) {
      pheno_data$p21003_i0 <- NULL
    }

    pheno_data <- as.data.frame(merge(data3, pheno_data, by = "eid", all.x = TRUE))
    columns <- 1:ncol(pheno_data)

    if (args[3] == "no") {
      res <- sapply(columns, function(column) phen_fast_uni(pheno_data, column, args[4]))
    }

    if (args[3] == "adj") {
      res <- sapply(columns, function(column) phen_fast_ukbb(pheno_data, column, args[4]))
    }

    res <- as.data.frame(t(res))
    RES <- rbind(RES, res)
  }

  pheno <- as.data.frame(fread(ukbb_dictionary_file, header = TRUE))
  pheno <- pheno %>% select(pheno, description, linkout)
  colnames(pheno)[1] <- "names"

  RES <- RES %>% filter(!is.na(pp)) %>% dplyr::arrange(pp) %>% filter(names %!in% c(args[4], "SEX", "BL_AGE"))
  RES <- as.data.frame(merge(RES, pheno, by = "names", all.x = TRUE))
  RES$pp <- as.numeric(RES$pp)
  RES <- RES %>% dplyr::arrange(pp)

  write.table(RES, make_output_file(), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}
