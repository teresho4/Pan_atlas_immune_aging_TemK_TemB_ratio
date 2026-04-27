library(data.table)
library(dplyr)
library(tidyr)
library(glmnet)
library(readxl)
library(ggplot2)
library(caret)
library(sva)
library(scales)
library(gridExtra)

cli_defaults <- list(
  seattle_pct = "",
  seattle_olink = "",
  stl_annotations = "",
  stl_olink = "",
  kcl_annotations = "",
  kcl_olink = "",
  ra_annotations = "",
  ra_olink = "",
  ukbb_metadata = "",
  ukbb_excluded = "",
  ukbb_olink_file = "",
  ukbb_olink_instance_1 = "",
  ukbb_olink_instance_2 = "",
  ukbb_olink_instance_3 = "",
  ukbb_olink_instance_4 = "",
  ukbb_olink_instance_5 = "",
  ukbb_olink_instance_6 = "",
  out_union_cor1_f = "outputs/UNION_cor1_f.txt",
  out_ukbb_predictions = "outputs/ukbb_pheb_predictions.txt",
  out_r2_metrics = "outputs/model_r2_metrics.txt",
  out_coefficients = "outputs/model_coefficients.txt",
  out_cross_validation = "outputs/cross_validation_r2.txt",
  out_lasso_model = "outputs/protein_lasso_model.rds",
  out_pca_plot = "outputs/pca_plot.pdf",
  out_r2_plot = "outputs/r2_plot.pdf",
  out_crossval_plot = "outputs/cross_validation_r2_plot.pdf"
)

cli_required <- c(
  "seattle_pct",
  "seattle_olink",
  "stl_annotations",
  "stl_olink",
  "kcl_annotations",
  "kcl_olink",
  "ra_annotations",
  "ra_olink"
)

ukbb_metadata_args <- c(
  "ukbb_metadata",
  "ukbb_excluded"
)

ukbb_single_file_args <- c(
  "ukbb_olink_file"
)

ukbb_instance_args <- c(
  "ukbb_olink_instance_1",
  "ukbb_olink_instance_2",
  "ukbb_olink_instance_3",
  "ukbb_olink_instance_4",
  "ukbb_olink_instance_5",
  "ukbb_olink_instance_6"
)

ukbb_input_args <- c(
  ukbb_metadata_args,
  ukbb_single_file_args,
  ukbb_instance_args
)

print_usage <- function(defaults, required) {
  required_lines <- paste0("  --", gsub("_", "-", required), " <path>")
  ukbb_metadata_lines <- paste0("  --", gsub("_", "-", ukbb_metadata_args), " <path>")
  ukbb_single_lines <- paste0("  --", gsub("_", "-", ukbb_single_file_args), " <path>")
  ukbb_instance_lines <- paste0("  --", gsub("_", "-", ukbb_instance_args), " <path>")
  output_args <- grep("^out_", names(defaults), value = TRUE)
  output_lines <- paste0(
    "  --", gsub("_", "-", output_args), " <path>",
    "    default: ", unlist(defaults[output_args], use.names = FALSE)
  )
  cat(paste(c(
    "Usage:",
    "  Rscript protein_model_public_ukbb_1_or_6.R [arguments]",
    "",
    "Required input files:",
    required_lines,
    "",
    "Optional UKBB metadata files required for either UKBB proteome mode:",
    ukbb_metadata_lines,
    "",
    "Optional UKBB proteome input, mode 1: one pre-merged proteome file:",
    ukbb_single_lines,
    "",
    "Optional UKBB proteome input, mode 2: six Olink instance files:",
    ukbb_instance_lines,
    "",
    "Use either --ukbb-olink-file or all six --ukbb-olink-instance-* arguments. Do not use both modes together.",
    "",
    "Optional output files:",
    output_lines,
    "",
    "Arguments can be passed as --name value or --name=value. Hyphens and underscores are equivalent in argument names."
  ), collapse = "\n"), "\n")
}

parse_cli_args <- function(defaults, required) {
  raw_args <- commandArgs(trailingOnly = TRUE)

  if (any(raw_args %in% c("-h", "--help"))) {
    print_usage(defaults, required)
    quit(status = 0)
  }

  if (length(raw_args) == 0) {
    print_usage(defaults, required)
    stop("Missing required arguments.", call. = FALSE)
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

  missing <- required[!nzchar(unlist(values[required], use.names = FALSE))]
  if (length(missing) > 0) {
    stop(paste0("Missing required arguments: ", paste(paste0("--", gsub("_", "-", missing)), collapse = ", ")), call. = FALSE)
  }

  values
}

ensure_parent_dir <- function(path) {
  parent <- dirname(path)
  if (!dir.exists(parent)) {
    dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  }
}

write_tsv <- function(x, path, row_names = FALSE) {
  ensure_parent_dir(path)
  write.table(x, path, sep = "\t", quote = FALSE, row.names = row_names, col.names = TRUE)
}

coalesce_xy_columns <- function(df, columns) {
  for (column in columns) {
    x_col <- paste0(column, ".x")
    y_col <- paste0(column, ".y")

    if (x_col %in% names(df) && y_col %in% names(df)) {
      df[[column]] <- dplyr::coalesce(df[[x_col]], df[[y_col]])
      df[[x_col]] <- NULL
      df[[y_col]] <- NULL
    } else if (x_col %in% names(df) && !column %in% names(df)) {
      names(df)[names(df) == x_col] <- column
    } else if (y_col %in% names(df) && !column %in% names(df)) {
      names(df)[names(df) == y_col] <- column
    }
  }

  df
}

safe_r2 <- function(pred, obs) {
  ok <- is.finite(pred) & is.finite(obs)
  if (sum(ok) < 2) {
    return(NA_real_)
  }
  round(caret::R2(pred[ok], obs[ok]), 2)
}


cli <- parse_cli_args(cli_defaults, cli_required)

ukbb_metadata_values <- unlist(cli[ukbb_metadata_args], use.names = FALSE)
ukbb_instance_values <- unlist(cli[ukbb_instance_args], use.names = FALSE)

has_any_ukbb_metadata <- any(nzchar(ukbb_metadata_values))
has_all_ukbb_metadata <- all(nzchar(ukbb_metadata_values))
has_single_ukbb_file <- nzchar(cli$ukbb_olink_file)
has_any_ukbb_instances <- any(nzchar(ukbb_instance_values))
has_all_ukbb_instances <- all(nzchar(ukbb_instance_values))

if (has_single_ukbb_file && has_any_ukbb_instances) {
  stop(
    "Provide either --ukbb-olink-file or the six --ukbb-olink-instance-* files, not both.",
    call. = FALSE
  )
}

if (has_any_ukbb_instances && !has_all_ukbb_instances) {
  missing_ukbb_instances <- ukbb_instance_args[!nzchar(ukbb_instance_values)]
  stop(
    paste0(
      "The six-file UKBB proteome mode requires all six instance files. Missing: ",
      paste(paste0("--", gsub("_", "-", missing_ukbb_instances)), collapse = ", ")
    ),
    call. = FALSE
  )
}

has_ukbb_proteome <- has_single_ukbb_file || has_all_ukbb_instances
has_any_ukbb <- has_any_ukbb_metadata || has_ukbb_proteome || has_any_ukbb_instances

if (has_any_ukbb && !has_all_ukbb_metadata) {
  missing_ukbb_metadata <- ukbb_metadata_args[!nzchar(ukbb_metadata_values)]
  stop(
    paste0(
      "UKBB metadata files are required when UKBB proteome input is provided. Missing: ",
      paste(paste0("--", gsub("_", "-", missing_ukbb_metadata)), collapse = ", ")
    ),
    call. = FALSE
  )
}

if (has_all_ukbb_metadata && !has_ukbb_proteome) {
  stop(
    paste0(
      "UKBB metadata were provided, but no UKBB proteome input was provided. ",
      "Use either --ukbb-olink-file or all six --ukbb-olink-instance-* arguments."
    ),
    call. = FALSE
  )
}

has_ukbb <- has_all_ukbb_metadata && has_ukbb_proteome
ukbb_proteome_mode <- if (has_single_ukbb_file) {
  "single_file"
} else if (has_all_ukbb_instances) {
  "six_instances"
} else {
  "none"
}

if (has_ukbb) {
  message("Using UKBB proteome input mode: ", ukbb_proteome_mode, ".")
} else {
  message("UKBB inputs were not provided; skipping UKBB processing and UKBB prediction output.")
}


.r2_hist <- function(df, var, title,
                     fill = "#D55E00",
                     line_col = "#222222",
                     sd_col = "grey50") {
  v <- rlang::sym(var)
  m <- mean(df[[var]], na.rm = TRUE)
  sd_ <- sd(df[[var]], na.rm = TRUE)

  ggplot(df, aes(x = !!v)) +
    geom_histogram(
      bins = 30,
      fill = fill,
      colour = "white",
      linewidth = 0.3
    ) +
    geom_vline(
      xintercept = m,
      linetype = "dashed",
      colour = line_col
    ) +
    geom_vline(xintercept = m - sd_, linetype = "dotted", colour = sd_col) +
    geom_vline(xintercept = m + sd_, linetype = "dotted", colour = sd_col) +
    annotate(
      "text",
      x = m,
      y = Inf,
      vjust = 1.4,
      hjust = -0.05,
      label = sprintf("mean = %.3f\nSD = %.3f", m, sd_),
      colour = line_col,
      size = 3.5
    ) +
    scale_x_continuous(
      labels = number_format(accuracy = 0.01),
      limits = c(0, 0.6),
      expand = c(0, 0)
    ) +
    labs(
      x = expression(R^2),
      y = "Count",
      title = title
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title.position = "plot",
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )
}

get_cross_val_plot <- function(d) {
  g1 <- .r2_hist(d, "R2_all", "ALL cohorts")
  g2 <- .r2_hist(d, "R2_stl", "ABF300 test", fill = "#56B4E9")
  g3 <- .r2_hist(d, "R2_sea", "Sound Life", fill = "#E69F00")
  g4 <- .r2_hist(d, "R2_ra", "AIFI RA test", fill = "#999999")
  g5 <- .r2_hist(d, "R2_kcl", "KCL test", fill = "#234567")

  patchwork::wrap_plots(
    g1, g2, g3,
    g4, g5, patchwork::plot_spacer(),
    ncol = 3, byrow = TRUE
  ) +
    patchwork::plot_annotation(
      title = "Distribution of cross-validated R\u00b2 values (mean \u00b1 SD)",
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(size = 14, face = "bold")
      )
    )
}

cross_valid <- function(UNION, train, train1, train2, train3,
                        test, test1, mod, k, proc,
                        alpha_grid = seq(0, 1, by = 0.1),
                        pheno) {
  DF <- data.frame(matrix(ncol = 8, nrow = 0))
  colnames(DF) <- c("alpha", "lambda", "R2_all", "R2_stl", "R2_sea", "R2_kcl", "R2_sea_stl", "R2_ra")

  for (a in alpha_grid) {
    for (i in 1:k) {
      print(paste("alpha =", a, "fold =", i))

      f <- UNION %>% filter(COHORT == test | COHORT == test1)

      ff <- UNION %>% filter(COHORT == train | COHORT == train1 | COHORT == train2 | COHORT == train3)
      if (test == "RAH") {
        ff <- ff %>% filter(ID %!in% RA_unhealthy$ID)
      }
      UNION_val <- ff %>% sample_n(nrow(ff) * proc)

      SEA <- nrow(UNION_val %>% filter(COHORT == "SEA"))
      STL <- nrow(UNION_val %>% filter(COHORT == "STL"))
      KCL <- nrow(UNION_val %>% filter(COHORT == "KCL"))

      while (SEA < 15 & STL < 17 & train == "SEA" & train1 == "STL") {
        UNION_val <- ff %>% sample_n(nrow(ff) * proc)
        SEA <- nrow(UNION_val %>% filter(COHORT == "SEA"))
        STL <- nrow(UNION_val %>% filter(COHORT == "STL"))
        KCL <- nrow(UNION_val %>% filter(COHORT == "KCL"))
      }

      COHORT <- UNION_val$COHORT
      ID <- UNION_val$ID
      yy <- UNION_val[[pheno]]

      UNION_train <- ff %>% filter(ID %!in% UNION_val$ID)

      UNION_val <- UNION_val %>%
        select(-c("pheb", "phek", "ID", "COHORT")) %>%
        select_if(~ !any(is.na(.)))

      y <- UNION_train$pheb
      UNION_train <- UNION_train %>%
        select(-c("pheb", "phek", "ID", "COHORT")) %>%
        select_if(~ !any(is.na(.)))

      lasso_model <- cv.glmnet(as.matrix(UNION_train), y, alpha = a, nfolds = 10)
      if (mod == "min") {
        best_lambda <- lasso_model$lambda.min
      } else {
        best_lambda <- lasso_model$lambda.1se
      }

      final_model <- glmnet(as.matrix(UNION_train), y, alpha = a, lambda = best_lambda)

      selected_features <- which(coef(final_model) != 0)
      selected_features_names <- rownames(coef(final_model))[selected_features]

      pred1 <- as.numeric(predict(final_model, as.matrix(UNION_val)))

      coef_df <- as.data.frame(as.matrix(coef(final_model)))
      coef_df <- subset(coef_df, `s0` != 0)

      val <- cbind(UNION_val, pred1)
      val <- cbind(val, yy)
      val <- cbind(val, COHORT)
      val <- cbind(val, ID)

      val_stl <- val %>% filter(COHORT == "STL")
      val_sea <- val %>% filter(COHORT == "SEA")
      val_sea_stl <- val %>% filter(COHORT == "SEA" | COHORT == "STL")
      val_ra <- val %>% filter(COHORT == "RA")
      val_kcl <- val %>% filter(COHORT == "KCL")

      r2 <- safe_r2(val$pred1, val$yy)
      r2_stl <- safe_r2(val_stl$pred1, val_stl$yy)
      r2_sea <- safe_r2(val_sea$pred1, val_sea$yy)
      r2_sea_stl <- safe_r2(val_sea_stl$pred1, val_sea_stl$yy)
      r2_ra <- safe_r2(val_ra$pred1, val_ra$yy)
      r2_kcl <- safe_r2(val_kcl$pred1, val_kcl$yy)

      df <- data.frame(
        alpha = a,
        lambda = best_lambda,
        R2_all = r2,
        R2_stl = r2_stl,
        R2_sea = r2_sea,
        R2_kcl = r2_kcl,
        R2_sea_stl = r2_sea_stl,
        R2_ra = r2_ra
      )

      DF <- rbind(DF, df)
    }
  }

  return(DF)
}

get_R2_plots <- function(pred_lasso, sets = "r2_stl_sea", train, pheno) {
  d <- pred_lasso[[2]]

  g1 <- d %>%
    filter(COHORT == "STL") %>%
    ggplot(aes(yy, pred1, col = COHORT)) +
    geom_point() +
    geom_smooth(method = "lm") +
    theme_classic() +
    ylab(paste0("Predicted ", pheno)) +
    xlab(paste0(pheno)) +
    labs(title = "STL") +
    annotate(
      "text",
      x = mean(d[d$COHORT == "STL", ]$pred1, na.rm = T),
      y = max(d[d$COHORT == "STL", ]$pred1, na.rm = T) * 1,
      label = paste("R2: ", pred_lasso[[3]][["r2_stl"]], sep = ""),
      color = "red"
    )

  g2 <- d %>%
    filter(COHORT == "SEA") %>%
    ggplot(aes(yy, pred1, col = COHORT)) +
    geom_point() +
    geom_smooth(method = "lm") +
    theme_classic() +
    ylab(paste0("Predicted ", pheno)) +
    xlab(paste0(pheno)) +
    labs(title = "SEA") +
    annotate(
      "text",
      x = mean(d[d$COHORT == "SEA", ]$pred1, na.rm = T),
      y = max(d[d$COHORT == "SEA", ]$pred1, na.rm = T) * 1,
      label = paste("R2: ", pred_lasso[[3]][["r2_sea"]], sep = ""),
      color = "red"
    )

  g4 <- d %>%
    filter(COHORT == "RA") %>%
    filter(ID %in% RA_healthy$ID) %>%
    ggplot(aes(yy, pred1, col = COHORT)) +
    geom_point() +
    geom_smooth(method = "lm") +
    theme_classic() +
    ylab(paste0("Predicted ", pheno)) +
    xlab(paste0(pheno)) +
    labs(title = "RA healthy") +
    annotate(
      "text",
      x = mean(d[d$COHORT == "RA", ]$pred1, na.rm = T),
      y = max(d[d$COHORT == "RA", ]$pred1, na.rm = T) * 1,
      label = paste("R2: ", pred_lasso[[3]][["r2_rah"]], sep = ""),
      color = "red"
    )

  g6 <- d %>%
    filter(COHORT == "KCL") %>%
    ggplot(aes(yy, pred1, col = COHORT)) +
    geom_point() +
    geom_smooth(method = "lm") +
    theme_classic() +
    ylab(paste0("Predicted ", pheno)) +
    xlab(paste0(pheno)) +
    labs(title = "KCL") +
    annotate(
      "text",
      x = mean(d[d$COHORT == "KCL", ]$pred1, na.rm = T),
      y = max(d[d$COHORT == "KCL", ]$pred1, na.rm = T) * 1,
      label = paste("R2: ", pred_lasso[[3]][["r2_kcl"]], sep = ""),
      color = "red"
    )

  if (sets == "r2") {

    g7 <- d %>%
      filter(COHORT %in% train) %>%
      ggplot(aes(yy, pred1)) +
      geom_point() +
      geom_smooth(method = "lm") +
      theme_classic() +
      ylab(paste0("Predicted ", pheno)) +
      xlab(paste0(pheno)) +
      labs(title = "all training cohorts") +
      annotate(
        "text",
        x = mean(d[d$COHORT == "STL", ]$pred1, na.rm = T),
        y = max(d[d$COHORT == "STL", ]$pred1, na.rm = T) * 1,
        label = paste("R2: ", pred_lasso[[3]][[sets]], sep = ""),
        color = "red"
      )
  } else {
    g7 <- d %>%
      filter(COHORT %in% c("SEA", "STL", "RA")) %>%
      ggplot(aes(yy, pred1)) +
      geom_point() +
      geom_smooth(method = "lm") +
      theme_classic() +
      ylab(paste0("Predicted ", pheno)) +
      xlab(paste0(pheno)) +
      labs(title = "all training cohorts") +
      annotate(
        "text",
        x = mean(d[d$COHORT == "STL", ]$pred1, na.rm = T),
        y = max(d[d$COHORT == "STL", ]$pred1, na.rm = T) * 1,
        label = paste("R2: ", pred_lasso[[3]][[sets]], sep = ""),
        color = "red"
      )
  }

  gridExtra::grid.arrange(g1, g2, g4, g6, g7, ncol = 3, nrow = 2)
}

get_lasso_model <- function(UNION_cor1_f, train_, N, pheno,
                            alpha_grid = seq(0, 1, by = 0.1),
                            pick_metric = "r2",
                            do_scale = FALSE) {
  ff <- UNION_cor1_f %>% filter(COHORT %in% train_)
  if ("RAH" %in% train_) {
    ff <- UNION_cor1_f %>% filter(COHORT %in% c(train_, "RA"))
    ff <- ff %>% filter(ID %!in% RA_unhealthy$ID)
  }

  y <- ff[[pheno]]
  yy <- UNION_cor1_f[[pheno]]
  COHORT <- UNION_cor1_f$COHORT
  ID <- UNION_cor1_f$ID

  UNION_val <- UNION_cor1_f %>%
    select(-c("pheb", "phek", "ID", "COHORT")) %>%
    select_if(~ !any(is.na(.)))

  UNION_train <- ff %>%
    select(-c("pheb", "phek", "ID", "COHORT")) %>%
    select_if(~ !any(is.na(.)))

  scale_params <- NULL

  UNION_train_mat <- as.matrix(UNION_train)
  UNION_val_mat <- as.matrix(UNION_val)

  if (do_scale) {
    center <- colMeans(UNION_train_mat)
    sc <- apply(UNION_train_mat, 2, sd)

    sc[!is.finite(sc) | sc == 0] <- 1

    UNION_train_mat <- sweep(UNION_train_mat, 2, center, FUN = "-")
    UNION_train_mat <- sweep(UNION_train_mat, 2, sc, FUN = "/")

    UNION_val_mat <- sweep(UNION_val_mat, 2, center, FUN = "-")
    UNION_val_mat <- sweep(UNION_val_mat, 2, sc, FUN = "/")

    scale_params <- list(center = center, scale = sc)
  }

  best_alpha <- NA_real_
  best_metric <- -Inf
  best_lambda <- NA_real_
  best_final_model <- NULL

  all_alpha_results <- list()

  for (a in alpha_grid) {
    lasso_model <- cv.glmnet(
      x = UNION_train_mat, y = y,
      alpha = a, nfolds = 10,
      standardize = !do_scale
    )

    best_lambda_a <- lasso_model$lambda.min
    final_model_a <- glmnet(
      x = UNION_train_mat, y = y,
      alpha = a, lambda = best_lambda_a,
      standardize = !do_scale
    )

    pred1_a <- as.numeric(predict(final_model_a, UNION_val_mat))

    val_a <- cbind(UNION_val, pred1 = pred1_a)
    val_a <- cbind(val_a, yy = yy)
    val_a <- cbind(val_a, COHORT = COHORT)
    val_a <- cbind(val_a, ID = ID)

    val_ra <- val_a %>% filter(COHORT == "RA")
    val_stl <- val_a %>% filter(COHORT == "STL")
    val_sea <- val_a %>% filter(COHORT == "SEA")
    val_kcl <- val_a %>% filter(COHORT == "KCL")
    val_stl_sea <- val_a %>% filter(COHORT %in% c("SEA", "STL"))
    val_rah <- val_a %>% filter(COHORT == "RA") %>% filter(ID %in% RA_healthy$ID)
    val_rauh <- val_a %>% filter(COHORT == "RA") %>% filter(ID %in% RA_unhealthy$ID)
    val_all <- val_a %>% filter(COHORT %in% c("SEA", "STL", "RA", "KCL"))
    val_stl1 <- val_a %>% filter(COHORT %in% c("STL1"))

    r2 <- safe_r2(val_all$pred1, val_all$yy)
    r2_rauh <- safe_r2(val_rauh$pred1, val_rauh$yy)
    r2_rah <- safe_r2(val_rah$pred1, val_rah$yy)
    r2_ra <- safe_r2(val_ra$pred1, val_ra$yy)
    r2_stl <- safe_r2(val_stl$pred1, val_stl$yy)
    r2_sea <- safe_r2(val_sea$pred1, val_sea$yy)
    r2_kcl <- safe_r2(val_kcl$pred1, val_kcl$yy)
    r2_stl_sea <- safe_r2(val_stl_sea$pred1, val_stl_sea$yy)
    r2_stl1 <- safe_r2(val_stl1$pred1, val_stl1$yy)

    r2_ <- c(r2_rah, r2_rauh, r2_ra, r2_stl, r2_sea, r2_kcl, r2_stl_sea, r2, r2_stl1)
    names(r2_) <- c("r2_rah", "r2_rauh", "r2_ra", "r2_stl", "r2_sea", "r2_kcl", "r2_stl_sea", "r2", "r2_stl1")

    metric_val <- unname(r2_[[pick_metric]])
    all_alpha_results[[paste0("alpha_", a)]] <- list(alpha = a, lambda = best_lambda_a, r2_ = r2_)

    if (!is.na(metric_val) && metric_val > best_metric) {
      best_metric <- metric_val
      best_alpha <- a
      best_lambda <- best_lambda_a
      best_final_model <- final_model_a
    }
  }

  final_model <- best_final_model
  cat("Best alpha:", best_alpha, "Best lambda:", best_lambda, "Best", pick_metric, "=", best_metric, "\n")

  selected_features <- which(coef(final_model) != 0)
  selected_features_names <- rownames(coef(final_model))[selected_features]


  pred1 <- as.numeric(predict(final_model, UNION_val_mat))

  coef_df <- as.data.frame(as.matrix(coef(final_model)))
  coef_df <- subset(coef_df, `s0` != 0)


  coef_dff <- coef_df %>% mutate(s0 = abs(s0))

  val <- cbind(UNION_val, pred1)
  val <- cbind(val, yy)
  val <- cbind(val, COHORT)
  val <- cbind(val, ID)

  val_ra <- val %>% filter(COHORT == "RA")
  val_stl <- val %>% filter(COHORT == "STL")
  val_sea <- val %>% filter(COHORT == "SEA")
  val_kcl <- val %>% filter(COHORT == "KCL")
  val_stl_sea <- val %>% filter(COHORT %in% c("SEA", "STL"))
  val_rah <- val %>% filter(COHORT == "RA") %>% filter(ID %in% RA_healthy$ID)
  val_rauh <- val %>% filter(COHORT == "RA") %>% filter(ID %in% RA_unhealthy$ID)
  val_all <- val %>% filter(COHORT %in% c("SEA", "STL", "RA", "KCL"))
  val_stl1 <- val %>% filter(COHORT %in% c("STL1"))

  r2 <- safe_r2(val_all$pred1, val_all$yy)
  r2_rauh <- safe_r2(val_rauh$pred1, val_rauh$yy)
  r2_rah <- safe_r2(val_rah$pred1, val_rah$yy)
  r2_ra <- safe_r2(val_ra$pred1, val_ra$yy)
  r2_stl <- safe_r2(val_stl$pred1, val_stl$yy)
  r2_sea <- safe_r2(val_sea$pred1, val_sea$yy)
  r2_kcl <- safe_r2(val_kcl$pred1, val_kcl$yy)
  r2_stl_sea <- safe_r2(val_stl_sea$pred1, val_stl_sea$yy)
  r2_stl1 <- safe_r2(val_stl1$pred1, val_stl1$yy)

  r2_ <- c(r2_rah, r2_rauh, r2_ra, r2_stl, r2_sea, r2_kcl, r2_stl_sea, r2, r2_stl1)
  names(r2_) <- c("r2_rah", "r2_rauh", "r2_ra", "r2_stl", "r2_sea", "r2_kcl", "r2_stl_sea", "r2", "r2_stl1")

  pheb_pred <- val %>% filter(COHORT == "UKBB") %>% select(ID, pred1)
  colnames(pheb_pred) <- c("eid", N)

  return(list(
    pheb_pred, val, r2_, coef_df,
    best_alpha = best_alpha, best_lambda = best_lambda,
    all_alpha_results = all_alpha_results,
    scale_params = scale_params,
    do_scale = do_scale
  ))
}

cb <- function(d) {
  prot_cols <- setdiff(names(d), c("ID", "COHORT", "SEX", "BL_AGE", "pheb", "phek", "rat"))
  mat <- t(as.matrix(d[, prot_cols]))
  batch <- d$COHORT
  combat_mat <- ComBat(
    dat = mat,
    batch = batch,
    par.prior = TRUE
  )

  d_corrected <- d
  d_corrected[, prot_cols] <- t(combat_mat)
  return(d_corrected)
}

PCA <- function(data) {
  dd <- data[, 2:(ncol(data) - 1)]

  pca <- prcomp(dd, center = TRUE, scale = TRUE)
  pca_ <- as.data.frame(pca$x[, 1:10])
  imp <- summary(pca)$importance[2, ]

  data <- cbind(data, pca_)
  return(list(data, imp))
}

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

plot_pca <- function(df, ukbb_n = 80, seed = 42,
                     palette = c(
                       SEA = "#E69F00",
                       STL = "#56B4E9",
                       RA = "#009E73",
                       UKBB = "grey50",
                       KCL = "#CC79A7"
                     )) {
  set.seed(seed)

  non_ukbb_df <- df %>% filter(COHORT != "UKBB")
  ukbb_df <- df %>% filter(COHORT == "UKBB")

  if (nrow(ukbb_df) > 0) {
    ukbb_df <- ukbb_df %>% slice_sample(n = min(ukbb_n, nrow(ukbb_df)))
    pca_df <- bind_rows(non_ukbb_df, ukbb_df)
  } else {
    pca_df <- non_ukbb_df
  }

  ggplot(pca_df, aes(PC1, PC2, colour = COHORT)) +
    geom_point(size = 2, alpha = 0.8) +
    scale_colour_manual(values = palette, name = NULL) +
    labs(x = "PC 1", y = "PC 2") +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "right",
      legend.key.width = unit(0.8, "cm")
    )
}

union <- function(seattle_sc, proteom_stl, RA, KCL) {
  sea <- seattle_sc$ID
  stl <- proteom_stl$ID
  ra <- RA$ID
  kcl <- KCL$ID
  col <- intersect(colnames(seattle_sc), colnames(proteom_stl))
  col <- intersect(col, colnames(RA))
  col <- intersect(col, colnames(KCL))
  RA_U <- RA %>% select(all_of(col))
  seattle_sc_U <- seattle_sc %>% select(all_of(col))
  proteom_stl_U <- proteom_stl %>% select(all_of(col))
  proteom_kcl_U <- KCL %>% select(all_of(col))

  UNION <- bind_rows(proteom_stl_U, seattle_sc_U)
  UNION <- bind_rows(UNION, RA_U)
  UNION <- bind_rows(UNION, proteom_kcl_U)
  UNION <- UNION %>% mutate(COHORT = case_when(
    ID %in% c(sea) ~ "SEA",
    ID %in% c(stl) ~ "STL",
    ID %in% c(ra) ~ "RA",
    ID %in% c(kcl) ~ "KCL"
  ))
  return(UNION)
}

protein_call_rate <- function(proteom_UKBB) {
  d <- nrow(proteom_UKBB)
  na_by_col <- summarise(proteom_UKBB, across(everything(), ~ sum(is.na(.)))) %>%
    pivot_longer(
      cols = everything(),
      names_to = "variable",
      values_to = "n_na"
    ) %>%
    mutate(protein_call_rate = 1 - n_na / d)
  return(na_by_col)
}

sample_call_rate <- function(proteom_UKBB) {
  dd <- ncol(proteom_UKBB) - 1
  na_by_row <- proteom_UKBB %>%
    mutate(na_ctn_row = rowSums(is.na(across(everything())))) %>%
    select(ID, na_ctn_row) %>%
    mutate(sample_call_rate = 1 - na_ctn_row / dd)
  return(na_by_row)
}

union1 <- function(seattle_sc, proteom_stl, proteom_UKBB, RA, KCL) {
  sea <- seattle_sc$ID
  stl <- proteom_stl$ID
  ukbb <- proteom_UKBB$ID
  ra <- RA$ID
  kcl <- KCL$ID

  col <- intersect(colnames(seattle_sc), colnames(proteom_stl))
  col <- intersect(col, colnames(proteom_UKBB))
  col <- intersect(col, colnames(RA))
  col <- intersect(col, colnames(KCL))
  RA_U <- RA %>% select(all_of(col))
  seattle_sc_U <- seattle_sc %>% select(all_of(col))
  proteom_stl_U <- proteom_stl %>% select(all_of(col))
  proteom_UKBB_U <- proteom_UKBB %>% select(all_of(col))
  proteom_KCL <- KCL %>% select(all_of(col))

  proteom_UKBB_U$ID <- as.character(proteom_UKBB_U$ID)
  UNION <- bind_rows(proteom_stl_U, seattle_sc_U)
  UNION <- bind_rows(UNION, RA_U)
  UNION <- bind_rows(UNION, proteom_UKBB_U)
  UNION <- bind_rows(UNION, proteom_KCL)
  UNION <- UNION %>% mutate(COHORT = case_when(
    ID %in% c(sea) ~ "SEA",
    ID %in% c(stl) ~ "STL",
    ID %in% c(ukbb) ~ "UKBB",
    ID %in% c(ra) ~ "RA",
    ID %in% c(kcl) ~ "KCL"
  ))
  return(UNION)
}

"%!in%" <- function(x, y) !("%in%"(x, y))

seattle_sc <- as.data.frame(fread(cli$seattle_pct, header = T, sep = ","))
seattle_sc <- seattle_sc %>%
  select(donor_id, minor_celltypes, pct) %>%
  pivot_wider(
    names_from = minor_celltypes,
    values_from = pct,
    values_fill = NA
  )
seattle_sc <- seattle_sc %>% select(donor_id, `Tem GZMB+`, `Tem GZMK+`)
colnames(seattle_sc) <- c("donor_id", "pheb", "phek")

seattle_sc$rat <- seattle_sc$phek / seattle_sc$pheb
seattle_sc$ratio <- seattle_sc$pheb / seattle_sc$phek

seattle_proteom <- read_excel(cli$seattle_olink, sheet = 1)
SEX_AGE <- seattle_proteom %>% select(Subject, `Baseline Age`, sex)
SEX_AGE <- SEX_AGE[!duplicated(SEX_AGE), ]
colnames(SEX_AGE) <- c("donor_id", "BL_AGE", "SEX")
SEX_AGE$SEX <- ifelse(SEX_AGE$SEX == "male", 1, 0)
seattle_sc <- as.data.frame(merge(seattle_sc, SEX_AGE, by = "donor_id"))
colnames(seattle_sc)[1] <- "Subject"

seattle_proteom_F1D0 <- seattle_proteom %>%
  filter(Visit == "Flu Year 1 Day 0") %>%
  select(Subject, NPX_bridged, Assay) %>%
  pivot_wider(
    names_from = Assay,
    values_from = NPX_bridged,
    values_fn = mean,
    values_fill = NA
  )

seattle_proteom_F1D7 <- seattle_proteom %>%
  filter(Visit == "Flu Year 1 Day 7") %>%
  select(Subject, NPX_bridged, Assay) %>%
  pivot_wider(
    names_from = Assay,
    values_from = NPX_bridged,
    values_fn = mean,
    values_fill = NA
  )

seattle_proteom_F1D90 <- seattle_proteom %>%
  filter(Visit == "Flu Year 1 Day 90") %>%
  select(Subject, NPX_bridged, Assay) %>%
  pivot_wider(
    names_from = Assay,
    values_from = NPX_bridged,
    values_fn = mean,
    values_fill = NA
  )

seattle_proteom_F2D0 <- seattle_proteom %>%
  filter(Visit == "Flu Year 2 Day 0") %>%
  select(Subject, NPX_bridged, Assay) %>%
  pivot_wider(
    names_from = Assay,
    values_from = NPX_bridged,
    values_fn = mean,
    values_fill = NA
  )

seattle_proteom_F2D7 <- seattle_proteom %>%
  filter(Visit == "Flu Year 2 Day 7") %>%
  select(Subject, NPX_bridged, Assay) %>%
  pivot_wider(
    names_from = Assay,
    values_from = NPX_bridged,
    values_fn = mean,
    values_fill = NA
  )

seattle_proteom_F2D90 <- seattle_proteom %>%
  filter(Visit == "Flu Year 2 Day 90") %>%
  select(Subject, NPX_bridged, Assay) %>%
  pivot_wider(
    names_from = Assay,
    values_from = NPX_bridged,
    values_fn = mean,
    values_fill = NA
  )

seattle_proteom_F1SA <- seattle_proteom %>%
  filter(Visit == "Flu Year 1 Stand-Alone") %>%
  select(Subject, NPX_bridged, Assay) %>%
  pivot_wider(
    names_from = Assay,
    values_from = NPX_bridged,
    values_fn = mean,
    values_fill = NA
  )

seattle_proteom_F2SA <- seattle_proteom %>%
  filter(Visit == "Flu Year 2 Stand-Alone") %>%
  select(Subject, NPX_bridged, Assay) %>%
  pivot_wider(
    names_from = Assay,
    values_from = NPX_bridged,
    values_fn = mean,
    values_fill = NA
  )

seattle_proteom_F3SA <- seattle_proteom %>%
  filter(Visit == "Flu Year 3 Stand-Alone") %>%
  select(Subject, NPX_bridged, Assay) %>%
  pivot_wider(
    names_from = Assay,
    values_from = NPX_bridged,
    values_fn = mean,
    values_fill = NA
  )

seattle_proteom_OTH <- seattle_proteom %>%
  filter(Visit == "Other - Non-Flu") %>%
  select(Subject, NPX_bridged, Assay) %>%
  pivot_wider(
    names_from = Assay,
    values_from = NPX_bridged,
    values_fn = mean,
    values_fill = NA
  )

dop1 <- seattle_proteom %>%
  filter(Visit == "Other - Non-Flu" & Subject %in% c("BR1029", "BR1035", "BR2004")) %>%
  group_by(Subject, Assay) %>%
  dplyr::summarise(NPX_bridged = mean(NPX_bridged)) %>%
  select(Subject, NPX_bridged, Assay) %>%
  pivot_wider(
    names_from = Assay,
    values_from = NPX_bridged,
    values_fn = mean,
    values_fill = NA
  )
dop2 <- seattle_proteom %>%
  filter(Visit == "Flu Year 1 Stand-Alone" & Subject %in% c("BR2027")) %>%
  group_by(Subject, Assay) %>%
  dplyr::summarise(NPX_bridged = mean(NPX_bridged)) %>%
  select(Subject, NPX_bridged, Assay) %>%
  pivot_wider(
    names_from = Assay,
    values_from = NPX_bridged,
    values_fn = mean,
    values_fill = NA
  )
dop <- rbind(dop1, dop2)

dop_sc <- as.data.frame(merge(seattle_sc, dop, by = "Subject"))

seattle_sc <- as.data.frame(merge(seattle_sc, seattle_proteom_F1D0, by = "Subject"))
seattle_sc <- rbind(seattle_sc, dop_sc)
colnames(seattle_sc)[1] <- "ID"
seattle_sc$ID <- as.character(seattle_sc$ID)

na_by_col <- summarise(seattle_sc, across(everything(), ~ sum(is.na(.)))) %>%
  pivot_longer(
    cols = everything(),
    names_to = "variable",
    values_to = "n_na"
  )
p_ex <- na_by_col[na_by_col$n_na > 1, ]$variable
na_by_row <- seattle_sc %>%
  mutate(na_ctn_row = rowSums(is.na(across(everything())))) %>%
  select(ID, na_ctn_row)
seattle_sc <- seattle_sc %>% select(-all_of(p_ex)) %>% filter(ID != "BR2036")
colnames(seattle_sc) <- gsub("-", "_", colnames(seattle_sc))

sc <- as.data.frame(fread(cli$stl_annotations, header = T, sep = ","))
colnames(sc) <- c("ID", "pheb", "phek", "temra", "BL_AGE", "SEX", "rat", "ratio")
sc$ID <- as.character(sc$ID)
sc$SEX <- ifelse(sc$SEX == "Male", 1, 0)
sc$phebtemra <- sc$pheb + sc$temra
sc$rat2 <- sc$phek / sc$phebtemra
sc$ratio2 <- sc$phebtemra / sc$phek

for (name in c("pheb", "phek", "temra", "rat", "rat2", "phebtemra", "ratio", "ratio2")) {
  ranks <- rank(sc[, name], ties.method = "average")
  sc[, paste(name, "_IRNT", sep = "")] <- qnorm((ranks - 0.5) / nrow(sc))
}
proteom_stl <- read.csv(cli$stl_olink)
proteom_stl <- proteom_stl %>% filter(Assay %!in% c(
  "Extension control 1", "Incubation control 1", "Amplification control 1",
  "Extension control 2", "Incubation control 2", "Amplification control 2",
  "Extension control 3", "Incubation control 3", "Amplification control 3",
  "Extension control 4", "Incubation control 4", "Amplification control 4"
))

proteom_stl1 <- proteom_stl %>%
  select(SampleID, SampleType, Assay, NPX) %>%
  pivot_wider(
    names_from = Assay,
    values_from = NPX,
    values_fn = mean,
    values_fill = NA
  )
proteom_stl1 <- proteom_stl1 %>% select_if(~ !any(is.na(.)))
colnames(proteom_stl1)[1] <- "ID"
proteom_stl1 <- as.data.frame(merge(sc, proteom_stl1, by = "ID"))

dd <- proteom_stl1[, 21:(ncol(proteom_stl1))]
dd[is.infinite(as.matrix(dd))] <- NA
dd <- dd[, sapply(dd, function(x) var(x, na.rm = TRUE) != 0)]

pca <- prcomp(dd, center = TRUE, scale = TRUE)
pca_ <- as.data.frame(pca$x[, 1:10])

dd <- cbind(proteom_stl1, pca_)

SUS <- c("E12", "FA05", "FE01", "FE06", "FE07", "FH04", "G04", "H18", "M14", "P15")
proteom_stl <- proteom_stl %>% filter(SampleID %!in% SUS)
proteom_stl <- proteom_stl %>%
  select(SampleID, SampleType, Assay, NPX) %>%
  pivot_wider(
    names_from = Assay,
    values_from = NPX,
    values_fn = mean,
    values_fill = NA
  )
setdiff(proteom_stl$SampleID, sc$ID)

colnames(proteom_stl)[1] <- "ID"
proteom_stl$ID <- as.character(proteom_stl$ID)
proteom_stl <- as.data.frame(merge(sc, proteom_stl, by = "ID"))
proteom_stl <- proteom_stl %>% select_if(~ !any(is.na(.)))
colnames(proteom_stl) <- gsub("-", "_", colnames(proteom_stl))

data <- as.data.frame(fread(cli$kcl_annotations, header = T, sep = ","))
data <- data %>% filter(Dataset == "KCL")

data <- data %>%
  select(donor_id, Annotation, percent, Sex, Ethnicity, Age) %>%
  pivot_wider(
    names_from = Annotation,
    values_from = percent,
    values_fn = mean,
    values_fill = NA
  )
data <- data %>%
  select(donor_id, `CD8 Tem GZMB+`, `CD8 Tem GZMK+`, Age, Sex) %>%
  mutate(rat = `CD8 Tem GZMK+` / `CD8 Tem GZMB+`) %>%
  select(
    donor_id,
    `CD8 Tem GZMB+`, `CD8 Tem GZMK+`, rat,
    Age, Sex
  ) %>%
  mutate(Sex = case_when(
    Sex == "Male" ~ 1,
    Sex == "Female" ~ 0
  ))
colnames(data) <- c("ID", "pheb", "phek", "rat", "BL_AGE", "SEX")
data$ID <- as.character(data$ID)

data_prot <- as.data.frame(fread(cli$kcl_olink, header = T, sep = ","))

data_prot <- data_prot %>% filter(Assay %!in% c(
  "Extension control 1", "Incubation control 1", "Amplification control 1",
  "Extension control 2", "Incubation control 2", "Amplification control 2",
  "Extension control 3", "Incubation control 3", "Amplification control 3",
  "Extension control 4", "Incubation control 4", "Amplification control 4",
  "Extension control 5", "Incubation control 5", "Amplification control 5",
  "Extension control 6", "Incubation control 6", "Amplification control 6",
  "Extension control 7", "Incubation control 7", "Amplification control 7",
  "Extension control 8", "Incubation control 8", "Amplification control 8"
))


data_prot <- data_prot %>%
  select(SampleID, SampleType, Assay, NPX) %>%
  pivot_wider(
    names_from = Assay,
    values_from = NPX,
    values_fn = mean,
    values_fill = NA
  )

names(data_prot)[names(data_prot) == "MICA_MICB"] <- "MICB_MICA"
colnames(data_prot) <- gsub("-", "_", colnames(data_prot))


data_prot <- data_prot %>% filter(SampleType == "SAMPLE")
data_prot <- data_prot %>% select_if(~ !any(is.na(.)))
colnames(data_prot)[1] <- "ID"
data_prot$ID <- as.character(data_prot$ID)

data_prot$SampleType <- NULL
data_prot$COHORT <- "KCL"
data_prot <- as.data.frame(merge(data, data_prot, by = "ID"))
dd <- data_prot[, 8:(ncol(data_prot)) - 1]
dd[is.infinite(as.matrix(dd))] <- NA
dd <- dd[, sapply(dd, function(x) var(x, na.rm = TRUE) != 0)]

pca <- prcomp(dd, center = TRUE, scale = TRUE)
pca_ <- as.data.frame(pca$x[, 1:10])

dd <- cbind(data_prot, pca_)

excl <- c(
  "A023", "A004", "A027", "A019", "A133", "A005", "A142", "A137", "A136", "A130", "A013", "A003", "A020", "A007",
  "A002", "A026", "A008", "A132", "A022", "A009", "A135", "A018", "A011", "A138", "A006"
)

dd %>%
  mutate(excl = case_when(
    ID %in% excl ~ "yes",
    TRUE ~ "no"
  )) %>%
  ggplot(aes(PC1, PC2, col = excl)) +
  geom_point() +
  scale_colour_manual(values = c(yes = "steelblue", no = "grey")) +
  theme_classic()

data_prot <- data_prot %>% filter(ID %!in% excl)
dd <- data_prot[, 3:(ncol(data_prot)) - 1]
dd[is.infinite(as.matrix(dd))] <- NA
dd <- dd[, sapply(dd, function(x) var(x, na.rm = TRUE) != 0)]

pca <- prcomp(dd, center = TRUE, scale = TRUE)
pca_ <- as.data.frame(pca$x[, 1:10])

dd <- cbind(data_prot, pca_)
dd %>%
  mutate(excl = case_when(
    ID %in% excl ~ "yes",
    TRUE ~ "no"
  )) %>%
  ggplot(aes(PC1, PC2, col = excl)) +
  geom_point() +
  scale_colour_manual(values = c(yes = "steelblue", no = "grey")) +
  theme_classic()

excl1 <- (dd %>% filter(PC2 < -50))$ID

data_prot <- data_prot %>% filter(ID %!in% excl1)
dd <- data_prot[, 3:(ncol(data_prot)) - 1]
dd[is.infinite(as.matrix(dd))] <- NA
dd <- dd[, sapply(dd, function(x) var(x, na.rm = TRUE) != 0)]

pca <- prcomp(dd, center = TRUE, scale = TRUE)
pca_ <- as.data.frame(pca$x[, 1:10])

dd <- cbind(data_prot, pca_)

dd %>%
  mutate(excl = case_when(
    ID %in% excl ~ "yes",
    TRUE ~ "no"
  )) %>%
  ggplot(aes(PC1, PC2, col = excl)) +
  geom_point() +
  scale_colour_manual(values = c(yes = "steelblue", no = "grey")) +
  theme_classic()

RA_sc <- as.data.frame(fread(cli$ra_annotations, header = T, sep = ","))
RA_proteom <- as.data.frame(fread(cli$ra_olink, header = T, sep = ","))
RA_proteom <- RA_proteom %>%
  select(sample.sampleKitGuid, subject.biologicalSex, sample.subjectAgeAtDraw, olink.assay, olink.NPX_norm) %>%
  pivot_wider(
    names_from = olink.assay,
    values_from = olink.NPX_norm,
    values_fn = mean,
    values_fill = NA
  )
RA <- as.data.frame(merge(RA_sc, RA_proteom, by = "sample.sampleKitGuid"))
colnames(RA)[1] <- "ID"
RA$ID <- as.character(RA$ID)
colnames(RA)[2] <- "phek"
colnames(RA)[3] <- "pheb"
colnames(RA)[4] <- "rat"
colnames(RA)[5] <- "disease"
colnames(RA)[9] <- "SEX"
colnames(RA)[10] <- "BL_AGE"
RA$SEX <- ifelse(RA$SEX == "Male", 1, 0)

RA <- RA[, colSums(!is.na(RA)) > 0]
colnames(RA) <- gsub("-", "_", colnames(RA))

RA_healthy <- RA %>% filter(disease == "Control (HC1)")
RA_healthy <- RA_healthy[, colSums(!is.na(RA_healthy)) > 0]
RA_unhealthy <- RA %>% filter(disease != "Control (HC1)")
RA_unhealthy <- RA_unhealthy[, colSums(!is.na(RA_unhealthy)) > 0]
RA_healthy$ID <- as.character(RA_healthy$ID)
RA_unhealthy$ID <- as.character(RA_unhealthy$ID)

RA <- RA %>% select_if(~ !any(is.na(.)))
RA_healthy <- RA_healthy %>% select_if(~ !any(is.na(.)))

UNION_STLSEARAKCL <- union(seattle_sc, proteom_stl, RA_healthy, data_prot)
PROTEINS <- colnames(UNION_STLSEARAKCL[, 7:ncol(UNION_STLSEARAKCL)])


if (has_ukbb) {
  data4 <- as.data.frame(fread(cli$ukbb_metadata, header = T, select = c("eid", "SEX", "BL_AGE")))
  colnames(data4)[1] <- "ID"
  data4$ID <- as.character(data4$ID)
  EXCLUDED <- read.table(cli$ukbb_excluded, sep = "\t", header = F)
  colnames(EXCLUDED)[1] <- "eid"
  EXCLUDED$eid <- as.character(EXCLUDED$eid)
  data4 <- data4 %>% filter(ID %!in% EXCLUDED$eid)

  if (ukbb_proteome_mode == "single_file") {
    proteom <- as.data.frame(fread(cli$ukbb_olink_file, header = TRUE))
  } else if (ukbb_proteome_mode == "six_instances") {
    proteom <- fread(cli$ukbb_olink_instance_4, header = TRUE)
    proteom1 <- fread(cli$ukbb_olink_instance_3, header = TRUE)
    proteom2 <- fread(cli$ukbb_olink_instance_1, header = TRUE)
    proteom3 <- fread(cli$ukbb_olink_instance_2, header = TRUE)
    proteom4 <- fread(cli$ukbb_olink_instance_5, header = TRUE)
    proteom5 <- fread(cli$ukbb_olink_instance_6, header = TRUE)
    proteom <- as.data.frame(merge(proteom, proteom1, by = "eid"))
    proteom <- as.data.frame(merge(proteom, proteom2, by = "eid"))
    proteom <- as.data.frame(merge(proteom, proteom3, by = "eid"))
    proteom <- as.data.frame(merge(proteom, proteom4, by = "eid"))
    proteom <- as.data.frame(merge(proteom, proteom5, by = "eid"))
  } else {
    stop("Internal error: unknown UKBB proteome mode.", call. = FALSE)
  }

  colnames(proteom) <- toupper(colnames(proteom))
  names(proteom)[names(proteom) == "EID"] <- "ID"
  if (!"ID" %in% names(proteom)) {
    colnames(proteom)[1] <- "ID"
  }
  names(proteom)[names(proteom) == "C19ORF12"] <- "C19orf12"
  proteom$ID <- as.character(proteom$ID)
  proteom$COHORT <- "UKBB"

  proteom <- proteom %>% filter(ID %!in% EXCLUDED$eid)

  PROTEINS_UKBB <- intersect(PROTEINS, colnames(proteom))

  proteom_UKBB <- proteom %>% select("ID", all_of(PROTEINS_UKBB))
  proteom_UKBB$COHORT <- NULL

  protein_qc <- protein_call_rate(proteom_UKBB)
  sample_qc <- sample_call_rate(proteom_UKBB)

  d <- sample_qc %>% filter(ID %in% EXCLUDED$eid)

  group1 <- sample_qc %>% filter(sample_call_rate > 0.9)
  prot1 <- protein_qc %>% filter(protein_call_rate > 0.9)


  med_cr <- median(sample_qc$sample_call_rate, na.rm = TRUE)

  ggplot(sample_qc, aes(sample_call_rate)) +
    geom_histogram(
      aes(y = after_stat(density)),
      bins = 40,
      fill = "#4DBBD5FF",
      colour = "white",
      linewidth = 0.3
    ) +
    geom_density(colour = "#1A2C42", linewidth = 0.8, adjust = 1.3) +
    geom_vline(
      xintercept = med_cr,
      linetype = "dashed",
      colour = "#1A2C42"
    ) +
    annotate(
      "text",
      x = med_cr,
      y = Inf, vjust = 1.4, hjust = -0.05,
      label = paste0("median = ", round(med_cr, 3)),
      colour = "#1A2C42",
      size = 3.5
    ) +
    scale_x_continuous(
      labels = percent_format(accuracy = 1),
      limits = c(0.90, 1.00), expand = c(0, 0)
    ) +
    labs(
      x = "Sample non-missingness rate",
      y = "Density",
      title = "Distribution of sample call rate"
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title.position = "plot",
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )

  prot_df <- protein_qc %>% filter(variable != "ID")
  med_pcr <- median(prot_df$protein_call_rate, na.rm = TRUE)

  ggplot(prot_df, aes(protein_call_rate)) +
    geom_histogram(
      aes(y = after_stat(density)),
      bins = 40,
      fill = "#4DBBD5FF",
      colour = "white",
      linewidth = 0.3
    ) +
    geom_density(colour = "#1A2C42", linewidth = 0.8, adjust = 1.3) +
    geom_vline(
      xintercept = med_pcr,
      linetype = "dashed",
      colour = "#1A2C42"
    ) +
    annotate(
      "text",
      x = med_pcr,
      y = Inf,
      vjust = 1.4,
      hjust = -0.05,
      label = paste0("median = ", round(med_pcr, 3)),
      colour = "#1A2C42",
      size = 3.5
    ) +
    scale_x_continuous(
      labels = percent_format(accuracy = 1),
      limits = c(0.90, 1.00),
      expand = c(0, 0)
    ) +
    labs(
      x = "Protein non-missingness rate",
      y = "Density",
      title = "Distribution of protein call rate (excluding ID)"
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title.position = "plot",
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )

  proteom_UKBB1 <- proteom_UKBB %>% filter(ID %in% group1$ID)
  proteom_UKBB1 <- proteom_UKBB1 %>% select(all_of(prot1$variable))

  d <- proteom_UKBB1 %>% filter(ID %in% EXCLUDED$eid)

  protein_qc1 <- protein_call_rate(proteom_UKBB1)
  sample_qc1 <- sample_call_rate(proteom_UKBB1)

  sample_qc1 %>%
    ggplot(aes(sample_call_rate)) +
    geom_histogram() +
    theme_classic()

  protein_qc1 %>%
    filter(variable != "ID") %>%
    ggplot(aes(protein_call_rate)) +
    geom_histogram() +
    theme_classic()

  proteom_UKBB1_imp <- proteom_UKBB1 %>%
    mutate(across(
      where(is.numeric),
      ~ replace_na(., mean(., na.rm = TRUE))
    ))
} else {
  data4 <- NULL
  proteom_UKBB1_imp <- NULL
}

if (has_ukbb) {
  UNION1 <- union1(seattle_sc, proteom_stl, proteom_UKBB1_imp, RA_healthy, data_prot)
} else {
  UNION1 <- UNION_STLSEARAKCL %>% select("ID", all_of(PROTEINS))
  if (!"COHORT" %in% names(UNION1)) {
    UNION1$COHORT <- UNION_STLSEARAKCL$COHORT
  }
}


set.seed(42)
UNION_cor1 <- cb(UNION1)

pca_cor1_80 <- PCA(UNION_cor1)
pca_plot <- plot_pca(pca_cor1_80[[1]])
ensure_parent_dir(cli$out_pca_plot)
ggsave(filename = cli$out_pca_plot, plot = pca_plot, width = 7, height = 5)
pca_plot
message("PCA variance explained, first 10 PCs:")


message("Building model input table.")
SINGLE_CELL <- bind_rows(sc[, c("ID", "pheb", "phek", "BL_AGE", "SEX")], seattle_sc[, c("ID", "pheb", "phek", "BL_AGE", "SEX")])
SINGLE_CELL <- bind_rows(SINGLE_CELL, RA[, c("ID", "pheb", "phek", "BL_AGE", "SEX")])
SINGLE_CELL <- bind_rows(SINGLE_CELL, data_prot[, c("ID", "pheb", "phek", "BL_AGE", "SEX")])
message("Merging single-cell annotations with ComBat-corrected proteome.")
UNION_cor1_f <- as.data.frame(merge(SINGLE_CELL, UNION_cor1, by = "ID", all.y = T))


if (has_ukbb) {
  message("Adding UKBB metadata.")
  data4$SEX <- ifelse(data4$SEX == 2, 0, 1)
  UNION_cor1_f <- as.data.frame(merge(UNION_cor1_f, data4, by = "ID", all.x = T))
}

message("Coalescing phenotype and covariate columns.")
UNION_cor1_f <- coalesce_xy_columns(UNION_cor1_f, c("pheb", "phek", "BL_AGE", "SEX", "rat"))


write_tsv(UNION_cor1_f, cli$out_union_cor1_f)

message("Fitting lasso model.")
set.seed(42)
pred_lasso_all <- get_lasso_model(UNION_cor1_f, c("SEA", "STL", "KCL", "RAH"), "pheb_combat_1_lasso_stlseakclrah_v3", "pheb", alpha_grid = 1)
if (has_ukbb) {
  write_tsv(pred_lasso_all[[1]], cli$out_ukbb_predictions)
} else {
  message("Skipping UKBB prediction output because UKBB inputs were not provided.")
}
r2_metrics <- as.data.frame(t(pred_lasso_all[[3]]))
write_tsv(r2_metrics, cli$out_r2_metrics)
coefficient_table <- pred_lasso_all[[4]]
coefficient_table$feature <- rownames(coefficient_table)
coefficient_table <- coefficient_table[, c("feature", setdiff(names(coefficient_table), "feature"))]
write_tsv(coefficient_table, cli$out_coefficients)
ensure_parent_dir(cli$out_lasso_model)
saveRDS(pred_lasso_all, cli$out_lasso_model)
r2_plot <- get_R2_plots(pred_lasso_all, "r2", c("KCL", "STL", "SEA", "RAH"), "GZMB")
ensure_parent_dir(cli$out_r2_plot)
ggsave(filename = cli$out_r2_plot, plot = r2_plot, width = 10, height = 7)

set.seed(42)
dd <- cross_valid(UNION_cor1_f %>% filter(COHORT %in% c("SEA", "STL", "RAH", "KCL")), "SEA", "STL", "KCL", "RAH", "E", "E", "min", 100, 0.3, alpha_grid = 1, "pheb")
write_tsv(dd, cli$out_cross_validation)
cross_validation_plot <- get_cross_val_plot(dd)
ensure_parent_dir(cli$out_crossval_plot)
ggsave(filename = cli$out_crossval_plot, plot = cross_validation_plot, width = 10, height = 7)
cross_validation_plot
