library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(ggalluvial)
library(ggpubr)
library(ggrepel)
library(ggbreak)
library(ggrastr)
library(rstatix)
library(broom)
library(variancePartition)
library(reshape2)
library(factoextra)
library(forcats)
library(tidyverse)
library(mascarade)
# -----------------------------
# Configuration
# -----------------------------

dataset_levels = c( "OneK1K", "AIDA", "CIMA", "KCL", "ABF300", "UCSF", "SoundLife", "SPAC")
dataset_colors <- c("OneK1K" = "#6BA6A8", "AIDA" = "#E69A45", "SoundLife" = "#C4B454", "SPAC" = "#A78CC1",
                    "UCSF" = "#C15B58", "ABF300" = "#8FB98B","CIMA" = "#A68A74", "KCL" = "#BF4294")
density_colors <- c("OneK1K" = "#487173", "AIDA" = "#c97e28", "SoundLife" = "#968833", "SPAC" = "#5a3082",
                    "UCSF" = "#9c322f", "ABF300" = "#31732a", "KCL" = "#8c2a6a", "CIMA" = "#75573f")
celltype_colors_major <- rev(c( "#00A087", "#4DBBD5", "#4DBBD5", "#7CAE00",  "#c7345d", "#8491B4",
                                "#F4CAE4", "#E64B35","#3C5488", "#a39b83", "#EEC900", "#B15928"))

data <- read.csv('ALL_DATASETS_ANNOTATION.csv')

# -----------------------------
# Figure 1, S1
# -----------------------------
plot_donor_sex_by_dataset <- function(data) {
  donor_n_text <- data %>%
    distinct(Dataset, donor_id, sex) %>%
    count(Dataset, sex) %>%
    mutate(label = paste0("n = ", n))
  
  donor_n_text$Dataset <- factor(
    donor_n_text$Dataset,
    levels = dataset_levels
  )
  
  ggplot(donor_n_text, aes(x = Dataset, y=n, fill=Dataset)) +  
    geom_col_pattern(aes(pattern = sex), color='black') +
    theme_minimal() + ylab('') +  
    xlab('') + theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black")) +
    scale_fill_manual(values=dataset_colors)
}

plot_total_cells_by_dataset <- function(data) {
  n_cells <- as.data.frame(table(data$Dataset))
  colnames(n_cells) <- c("Dataset", "N_cells")
  
  n_cells$Dataset <- factor(
    n_cells$Dataset,
    levels = dataset_levels
  )
  
  ggplot(n_cells, aes(x = Dataset, y = N_cells, fill = Dataset)) +
    geom_col(color = "black") +
    scale_y_break(c(2000000, 5000000), scales = 1) +
    scale_fill_manual(values = dataset_colors) +
    theme_minimal(base_size = 11) +
    labs(x = NULL, y = NULL) +
    theme(
      axis.text.x = element_text(colour = "black", angle = 45, hjust = 1),
      axis.text.y = element_text(colour = "black"),
      legend.key.size = unit(0.7, "cm")
    )
}

plot_age_histograms_by_dataset <- function(data) {
  df <- data %>%
    distinct(Dataset, donor_id, age)
  
  x_limits <- c(18, 99)
  
  make_hist <- function(dataset_name) {
    ggplot(df %>% filter(Dataset == dataset_name), aes(x = age)) +
      geom_histogram(
        binwidth = 3,
        fill = dataset_colors[[dataset_name]],
        color = "black",
        boundary = 0,
        closed = "left"
      ) +
      geom_density(aes(y = after_stat(count * 3)),
                   color = density_colors[[dataset_name]],
                   linewidth = 1.1) +
      scale_x_continuous(limits = x_limits, breaks = seq(18, 98, 5)) +
      theme_minimal(base_size = 11) +
      labs(x = NULL, y = NULL, title = dataset_name)
  }
  
  plots <- lapply(names(dataset_colors), make_hist)
  patchwork::wrap_plots(plots, ncol = 1)
}

plot_ethnicity_distribution <- function(data) {
  df <- data %>%
    distinct(donor_id, ethnicity) %>%
    count(ethnicity, name = "Count") %>%
    rename(Ethnicity = ethnicity)
  
  ggplot(df, aes(x = Ethnicity, y = Count, fill = Ethnicity)) +
    geom_col(color = "black") +
    scale_y_break(c(50, 600), scales = 1) +
    theme_minimal(base_size = 11) +
    labs(title = "Ancestry Distribution", x = NULL, y = "Count") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
      legend.position = "none"
    )
}

plot_bad_quality <- function(data) {
  df_clean <- data %>%
    mutate(quality = ifelse(major_celltypes %in% c("Bad quality"), "Bad quality cells", "Rest"))
  
  df_percent <- df_clean %>%
    group_by(Dataset, quality) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(Dataset) %>%
    mutate(percentage = n / sum(n) * 100)
  
  ggplot(df_percent, aes(x = Dataset, y = percentage, fill = quality)) +
    geom_col(position = "stack", color = "black") +
    ylab("Percentage") +
    scale_fill_manual(values = c("Rest" = "#00bf7dff", "Bad quality cells" = "#ffe680ff")) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(colour = "black", angle = 45, hjust = 1))
}

plot_annotation_alluvial <- function(data, desired_order, celltype_colors) {
  df <- data %>% filter(major_celltypes != "Bad quality")
  df <- df %>% filter(Dataset %in% c( "OneK1K", "AIDA", "CIMA", "UCSF", "SoundLife", "SPAC"))
  
  df_count <- df %>%
    count(author_celltype, minor_celltypes)
  
  mapping <- df_count %>%
    group_by(author_celltype) %>%
    slice_max(order_by = n, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  ordered_author <- mapping %>%
    arrange(minor_celltypes) %>%
    pull(author_celltype)
  
  df$author_celltype <- factor(df$author_celltype, levels = ordered_author)
  
  df_alluv <- df %>%
    count(Dataset, minor_celltypes, author_celltype) %>%
    filter(n > 10)
  
  ggplot(
    df_alluv,
    aes(axis1 = minor_celltypes, axis2 = author_celltype, y = n)
  ) +
    geom_alluvium(aes(fill = minor_celltypes), width = 1 / 12) +
    geom_stratum(width = 1 / 12, fill = "gray80", color = "black") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 2.5) +
    scale_x_discrete(
      limits = c("New", "Author cell type"),
      expand = c(0.05, 0.05)
    ) +
    facet_wrap(~Dataset, scales = "free") +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      strip.text = element_text(size = 10),
      axis.text.x = element_text(size = 10, color = "black")
    ) +
    labs(y = "Cell count", x = NULL)
}


# -----------------------------
# Figure 2, S2, S3
# -----------------------------

# -----------------------------
# Calculate major population percentages
# % from total PBMC per donor
# -----------------------------

calculate_major_percentages <- function(data) {
  
  df <- data %>%
    mutate(Annotation = major_celltypes) %>%
    filter(Annotation != "Bad quality")
  
  full_grid <- df %>%
    distinct(Dataset, donor_id) %>%
    crossing(Annotation = unique(df$Annotation))
  
  cell_counts <- df %>%
    group_by(Dataset, donor_id, Annotation) %>%
    summarise(n_cells = n(), .groups = "drop")
  
  total_cells <- df %>%
    group_by(Dataset, donor_id) %>%
    summarise(total_cells = n(), .groups = "drop")
  
  major_percentages <- full_grid %>%
    left_join(cell_counts, by = c("Dataset", "donor_id", "Annotation")) %>%
    left_join(total_cells, by = c("Dataset", "donor_id")) %>%
    mutate(
      n_cells = replace_na(n_cells, 0),
      percent = 100 * n_cells / total_cells
    )
  
  return(major_percentages)
}

# -----------------------------
# Calculate minor population percentages
# % from corresponding major population per donor
# -----------------------------
`%ni%` = Negate(`%in%`)

calculate_minor_percentages <- function(data) {
  
  df <- data %>% 
    mutate(
      Major = major_celltypes,
      Annotation = minor_celltypes
    ) %>%
    filter(
      Major %ni% c("Bad quality", "ILCs", "MAIT", "DN T", "Progenitor"),
      Annotation != "Bad quality"
    )
  
  full_grid <- df %>%
    distinct(Dataset, donor_id, Major) %>%
    left_join(
      df %>% distinct(Major, Annotation),
      by = "Major"
    )
  
  cell_counts <- df %>%
    group_by(Dataset, donor_id, Major, Annotation) %>%
    summarise(n_cells = n(), .groups = "drop")
  
  total_major_cells <- df %>%
    group_by(Dataset, donor_id, Major) %>%
    summarise(total_major_cells = n(), .groups = "drop")
  
  minor_percentages <- full_grid %>%
    left_join(cell_counts, by = c("Dataset", "donor_id", "Major", "Annotation")) %>%
    left_join(total_major_cells, by = c("Dataset", "donor_id", "Major")) %>%
    mutate(
      n_cells = replace_na(n_cells, 0),
      percent = 100 * n_cells / total_major_cells
    )
  
  return(minor_percentages)
}

major_percentages <- calculate_major_percentages(data)
minor_percentages <- calculate_minor_percentages(data)

### Heatmap ###

prepare_age_association_df <- function(data, major_percentages) {
  
  df_summary <- major_percentages
  
  make_ratio <- function(df, labels, annotation_name, denominator_cols) {
    
    df_wide <- df %>%
      filter(Annotation %in% names(labels)) %>%
      mutate(Type = recode(Annotation, !!!labels)) %>%
      group_by(Dataset, donor_id, Type) %>%
      summarise(n = sum(n_cells), .groups = "drop") %>%
      pivot_wider(names_from = Type, values_from = n, values_fill = 0)
    
    df_wide %>%
      mutate(
        denominator = rowSums(across(all_of(denominator_cols))),
        percent = ifelse(denominator == 0, NA, CD4 / denominator),
        Annotation = annotation_name
      ) %>%
      select(Dataset, donor_id, Annotation, percent)
  }
  
  ratio_rows <- make_ratio(
    df_summary,
    labels = c("CD4 T" = "CD4", "CD8 T" = "CD8", "MAIT" = "MAIT"),
    annotation_name = "CD4/CD8 T cells",
    denominator_cols = c("CD8", "MAIT")
  )
  
  ratio_rows_2 <- make_ratio(
    df_summary,
    labels = c("CD4 T" = "CD4", "CD8 T" = "CD8"),
    annotation_name = "CD4/TRAV1-2- CD8 T cells",
    denominator_cols = c("CD8")
  )
  
  meta <- data %>%
    select(donor_id, sex, age, ethnicity) %>%
    distinct()
  
  bind_rows(df_summary, ratio_rows, ratio_rows_2) %>%
    left_join(meta, by = "donor_id") %>%
    mutate(
      Dataset = factor(Dataset),
      Sex = factor(sex),
      Ethnicity = factor(ethnicity),
      Age = age
    ) %>%
    as.data.frame()
}

calculate_age_associations <- function(df) {
  model_results <- df %>%
    group_by(Dataset, Annotation) %>%
    do({
      subdf <- .
      
      r <- if (
        length(unique(subdf$Age)) > 1 &&
        length(unique(subdf$percent)) > 1
      ) {
        cor(subdf$Age, subdf$percent, method = "pearson", use = "complete.obs")
      } else {
        NA
      }
      
      if (n_distinct(subdf$Ethnicity) > 1) {
        model <- lm(percent ~ Age + Sex + Ethnicity, data = subdf)
      } else {
        model <- lm(percent ~ Age + Sex, data = subdf)
      }
      
      out <- tidy(model)
      out$r <- r
      out
    }) %>%
    filter(term == "Age") %>%
    ungroup() %>%
    group_by(Dataset) %>%
    mutate(p.adj = p.adjust(p.value, method = "BH")) %>%
    ungroup() %>%
    mutate(sig = case_when(
      is.na(p.adj) ~ " ",
      p.adj < 0.001 & estimate > 0 ~ "+++",
      p.adj < 0.001 & estimate < 0 ~ "---",
      p.adj < 0.01  & estimate > 0 ~ "++",
      p.adj < 0.01  & estimate < 0 ~ "--",
      p.adj < 0.05  & estimate > 0 ~ "+",
      p.adj < 0.05  & estimate < 0 ~ "-",
      TRUE ~ "ns"
    )) %>%
    select(Dataset, Annotation, estimate, std.error, statistic, p.value, p.adj, r, sig)
  
  model_results$Annotation[model_results$Annotation == "CD8 T"] <- "TRAV1-2- CD8 T cells"
  
  
  model_results$Dataset <- factor(
    model_results$Dataset,
    levels = dataset_levels
  )
  
  model_results$Annotation <- ifelse(
    !grepl("cells", model_results$Annotation) &
      model_results$Annotation != "ILCs",
    paste(model_results$Annotation, "cells"),
    model_results$Annotation
  )
  
  model_results$Annotation <- factor(
    model_results$Annotation,
    levels = rev(c(
      "MAIT cells", "TRAV1-2- CD8 T cells", "CD4 T cells", 
      "gd T cells", "DN T cells", "ILCs", "NK cells",
      "Myeloid cells", "B cells", "Progenitor cells",
      "CD4/CD8 T cells", "CD4/TRAV1-2- CD8 T cells"
    ))
  )
  
  model_results
}

plot_age_association_summary <- function(model_results, celltype_colors_major) {
  p_box <- ggplot(model_results, aes(y = Annotation, x = r, fill = Annotation)) +
    geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.6) +
    geom_point(position = position_jitterdodge(jitter.width = 0), size = 1, alpha = 0.8) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.4) +
    scale_x_continuous(name = "Correlation", limits = c(-0.5, 0.5)) +
    scale_fill_manual(values = celltype_colors_major) +
    theme_bw(base_size = 10) +
    theme(
      legend.position = "none",
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    ) +
    coord_fixed(ratio = 0.18)
  
  heatmap_data <- model_results %>%
    mutate(r_plot = ifelse(sig == "ns", 0, r))
  
  p_heatmap <- ggplot(heatmap_data, aes(x = Dataset, y = Annotation, fill = r_plot)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sig), size = 3) +
    scale_fill_gradient2(
      low = "#4575b4",
      mid = "white",
      high = "#d73027",
      midpoint = 0,
      name = "Pearson cor",
      limits = c(-0.5, 0.5),
      na.value = "white"
    ) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
      axis.text.y = element_text(color = "black"),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      legend.position = "top",
      legend.key.size = unit(0.5, "cm"),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7)
    ) +
    coord_fixed(ratio = 0.9)
  
  p_heatmap + p_box
}

df_age <- prepare_age_association_df(data = data, major_percentages = major_percentages)
model_results <- calculate_age_associations(df_age)
plot_age_association_summary(model_results, celltype_colors_major)

### Bar plot & Box plot major populations ###

plot_major_percent_boxplot <- function(
    major_percentages,
    annotation_levels = c("CD4 T", "CD8 T", "MAIT", "gd T",
                          "DN T", "ILCs", "NK", "Myeloid",
                          "B", "Progenitor")
) {
  
  df <- major_percentages %>%
    mutate(
      Dataset = factor(Dataset, levels = dataset_levels),
      Annotation = factor(Annotation, levels = annotation_levels)
    )
  
  donor_n_text <- df %>%
    distinct(Dataset, donor_id) %>%
    count(Dataset) %>%
    mutate(label = paste0("n = ", n))
  
  donor_caption <- paste(
    donor_n_text$Dataset,
    donor_n_text$label,
    sep = ": ",
    collapse = "   "
  )
  
  ggplot(df, aes(x = Dataset, y = percent, fill = Dataset, color = Dataset)) +
    geom_jitter(width = 0.15, size = 0.2, alpha = 0.8) +
    stat_boxplot(geom = "errorbar", width = 0.4) +
    geom_boxplot(width = 0.5, outlier.shape = NA, color = "black", alpha = 0.6) +
    scale_fill_manual(values = dataset_colors) +
    scale_color_manual(values = dataset_colors) +
    labs(y = "% of CD45+ cells", x = NULL, caption = donor_caption) +
    facet_wrap(~Annotation, scales = "free", ncol = 5) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "top",
      axis.ticks.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey40", linewidth = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black")
    )
}


plot_major_percent_alluvial <- function(
    major_percentages,
    annotation_levels = c("CD4 T", "CD8 T", "MAIT", "gd T",
                          "DN T", "ILCs", "NK", "Myeloid",
                          "B", "Progenitor")
) {
  
  df_plot <- major_percentages %>%
    mutate(
      Dataset = factor(Dataset, levels = dataset_levels),
      Annotation = factor(Annotation, levels = annotation_levels)
    ) %>%
    group_by(Dataset, Annotation) %>%
    summarise(mean_percent = mean(percent), .groups = "drop")
  
  ggplot(
    df_plot,
    aes(
      x = Dataset,
      stratum = Annotation,
      alluvium = Annotation,
      y = mean_percent,
      fill = Annotation,
      label = Annotation
    )
  ) +
    geom_flow(
      stat = "alluvium",
      lode.guidance = "forward",
      color = "black",
      alpha = 0.3,
      width = 0.7
    ) +
    geom_stratum(width = 0.7, color = "black", alpha = 0.6) +
    scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
    scale_fill_manual(values = celltype_colors_major) +
    xlab("Dataset") +
    ylab("") +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "right",
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
      panel.grid.major.y = element_line(color = "gray70")
    )
}

plot_major_percent_boxplot(major_percentages)
plot_major_percent_alluvial(major_percentages)

### Scatter plots ###

calculate_age_associations_minor <- function(minor_percentages) {
  
  minor_percentages %>%
    group_by(Major, Dataset, Annotation) %>%
    do({
      subdf <- .
      
      r <- if (length(unique(subdf$age)) > 1 &&
               length(unique(subdf$percent)) > 1) {
        cor(subdf$age, subdf$percent, method = "pearson", use = "complete.obs")
      } else {
        NA
      }
      
      model <- if (n_distinct(subdf$Ethnicity) > 1) {
        lm(percent ~ age + sex + Ethnicity, data = subdf)
      } else {
        lm(percent ~ age + sex, data = subdf)
      }
      
      out <- broom::tidy(model)
      out$r <- r
      out
    }) %>%
    filter(term == "age") %>%
    ungroup() %>%
    group_by(Major, Dataset) %>%
    mutate(p.adj = p.adjust(p.value, method = "BH")) %>%
    ungroup()
}
donor_meta = data[, c('donor_id', 'age', 'sex', 'ethnicity')]
donor_meta = donor_meta[!duplicated(donor_meta),]
minor_percentages = merge(minor_percentages, donor_meta, by='donor_id')
model_minor_results <- calculate_age_associations_minor(minor_percentages)

plot_minor_age_scatter <- function(
    minor_percentages,
    model_results,
    major = "CD8 T",
    annotations = c("RTE", "Naive", "Tcm CCR4+", "Tem GZMK+", 
                    "Tem GZMB+", "Tmem KLRC2+", "HLA-DR+")){
  
  plot_list <- list()
  
  for (subset in annotations) {
    
    df_plot <- minor_percentages %>%
      filter(Major == major, Annotation == subset) %>%
      left_join(
        model_results %>% 
          filter(Major == major) %>%
          select(Annotation, Dataset, p.adj),
        by = c("Annotation", "Dataset")
      ) %>%
      mutate(
        Dataset = factor(Dataset, levels = dataset_levels),
        p_label = ifelse(
          is.na(p.adj),
          "adj. p = NA",
          paste0("adj. p = ", signif(p.adj, 2))
        )
      )
    
    labels_df <- df_plot %>%
      group_by(Dataset) %>%
      slice(1) %>%
      ungroup()
    
    plot_list[[subset]] <- ggplot(df_plot, aes(x = age, y = percent, fill = Dataset)) +
      geom_point(alpha = 0.5, shape = 21) +
      geom_smooth(method = "lm", se = FALSE, color = "black") +
      facet_wrap(~ Dataset, scales = "free_y", ncol = 4) +
      geom_text(
        data = labels_df,
        aes(x = -Inf, y = Inf, label = p_label),
        hjust = -0.1,
        vjust = 1.1,
        inherit.aes = FALSE,
        size = 3.5
      ) +
      scale_fill_manual(values = dataset_colors) +
      theme_bw() +
      labs(
        x = "Age",
        y = paste0("% from ", major)
      ) +
      theme(
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)
      ) +
      ggtitle(subset)
  }
  
  plot_list
}

plot_list <- plot_minor_age_scatter(minor_percentages = minor_percentages, model_results = model_minor_results, major = "CD8 T")
plot_list$`Tem GZMK+`


### Oldest old ###

meta <- data %>%
  select(donor_id, sex, age, ethnicity, Dataset) %>%
  distinct()

donors_keep <- meta %>%
  mutate(Age_group = case_when(
    age < 35 ~ "< 35",
    age >= 60 & age <= 70 & Dataset %in% c("OneK1K", "SPAC", "KCL") ~ "60-70",
    age >= 85 ~ "85-97",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Age_group))

percentages_age <- minor_percentages %>%
  filter(donor_id %in% donors_keep$donor_id) %>%
  left_join(donors_keep %>% select(donor_id, Age_group), by = "donor_id") %>%
  mutate(Age_group = factor(Age_group, levels = c("< 35", "60-70", "85-97")))

dunn_results <- percentages_age %>%
  group_by(Annotation) %>%
  rstatix::dunn_test(percent ~ Age_group, p.adjust.method = "bonferroni") %>%
  filter(group1 %in% c("< 35", "60-70") & group2 %in% c("60-70", "85-97")) %>%
  mutate(y.position = 12.5)

ggplot(
  percentages_age, aes(x = Age_group, y = percent, fill = Age_group, color = Age_group)) +
  geom_jitter(width = 0.2, alpha = 0.8, size = 1) +
  geom_boxplot(outlier.shape = NA, width = 0.5, color = "black") +
  facet_wrap(~Annotation, scales = "free") +
  theme_minimal() +
  scale_fill_manual(values = c("< 35" = "#d9b782", "60-70" = "#8a9167", "85-97" = "#a96f98")) +
  scale_color_manual(values = c("< 35" = "#d9b782", "60-70" = "#8a9167", "85-97" = "#a96f98"))

### Selected sex & ancestry ###

percentages <- minor_percentages %>%
  filter(ethnicity %in% c("Asian", "European"))

percentages$Annotation = paste(percentages$Major, percentages$Annotation, sep=' ')

ethnicity_percentages <- percentages %>%
  filter(Annotation %in% c("CD4 T RTE", "CD8 T Tem GZMB+", "NK CD56bright", "CD8 T Trm"))

ggplot(ethnicity_percentages, aes(x = ethnicity, y = percent, fill = ethnicity, color = ethnicity)) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  geom_boxplot(outlier.shape = NA, width = 0.5, color = "black") +
  scale_fill_manual(values = c("Asian" = "#96af7e", "European" = "#bfbd7d")) +
  scale_color_manual(values = c("Asian" = "#96af7e", "European" = "#bfbd7d")) +
  facet_wrap(~Annotation, scales = "free") +
  theme_minimal() +
  labs(x = NULL, y = "%") +
  theme(axis.text.x = element_text(color = "black"), legend.position = "none")

ethnicity_age <- ethnicity_percentages %>%
  filter(age < 35 | age >= 55) %>%
  mutate(
    Age_group = ifelse(age < 35, "< 35", "55+"),
    Condition = paste(Age_group, ethnicity)
  )

ggplot(ethnicity_age, aes(x = Condition, y = percent, fill = Condition, color = Condition)) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  geom_boxplot(outlier.shape = NA, width = 0.5, color = "black") +
  scale_fill_manual(values = c("#96af7e", "#a699ae", "#bfbd7d", "#a96f98")) +
  scale_color_manual(values = c("#96af7e", "#a699ae", "#bfbd7d", "#a96f98")) +
  facet_wrap(~Annotation, scales = "free") +
  theme_minimal() +
  labs(x = NULL, y = "%") +
  theme(axis.text.x = element_text(color = "black"), legend.position = "none")


sex_percentages <- percentages %>%
  filter(Annotation %in% c("CD4 T Th2", "CD8 T RTE", "NK CD56bright", "CD4 Naive"))

ggplot(sex_percentages, aes(x = sex, y = percent, fill = sex, color = sex)) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  geom_boxplot(outlier.shape = NA, width = 0.5, color = "black") +
  scale_fill_manual(values = c("Female" = "#6794a0", "Male" = "#af7187")) +
  scale_color_manual(values = c("Female" = "#6794a0", "Male" = "#af7187")) +
  facet_wrap(~Annotation, scales = "free") +
  theme_minimal() +
  labs(x = NULL, y = "%") +
  theme(axis.text.x = element_text(color = "black"), legend.position = "none")

ethnicity_age <- ethnicity_percentages %>%
  filter(age < 35 | age >= 55) %>%
  mutate(
    Age_group = ifelse(age < 35, "< 35", "55+"),
    Condition = paste0(Age_group, ethnicity)
  )

ggplot(ethnicity_age, aes(x = Condition, y = percent, fill = Condition, color = Condition)) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  geom_boxplot(outlier.shape = NA, width = 0.5, color = "black") +
  scale_fill_manual(values = c("#96af7e", "#a699ae", "#bfbd7d", "#a96f98")) +
  scale_color_manual(values = c("#96af7e", "#a699ae", "#bfbd7d", "#a96f98")) +
  facet_wrap(~Annotation, scales = "free") +
  theme_minimal() +
  labs(x = NULL, y = "%") +
  theme(axis.text.x = element_text(color = "black"), legend.position = "none")

### Variance partition ###

data <- data %>%
  filter(ethnicity %in% c("Asian", "European"))

meta <- data %>%
  select(donor_id, sex, age, ethnicity, Dataset) %>%
  distinct() %>%
  rename(Sex = sex, Age = age, Ethnicity = ethnicity)

expr_df <- minor_percentages %>%
  filter(donor_id %in% meta$donor_id) %>%
  mutate(
    Annotation = paste(Major, Annotation, sep = "_"),
    Annotation = str_replace_all(Annotation, " T cells_", "_"),
    Annotation = str_replace_all(Annotation, " cells_", "_")
  ) %>%
  select(Annotation, donor_id, percent) %>%
  pivot_wider(names_from = donor_id, values_from = percent, values_fill = 0)

expr <- as.data.frame(expr_df)
rownames(expr) <- expr$Annotation
expr$Annotation <- NULL

rownames(meta) <- meta$donor_id
expr <- expr[, rownames(meta)]
expr[is.na(expr)] <- 0

form <- ~ (1 | Sex) + (1 | Ethnicity) + (1 | Dataset) + Age

varPart <- fitExtractVarPartModel(expr, form, meta)

vp <- sortCols(varPart)
plotVarPart(vp)

vp_long <- as.data.frame(varPart) %>%
  rownames_to_column("population") %>%
  pivot_longer(
    cols = c("Age", "Sex", "Ethnicity"),
    names_to = "Variable",
    values_to = "Variance"
  ) %>%
  mutate(
    Variance = Variance * 100,
    Cluster = case_when(
      str_detect(population, "^CD4 T_") ~ "CD4",
      str_detect(population, "^CD8 T_") ~ "CD8",
      str_detect(population, "^gd T_") ~ "gd",
      str_detect(population, "^NK_") ~ "NK",
      str_detect(population, "^Myeloid_") ~ "Myeloid",
      str_detect(population, "^B_") ~ "B",
      TRUE ~ "Other"
    ),
    CleanName = population %>%
      str_remove("^CD4 T_") %>%
      str_remove("^CD8 T_") %>%
      str_remove("^gd T_") %>%
      str_remove("^NK_") %>%
      str_remove("^Myeloid_") %>%
      str_remove("^B_"),
    CleanName = ifelse(population == "B_CD5+ B cells", "CD5+ B cells", paste(Cluster, CleanName))
  )

final_order <- vp_long %>%
  filter(Variable == "Age") %>%
  group_by(CleanName) %>%
  summarise(mean_var = mean(Variance, na.rm = TRUE), .groups = "drop") %>%
  arrange(mean_var) %>%
  pull(CleanName)

vp_bar <- vp_long %>%
  filter(Variance >= 1) %>%
  mutate(
    CleanName = factor(CleanName, levels = final_order),
    Variance_label = round(Variance, 1)
  )

ggplot(vp_bar, aes(x = Variance, y = CleanName, fill = Cluster)) +
  geom_col(width = 0.85, color = "black", alpha = 0.7) +
  geom_text(aes(label = Variance_label), hjust = -0.2, size = 3) +
  facet_grid(. ~ Variable, scales = "free_x", space = "free") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.25))) +
  scale_fill_manual(values = c(
    "B" = "#ff8899",
    "CD4" = "#eb883c",
    "CD8" = "#298283",
    "gd" = "#6f7dc6",
    "Myeloid" = "#94755e",
    "NK" = "#9d6db8"
  )) +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 8),
    axis.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(fill = "Cluster")

# -----------------------------
# Figure 3, S3
# -----------------------------

cd8 <- c("RTE", "Naive", "Naive-IFN", "Tcm CCR4-", "Tcm CCR4+",
         "Trm", "Tem GZMK+", "Tem GZMB+", "Temra", "NKT-like",
         "Tmem KLRC2+", "HLA-DR+", "Proliferative")

cd8_mat <- minor_percentages %>%
  filter(Major == "CD8 T", Annotation %in% cd8) %>%
  select(donor_id, Annotation, percent) %>%
  pivot_wider(names_from = Annotation, values_from = percent, values_fill = 0)

t_percentages <- cd8_mat %>%
  column_to_rownames("donor_id")

pca <- prcomp(t_percentages, scale. = FALSE, center = TRUE)
donor_meta = data[, c('donor_id', 'sex', 'age', 'ethnicity', 'Dataset')]
donor_meta = donor_meta[!duplicated(donor_meta),]

pca_df <- as.data.frame(pca$x) %>%
  rownames_to_column("donor_id") %>%
  left_join(
    data %>% select(donor_id, sex, age, ethnicity, Dataset, CMV) %>% distinct(),
    by = "donor_id"
  ) %>%
  left_join(cd8_mat, by = "donor_id")

ggplot(pca_df, aes(PC1, PC2, color = age)) +
  geom_point(size = 1.5, alpha = 0.8) +
  theme_minimal() +
  labs(
    x = paste0("PC1 (", round(summary(pca)$importance[2, 1] * 100, 1), "% variance)"),
    y = paste0("PC2 (", round(summary(pca)$importance[2, 2] * 100, 1), "% variance)")
  ) +
  theme(legend.position = "bottom") +
  scale_color_gradientn(colors = c("#5aa353", "#f09618", "#cc0418"))

ggplot(pca_df, aes(PC1, PC2, color = Naive)) +
  geom_point(size = 1.5, alpha = 0.8) +
  theme_minimal() +
  labs(
    x = paste0("PC1 (", round(summary(pca)$importance[2, 1] * 100, 1), "% variance)"),
    y = paste0("PC2 (", round(summary(pca)$importance[2, 2] * 100, 1), "% variance)")
  ) +
  theme(legend.position = "bottom") +
  scale_color_gradient2(low = "#569cb4", mid = "#b88084", high = "#ff001e", midpoint = 40)

fviz_eig(pca) +
  ggtitle("% of CD8 T cell populations for each individual")

fviz_pca_biplot(
  pca,
  geom = "point",
  col.var = "#1F78B4",
  col.ind = "#bdbbbb"
) +
  ggtitle("Loadings of CD8 T cell subpopulations")

fviz_contrib(pca, choice = "var", axes = 2, top = 5) +
  ggtitle("Contribution of variables to PC2")


global_hull <- pca_df %>%
  slice(chull(PC1, PC2))

global_hull_expanded <- bind_rows(
  lapply(unique(pca_df$Dataset), function(ds) {
    global_hull %>% mutate(Dataset = ds)
  })
)

pca_df <- pca_df %>% mutate(emK_emB = `Tem GZMK+` / `Tem GZMB+`)
pca_df[is.infinite(pca_df$emK_emB), 'emK_emB'] = max(pca_df[!is.infinite(pca_df$emK_emB), 'emK_emB'])

ggplot(pca_df, aes(x = PC1 * -1, y = PC2)) +
  geom_polygon(
    data = global_hull_expanded,
    aes(group = Dataset),
    fill = NA,
    color = "#696969",
    linewidth = 0.6
  ) +
  geom_point(aes(color = log1p(emK_emB)), size = 1.5, alpha = 0.8) +
  scale_color_gradient2(low = "#dadfe0", mid = "#f27466", high = "#c71602", midpoint = 3) +
  facet_wrap(~Dataset, ncol = 4) +
  labs(
    x = paste0("PC1 (", round(summary(pca)$importance[2, 1] * 100, 1), "% variance)"),
    y = paste0("PC2 (", round(summary(pca)$importance[2, 2] * 100, 1), "% variance)")
  ) +
  theme_minimal()

model_cd8_summary <- model_minor_results %>%
  filter(Major == "CD8 T") %>%
  group_by(Annotation) %>%
  summarise(
    mean_r = mean(r, na.rm = TRUE),
    mean_logp = mean(-log10(p.adj), na.rm = TRUE),
    n_sig = sum(p.adj < 0.05, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(model_cd8_summary, aes(x = mean_r, y = mean_logp)) +
  geom_point(aes(size = n_sig, color = mean_r), alpha = 0.8) +
  geom_text_repel(
    data = model_cd8_summary %>% filter(n_sig > 0),
    aes(label = Annotation),
    size = 3
  ) +
  scale_color_gradient2(low = "blue", mid = "grey85", high = "red", midpoint = 0) +
  scale_size_continuous(name = "N significant datasets", range = c(1.5, 6)) +
  theme_bw() +
  labs(
    x = "Mean Pearson r across datasets",
    y = "Mean -log10(adj. p)"
  )


# -----------------------------
# Figure 4, S4
# -----------------------------

cd8_ratio <- minor_percentages %>%
  filter(
    Major == "CD8 T",
    Annotation %in% c("Tem GZMK+", "Tem GZMB+")
  ) %>%
  select(donor_id, Annotation, n_cells, percent) %>%
  pivot_wider(
    names_from = Annotation,
    values_from = c(n_cells, percent),
    values_fill = 0
  ) %>%
  left_join(
    data %>%
      select(donor_id, sex, BMI, CMV, age, ethnicity, Dataset) %>%
      distinct(),
    by = "donor_id"
  ) %>%
  mutate(
    Ratio = `n_cells_Tem GZMK+` / `n_cells_Tem GZMB+`,
    Ratio = ifelse(is.infinite(Ratio),
                   max(Ratio[is.finite(Ratio)], na.rm = TRUE),
                   Ratio)
  )

ggplot(cd8_ratio, aes(x = sex, y = log1p(Ratio), fill = sex, color = sex)) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  geom_boxplot(outlier.shape = NA, width = 0.5, color = "black") +
  scale_fill_manual(values = c("#7293A0", "#A37986")) +
  scale_color_manual(values = c("#7293A0", "#A37986")) +
  xlab("") +
  ylab("log1p(Tem GZMK+/Tem GZMB+)") +
  theme_minimal() +
  theme(axis.text.x = element_text(color = "black"),
        legend.position = "none")

ratio_ethnicity <- cd8_ratio %>%
  filter(ethnicity %in% c("Asian", "European"))

ggplot(ratio_ethnicity, aes(x = ethnicity, y = log1p(Ratio), fill = ethnicity, color = ethnicity)) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  geom_boxplot(outlier.shape = NA, width = 0.5, color = "black") +
  scale_fill_manual(values = c("#96af7e", "#bfbd7d")) +
  scale_color_manual(values = c("#96af7e", "#bfbd7d")) +
  xlab("") +
  ylab("log1p(Tem GZMK+/Tem GZMB+)") +
  theme_minimal() +
  theme(axis.text.x = element_text(color = "black"),
        legend.position = "none")

ratio_bmi <- cd8_ratio %>%
  filter(BMI != "No_info") %>%
  mutate(
    BMI = as.numeric(gsub(".0", "", BMI, fixed = TRUE)),
    BMI_category = case_when(
      BMI < 18.5 ~ "Underweight",
      BMI < 25 ~ "Normal",
      BMI < 30 ~ "Overweight",
      BMI >= 30 ~ "Obese"
    )
  )

ratio_bmi$BMI_category = factor(ratio_bmi$BMI_category, levels = c( "Underweight", "Normal", "Overweight", "Obese"))
ggplot(ratio_bmi, aes(x = BMI_category, y = log1p(Ratio),
                      fill = BMI_category, color = BMI_category)) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  geom_boxplot(outlier.shape = NA, width = 0.5, color = "black") +
  scale_fill_manual(values = c("#b0c4b1", "#88a096", "#d9b382", "#a65e3f")) +
  scale_color_manual(values = c("#b0c4b1", "#88a096", "#d9b382", "#a65e3f")) +
  xlab("") +
  ylab("log1p(Tem GZMK+/Tem GZMB+)") +
  theme_minimal() +
  theme(axis.text.x = element_text(color = "black"),
        legend.position = "none")


ratio_age <- cd8_ratio %>%
  mutate(
    Age_group = case_when(
      age < 35 ~ "A",
      age < 45 ~ "B",
      age < 55 ~ "C",
      age < 65 ~ "D",
      age < 75 ~ "E",
      age >= 75 ~ "F"
    )
  )

ratio_age$Age_group = factor(ratio_age$Age_group, levels = c( "A", "B", "C", "D", "E", "F"))
ggplot(ratio_age, aes(x = Age_group, y = log1p(Ratio),
                      fill = Age_group, color = Age_group)) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  geom_boxplot(outlier.shape = NA, width = 0.5, color = "black") +
  scale_fill_manual(values = c(
    "#d4b483", "#a29dae", "#7d8ca3",
    "#b2836d", "#8e8c68", "#9c7a97"
  )) +
  scale_color_manual(values = c(
    "#d4b483", "#a29dae", "#7d8ca3",
    "#b2836d", "#8e8c68", "#9c7a97"
  )) +
  xlab("") +
  ylab("log1p(Tem GZMK+/Tem GZMB+)") +
  theme_minimal() +
  theme(axis.text.x = element_text(color = "black"),
        legend.position = "none")

ratio_cmv <- cd8_ratio %>%
  filter(CMV != "No_info")

ggplot(ratio_cmv, aes(x = CMV, y = log1p(Ratio), fill = CMV)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0, alpha = 0.7, size = 1.5,
              shape = 21, color = "#706e6e") +
  scale_fill_manual(values = c("#8fb996", "#c26a6a")) +
  xlab("") +
  ylab("log1p(Tem GZMK+/Tem GZMB+)") +
  theme_minimal() +
  theme(axis.text.x = element_text(color = "black"),
        legend.position = "none") + facet_wrap(~Dataset)

global_hull <- pca_df %>%
  slice(chull(PC1, PC2))

ggplot() +
  geom_point(data = pca_df %>% filter(CMV == 'No_info'), aes(x = PC1, y = PC2), color='#d4d4d4', size = 1.5) +
  geom_point(data = pca_df %>% filter(CMV != 'No_info'), aes(x = PC1, y = PC2, color=CMV), size = 1.5) +
  theme_minimal() +
  labs(
    x = paste0("PC1 (", round(summary(pca)$importance[2, 1] * 100, 1), "% variance)"),
    y = paste0("PC2 (", round(summary(pca)$importance[2, 2] * 100, 1), "% variance)")
  ) +
  theme(legend.position = "bottom")+
  scale_color_manual(values =   c('#8fb996',"#c26a6a")) +
  geom_polygon(data = global_hull, aes(x = PC1, y = PC2), fill = NA, color = "#696969", linewidth = 0.6)

cmv_df <- minor_percentages %>%
  filter(Dataset %in% c("SoundLife", "KCL")) %>%
  left_join(
    data %>% select(donor_id, CMV) %>% distinct(),
    by = "donor_id"
  ) %>%
  filter(CMV %in% c("Negative", "Positive")) %>%
  mutate(
    Age_group = case_when(
      age < 40 ~ "Young",
      age > 50 ~ "Old",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Age_group)) %>% filter(age < 75 & CMV != 'No_info') 

cmv_df$Major_Annotation = paste(cmv_df$Major, cmv_df$Annotation)

result <- cmv_df %>%
  group_by(Dataset, Age_group, Major_Annotation) %>%
  wilcox_test(percent ~ CMV) %>%
  mutate(p.adj = p * 55)

result_left <- result %>%
  group_by(Age_group, Major_Annotation) %>%
  summarise(
    n_sig = sum(p.adj < 0.05, na.rm = TRUE),
  ) %>%
  filter(n_sig >= 1)

result = result %>% filter(Major_Annotation %in% result_left$Major_Annotation)
result[result$p.adj > 1, 'p.adj'] = 1

plot_df <- result %>%
  mutate(
    neglog_padj = -log10(p.adj),
    neglog_padj = ifelse(is.infinite(neglog_padj), NA, neglog_padj)
  )

order_df <- plot_df %>%
  group_by(Major_Annotation) %>%
  summarise(max_p = max(neglog_padj, na.rm = TRUE), .groups = "drop")

plot_df <- plot_df %>%
  left_join(order_df, by = "Major_Annotation") %>%
  mutate(
    Major_Annotation = fct_reorder(Major_Annotation, max_p),
    Age_group = factor(Age_group, levels = c("Young", "Old"))
  )

bar_df <- plot_df %>%
  group_by(Major_Annotation, Age_group) %>%
  summarise(
    mean_neglog_padj = mean(neglog_padj, na.rm = TRUE),
    sem = sd(neglog_padj, na.rm = TRUE) / sqrt(sum(!is.na(neglog_padj))),
    upper = mean_neglog_padj + sem,
    .groups = "drop"
  )

ggplot(bar_df, aes(y = Major_Annotation, group = Age_group)) +
  geom_col(
    aes(x = mean_neglog_padj, fill = Age_group),
    position = position_dodge(width = 0.75),
    width = 0.65,
    color = "black",
    alpha = 0.8
  ) +
  geom_errorbarh(
    aes(xmin = mean_neglog_padj, xmax = upper),
    position = position_dodge(width = 0.75),
    height = 0.18,
    color = "black"
  ) +
  geom_point(
    data = plot_df,
    aes(x = neglog_padj, y = Major_Annotation, group = Age_group, fill=Age_group),
    position = position_dodge(width = 0.75),
    size = 2,
    alpha = 0.8,
    inherit.aes = FALSE, shape=21
  ) +
  scale_fill_manual(values = c("Young" = "#D9A75F", "Old" = "#8FA4BF")) +
  labs(
    title = "CMV- vs CMV+ Population Changes",
    x = expression(-log[10]("p.adj")),
    y = NULL,
    fill = "Age group"
  ) +
  theme_bw()


### Add validation Aquino et al. ###

aquino_percentages = read.csv('Aquino_immune_percentages.csv')
aquino_percentages = aquino_percentages %>% filter(Major_Annotation %in% result$Major_Annotation)

old = aquino_percentages %>% filter(Age > 45)
young = aquino_percentages %>% filter(Age < 46)

result_old <- old %>%
  group_by(Major_Annotation) %>%
  wilcox_test(percent ~ CMV) %>%
  adjust_pvalue(method = "bonferroni") %>%
  mutate(Age_group = "Old")

result_young <- young %>%
  group_by(Major_Annotation) %>%
  wilcox_test(percent ~ CMV) %>%
  adjust_pvalue(method = "bonferroni") %>%
  mutate(Age_group = "Young")

result_aquino <- bind_rows(
  result_old %>% select(Major_Annotation, p.adj, Age_group),
  result_young %>% select(Major_Annotation, p.adj, Age_group)
)
result_aquino$Dataset = 'Aquino'

result$Major_Annotation = result$Major_Annotation
result <- rbind(result[, colnames(result_aquino)], result_aquino)

plot_df <- result %>%
  mutate(
    neglog_padj = -log10(p.adj),
    neglog_padj = ifelse(is.infinite(neglog_padj), NA, neglog_padj)
  )

order_df <- plot_df %>%
  group_by(Major_Annotation) %>%
  summarise(max_p = max(neglog_padj, na.rm = TRUE), .groups = "drop")

plot_df <- plot_df %>%
  left_join(order_df, by = "Major_Annotation") %>%
  mutate(
    Major_Annotation = fct_reorder(Major_Annotation, max_p),
    Age_group = factor(Age_group, levels = c("Young", "Old"))
  )

bar_df <- plot_df %>%
  group_by(Major_Annotation, Age_group) %>%
  summarise(
    mean_neglog_padj = mean(neglog_padj, na.rm = TRUE),
    sem = sd(neglog_padj, na.rm = TRUE) / sqrt(sum(!is.na(neglog_padj))),
    upper = mean_neglog_padj + sem,
    .groups = "drop"
  )
bar_df$Age_group = factor(bar_df$Age_group, levels = c('Old', 'Young'))
plot_df$Age_group = factor(plot_df$Age_group, levels = c('Old', 'Young'))
ggplot(bar_df, aes(y = Major_Annotation, group = Age_group)) +
  geom_col(
    aes(x = mean_neglog_padj, fill = Age_group),
    position = position_dodge(width = 0.75),
    width = 0.65,
    color = "black"
  ) +
  geom_errorbarh(
    aes(xmin = mean_neglog_padj, xmax = upper),
    position = position_dodge(width = 0.75),
    height = 0.18,
    color = "black"
  ) +
  geom_point(
    data = plot_df,
    aes(x = neglog_padj, y = Major_Annotation, group = Age_group, fill=Age_group),
    position = position_dodge(width = 0.75),
    size = 1.5,
    inherit.aes = FALSE, shape=21
  ) +
  scale_fill_manual(values = c("Young" = "#D9A75F", "Old" = "#8FA4BF")) +
  labs(
    x = "-log10(p.adj)",
    y = NULL,
    fill = "Age group"
  ) +
  theme_bw()

ggplot(cmv_df %>% filter(Dataset == 'KCL' & Major_Annotation %in% result$Major_Annotation & Age_group == 'Young'), 
       aes(x = CMV, y = percent, fill = CMV, color = CMV)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, color = "black") +
  geom_jitter(width = 0.05, alpha = 0.7, size = 1.5, shape = 21) +
  facet_grid(Age_group ~ Major_Annotation, scales = "free") +
  scale_fill_manual(values = c("Negative" = "#D4B483", "Positive" = "#A29DAE")) +
  scale_color_manual(values = c("Negative" = "#a68e68", "Positive" = "#6b647d")) +
  theme_minimal() +
  xlab("") +
  ylab("% from corresponding major population") +
  theme(
    axis.text.x = element_text(color = "black"),
    legend.position = "none"
  ) + facet_wrap(~Major_Annotation, scales = 'free_y', ncol = 5)


# -----------------------------
# Figure 5, S5
# -----------------------------

df <- read.csv('Olink_proteins_all_cohorts.csv') 

results <- df %>%
  group_by(Dataset, PROTEIN) %>%
  do({
    data_subset <- .
    cor_test <- cor.test(data_subset$pheb, data_subset$NPX, method = "pearson")
    lm_fit <- lm(pheb ~ NPX + age + sex, data = data_subset)
    lm_summary <- tidy(lm_fit) %>% filter(term == "NPX")
    tibble(
      pearson_R = cor_test$estimate,
      pearson_p = cor_test$p.value,
      beta_adj  = lm_summary$estimate,
      p_adj     = lm_summary$p.value
    )
  }) %>%
  ungroup()

results <- results %>%
  group_by(Dataset) %>%
  mutate(signif = p_adj < 0.05) %>%
  ungroup()

results <- results %>%
  mutate(rank_metric = abs(pearson_R) * -log10(p_adj))

top_labels <- results %>%
  group_by(Dataset) %>%
  arrange(desc(rank_metric)) %>%
  filter(pearson_R > 0) %>%
  slice_head(n = 10) %>%
  bind_rows(
    results %>%
      group_by(Dataset) %>%
      arrange(desc(rank_metric)) %>%
      filter(pearson_R < 0) %>%
      slice_head(n = 10)
  ) %>%
  ungroup()

ggplot(results, aes(x = pearson_R, y = -log10(p_adj))) +
  geom_point(aes(fill = ifelse(signif, pearson_R, NA), color = ifelse(signif, pearson_R, NA)), 
             shape = 21, size = 2, alpha = 0.8) + 
  scale_fill_gradient(low="blue", high="red", na.value = 'darkgrey') +
  scale_color_gradient(low="black", high="black", na.value = 'darkgrey') +
  geom_text_repel(data = top_labels,
                  aes(label = PROTEIN),
                  size = 3,
                  max.overlaps = Inf) +
  facet_wrap(~Dataset) +
  theme_minimal(base_size = 10) +
  labs(x = "R",
       y = "-log10(p-value)") +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90"),
    legend.position = 'none'
  )

