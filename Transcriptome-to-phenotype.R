library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(pbapply)

# -------------------------
# Load and combine expression data (ignore sex)
# -------------------------
exp_female <- read_table("dgrp.array.exp.female.txt") %>%
  rename_with(~ str_replace(.x, "^line_", "DGRP"))

exp_male <- read_table("dgrp.array.exp.male.txt") %>%
  rename_with(~ str_replace(.x, "^line_", "DGRP"))

expression_matrix <- bind_rows(exp_female, exp_male)

expression_matrix_means <- expression_matrix %>%
  pivot_longer(cols = -gene, names_to = c("line", "rep"), names_sep = ":") %>%
  group_by(gene, line) %>%
  summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = line, values_from = mean_value)

expr_df <- as.data.frame(t(expression_matrix_means[-1]))
colnames(expr_df) <- expression_matrix_means$gene
rownames(expr_df) <- colnames(expression_matrix_means)[-1]
scaled_expr <- scale(expr_df)

# -------------------------
# Load phenotype data (raw, individual-level, both sexes)
# -------------------------
phenotype_matrix <- read_table("PhenoFM_trimmed.txt") %>%
  rename(line_id = `#IID`, phenotype = median_survival) %>%   # adjust colname if needed
  filter(phenotype != -9) %>%
  column_to_rownames("line_id")

# -------------------------
# Match lines
# -------------------------
common_lines <- intersect(rownames(scaled_expr), rownames(phenotype_matrix))
scaled_expr <- scaled_expr[common_lines, ]
phenotype_matrix <- phenotype_matrix[common_lines, , drop = FALSE]

# -------------------------
# Prepare phenotype vector
# -------------------------
y <- phenotype_matrix$phenotype

# -------------------------
# Simple linear model: phenotype ~ gene expression only
# -------------------------
run_simple_lm <- function(gene_name, y_vec, scaled_expr) {
  gene_expr <- scaled_expr[, gene_name]
  fit <- lm(y_vec ~ gene_expr)
  
  summary_fit <- summary(fit)
  
  data.frame(
    gene = gene_name,
    Beta = coef(summary_fit)["gene_expr", "Estimate"],
    SE = coef(summary_fit)["gene_expr", "Std. Error"],
    t_value = coef(summary_fit)["gene_expr", "t value"],
    P_value = coef(summary_fit)["gene_expr", "Pr(>|t|)"]
  )
}

gene_list <- colnames(scaled_expr)

simple_lm_results <- pbapply::pblapply(gene_list, function(g) {
  run_simple_lm(g, y, scaled_expr)
}) %>%
  bind_rows()

# -------------------------
# Save simple model results
# -------------------------
write_csv(simple_lm_results, "results_simple_lm_genes_nosex.csv")

print("Finished! Simple LM results saved.")


#### Manhattan plot

# -----------------------------
# Libraries
# -----------------------------
library(dplyr)
library(ggplot2)
library(readr)
library(rtracklayer)
library(patchwork)

# -----------------------------
# Load p-value files
# -----------------------------
pvalsF <- read_csv2("results_simple_lm_genes_nosex.csv", col_types = cols()) %>%
  mutate(P_value = as.numeric(P_value))

# -----------------------------
# Load GTF and extract gene positions
# -----------------------------
gtf <- rtracklayer::import("dm3_with_fb.gtf")
genes_gtf <- gtf[gtf$type == "gene"]

gene_positions <- data.frame(
  gene_id = genes_gtf$gene_id,
  gene_symbol = genes_gtf$gene_name,
  chr = as.character(seqnames(genes_gtf)),
  start = start(genes_gtf),
  end = end(genes_gtf)
)

# -----------------------------
# Prepare plot dataframe with cumulative positions
# -----------------------------
prepare_plot_df <- function(pvals, gene_positions) {
  df <- pvals %>%
    left_join(gene_positions, by = c("gene" = "gene_id")) %>%
    filter(!is.na(chr))
  
  # Main chromosomes
  main_chr <- c("2L","2R","3L","3R","4","X")
  
  # Collapse all others into "Others"
  df <- df %>%
    mutate(chr = ifelse(chr %in% main_chr, chr, "Others"))
  
  chr_order <- c(main_chr, "Others")
  df$chr <- factor(df$chr, levels = chr_order)
  
  # Compute chromosome lengths
  chr_lengths <- df %>%
    group_by(chr) %>%
    summarize(chr_len = max(end, na.rm = TRUE)) %>%
    arrange(factor(chr, levels = chr_order))
  
  # Cumulative positions
  chr_lengths <- chr_lengths %>%
    mutate(tot = cumsum(as.numeric(chr_len)) - chr_len)
  
  df <- df %>%
    left_join(chr_lengths[, c("chr", "tot")], by = "chr") %>%
    mutate(BPcum = start + tot,
           logP = -log10(P_value))
  
  return(df)
}

plotF <- prepare_plot_df(pvalsF, gene_positions)


# -----------------------------
# X-axis labels at chromosome centers
# -----------------------------
axisdf <- plotF %>%
  group_by(chr) %>%
  summarize(center = (min(BPcum) + max(BPcum)) / 2)

# -----------------------------
# Manhattan plot function
# -----------------------------
manhattan_plot <- function(df, title) {
  chr_levels <- levels(df$chr)
  
  # Chromosome shading
  xmin <- sapply(chr_levels, function(x) min(df$BPcum[df$chr == x]))
  xmax <- sapply(chr_levels, function(x) max(df$BPcum[df$chr == x]))
  shades <- data.frame(
    xmin = xmin, xmax = xmax,
    ymin = 0, ymax = Inf,
    fill = rep(c("#ffffff", "#ebebeb"), length.out = length(chr_levels))
  )
  
  # Bonferroni threshold
  N <- nrow(df)  # number of tests
  threshold <- 3
  
  ggplot() +
    geom_rect(data = shades, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill), alpha = 0.5) +
    scale_fill_identity() +
    geom_point(data = df, aes(x = BPcum, y = logP, color = chr), size = 1.5, alpha = 0.7) +
    # add this layer for red points
    geom_point(data = subset(df, logP >= 3),
               aes(x = BPcum, y = logP),
               color = "red", size = 1.5, alpha = 0.7) +
    geom_hline(yintercept = threshold, color = "red", linetype = "dashed", linewidth = 1) +
    scale_color_manual(values = rep(c("#4D4D4D", "#7FB3D5"), length.out = length(chr_levels))) +
    scale_x_continuous(breaks = axisdf$center, labels = axisdf$chr, expand = c(0, 0)) +
    scale_y_continuous(
      limits = c(0, max(df$logP, na.rm = TRUE) + 1),
      expand = c(0, 0)
    ) +
    labs(x = "Chromosome", y = "-log10(P-value)", title = title) +
    theme_minimal() +
    theme(
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  
}


# -----------------------------
# Plot side-by-side
# -----------------------------
manhattan_plot(plotF, "Mixed Expression")

# -----------------------------
# Optionally save
# -----------------------------
ggsave("manhattan_mixed.pdf", width = 14, height = 4)