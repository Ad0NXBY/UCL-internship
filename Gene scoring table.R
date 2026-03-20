scoring_df <- top_genes_per_locus %>%
  left_join(
    pops_preds %>%
      transmute(
        ENSGID = as.character(ENSGID),
        MAGMA_Y = as.numeric(Y)   
      ),
    by = "ENSGID"
  )

rescale_01 <- function(x) {
  r <- range(x, na.rm = TRUE)
  if (!is.finite(r[1]) || !is.finite(r[2]) || r[1] == r[2]) return(rep(0, length(x)))
  (x - r[1]) / (r[2] - r[1])
}

scoring_df <- scoring_df %>%
  group_by(locus_label) %>%
  mutate(
    pops_scaled_locus_01  = rescale_01(PoPS_score),
    magma_scaled_locus_01 = rescale_01(MAGMA_Y)
  ) %>%
  ungroup()

scoring_df <- scoring_df %>%
  left_join(
    comparison %>%
      select(locus_label, nearest_gene, pops_top_gene, same_gene),
    by = "locus_label"
  )



scoring_df <- scoring_df %>%
  mutate(
    is_nearest_gene = (GENE == nearest_gene),
    is_pops_top_gene = (GENE == pops_top_gene),
    
    # optional scores for these indicators
    nearest_gene_score = ifelse(is_nearest_gene, 1L, 0L),
    pops_top_gene_score = ifelse(is_pops_top_gene, 1L, 0L)
  )

scoring_df <- scoring_df %>%
  mutate(
    final_nearest_gene_score = ifelse(
      is_nearest_gene & is_pops_top_gene,
      1L, 0L
    )
  )


scoring_df <- scoring_df %>%
  mutate(
    combined_score = pops_scaled_locus_01 + magma_scaled_locus_01 + final_nearest_gene_score
  )

final_gene_scoring_table <- scoring_df %>%
  select(ENSGID, GENE, CHR, locus_label,
         pops_scaled_locus_01, magma_scaled_locus_01, 
         final_nearest_gene_score, combined_score) 


write.csv(
  final_gene_scoring_table,
  file = file.path(data_dir, "80k_gene_scoring_table.csv"),
  row.names = FALSE
)

#Final Scoring table heatmap====================================================
ncol_facets <- 4
nrow_facets <- 3

final_gene_scoring_table <- final_gene_scoring_table %>%
  mutate(
    locus_label = factor(
      locus_label,
      levels = unique(locus_label)
    )
  )

plot_df <- final_gene_scoring_table %>%
  select(locus_label, GENE,
         pops_scaled_locus_01,
         magma_scaled_locus_01,
         final_nearest_gene_score,
         combined_score) %>%
  mutate(
    locus_label = as.integer(locus_label),
    GENE = as.character(GENE)
  ) %>%
  pivot_longer(
    cols = c(pops_scaled_locus_01, magma_scaled_locus_01, final_nearest_gene_score),
    names_to = "metric",
    values_to = "score"
  ) %>%
  mutate(
    metric = recode(metric,
                    pops_scaled_locus_01     = "PoPS",
                    magma_scaled_locus_01    = "MAGMA",
                    final_nearest_gene_score = "Nearest gene"
    ),
    metric = factor(metric, levels = c("MAGMA", "Nearest gene", "PoPS"))
  )

# Helper to make ONE locus plot
make_one_locus_plot <- function(df_locus) {
  # order genes within this locus (top = highest combined_score)
  gene_levels <- df_locus %>%
    distinct(GENE, combined_score) %>%
    arrange(combined_score) %>%              # low->high
    pull(GENE)
  
  df_locus <- df_locus %>%
    mutate(GENE = factor(GENE, levels = gene_levels))
  
  ggplot(df_locus, aes(x = metric, y = GENE, fill = score)) +
    geom_tile(color = "black", linewidth = 0.25) +
    scale_fill_viridis_c(limits = c(0, 1), option = "C") +
    labs(title = paste0("Locus ", unique(df_locus$locus_label))) +
    theme_minimal(base_size = 9) +
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text.x = element_text(size = 7, face = "bold"),
      axis.text.y = element_text(size = 6),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 9),
      legend.position = "none"
    )
}

# Create list of 50 plots in locus order
plot_list <- plot_df %>%
  arrange(locus_label) %>%
  group_split(locus_label) %>%
  lapply(make_one_locus_plot)

# Stitch into one big figure
big_plot <- wrap_plots(plot_list, ncol = 10) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

# Add ONE shared legend (make it from a dummy plot)
legend_plot <- ggplot(plot_df, aes(metric, GENE, fill = score)) +
  geom_tile() +
  scale_fill_viridis_c(limits = c(0, 1), option = "C", name = "Score (0–1)") +
  theme_void() +
  theme(legend.position = "right")

big_plot <- big_plot + plot_layout(guides = "collect") &
  theme(legend.position = "right")

# Save
ggsave(
  filename = file.path(plot_dir, "gene_scoring_heatmaps_50panel_patchwork.png"),
  plot = big_plot,
  width = 28, height = 18, dpi = 300
)





