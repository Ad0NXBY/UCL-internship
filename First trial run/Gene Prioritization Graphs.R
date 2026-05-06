library(dplyr)
library(ggplot2)
library(ggtext)

window_bp <- 500000
top_n <- 4

get_top_genes_one_locus <- function(locus_chr, locus_bp, locus_label, top_n = 4, window_bp = 500000) {
  
  pops_gene_annot %>%
    mutate(CHR = as.integer(CHR),
           TSS = as.integer(TSS)) %>%
    inner_join(pops_preds %>% select(ENSGID, PoPS_score), by = "ENSGID") %>%
    filter(CHR == locus_chr, abs(TSS - locus_bp) <= window_bp) %>%
    mutate(dist = abs(TSS - locus_bp)) %>%
    arrange(desc(PoPS_score), dist) %>%
    group_by(ENSGID) %>%
    slice(1) %>%
    ungroup() %>%
    arrange(desc(PoPS_score)) %>%
    slice_head(n = top_n) %>%
    mutate(locus_label = .env$locus_label)   
}

lead_loci <- pops_lead_80k %>%
  transmute(
    locus_idx = row_number(),
    locus_label = paste0("Locus ", locus_idx),
    locus_chr = as.integer(CHR),
    locus_bp  = as.integer(BP),
    locus_snp = SNP
  ) %>%
  distinct()

top_genes_per_locus <- purrr::pmap_dfr(
  lead_loci,
  function(locus_idx, locus_label, locus_chr, locus_bp, locus_snp) {
    get_top_genes_one_locus(locus_chr, locus_bp, locus_label, top_n = top_n, window_bp = window_bp)
  }
)

# scale PoPS
rng <- range(top_genes_per_locus$PoPS_score, na.rm = TRUE)
top_genes_per_locus <- top_genes_per_locus %>%
  mutate(PoPS_scaled_01 = if (diff(rng) == 0) 0 else (PoPS_score - rng[1]) / diff(rng))

# rank within locus
plot_df <- top_genes_per_locus %>%
  group_by(locus_label) %>%
  arrange(desc(PoPS_scaled_01), .by_group = TRUE) %>%
  mutate(gene_rank = row_number()) %>%
  ungroup() %>%
  mutate(locus_number = readr::parse_number(locus_label)) %>%
  arrange(locus_number) %>%
  mutate(locus_label = factor(locus_label, levels = unique(locus_label)))

novel_snps <- novel_loci_for_80k$SNP

lead_loci_flags <- lead_loci %>%  
  transmute(
    locus_label,
    is_novel = locus_snp %in% novel_snps
  )

plot_df <- plot_df %>%
  left_join(lead_loci_flags, by = "locus_label") %>%
  mutate(
    is_novel = tidyr::replace_na(is_novel, FALSE),
    locus_label_md = ifelse(is_novel,
                            paste0("<b>", as.character(locus_label), "</b>"),
                            as.character(locus_label)),
    locus_label_md = factor(locus_label_md,
                            levels = unique(locus_label_md[order(locus_number)]))
  )
#Create Heatmap

p_big <- ggplot(plot_df, aes(x = gene_rank, y = 1, fill = PoPS_scaled_01)) +
  geom_tile(color = "black", linewidth = 0.25, width = 0.95, height = 0.95) +
  geom_text(aes(label = GENE), y = 1.6, size = 5, angle = 45, hjust = 0) +
  scale_fill_viridis_c(limits = c(0, 1), option = "C") +
  facet_wrap(~ locus_label_md, ncol = 10) +
  scale_y_continuous(limits = c(0.5, 2), breaks = NULL) +
  labs(
    x = NULL, y = NULL, fill = "PoPS (0–1)",
    title = "PoPS gene prioritisation per locus (top 3–4 genes)",
    subtitle = "Bolded loci indicate novel associations"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    strip.text = ggtext::element_markdown(face = "plain", size = 14),
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 18,
      margin = margin(b = 6)
    ),
    plot.subtitle = element_text(
      hjust = 0.5,
      size = 13,
      margin = margin(b = 10)
    )
  )


p_big

ggsave(
  filename = file.path(plot_dir, "PoPS_per_locus_top_genes.png"),
  plot = p_big,
  width = 40,   # increase width
  height = 20,   # increase height
  dpi = 300
)

ggsave(
  filename = file.path(plot_dir, "PoPS_per_locus_top_genes.pdf"),
  plot = p_big,
  width = 40,   # increase width
  height = 20   # increase height
)

#Loci table (Chr, RSID, effect allele, other allele, beta, p-value, prioritized)
top1_gene <- top_genes_per_locus %>%
  group_by(locus_label) %>%
  arrange(desc(PoPS_scaled_01), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  select(locus_label, prioritised_gene = GENE, PoPS_score, PoPS_scaled_01)

locus_table_80k <- pops_lead_80k %>%
  transmute(
    locus_label = paste0("Locus ", row_number()),
    chr = CHR,
    rsid = SNP,
    effect_allele = effect_allele,
    other_allele  = other_allele,
    beta = effect_size,
    p_value = P_raw
  ) %>%
  left_join(top1_gene, by = "locus_label") %>%
  select(chr, rsid, effect_allele, other_allele, beta, p_value, prioritised_gene)

# Save
write.csv(
  locus_table_80k,
  file = file.path(data_dir, "80k_loci_summary_with_PoPS_gene.csv"),
  row.names = FALSE
)
