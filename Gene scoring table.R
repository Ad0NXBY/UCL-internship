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
