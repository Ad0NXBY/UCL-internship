library(dplyr)

gene_annot <- pops_gene_annot %>%
  mutate(
    CHR = as.integer(CHR),
    TSS = as.integer(TSS),
    GENE = ifelse(is.na(GENE) | GENE == "", ENSGID, GENE)
  ) %>%
  filter(!is.na(CHR), !is.na(TSS))

lead_loci <- pops_lead_80k %>%
  transmute(
    locus_label = paste0("Locus ", row_number()),
    rsid = SNP,
    CHR = as.integer(CHR),
    BP  = as.integer(BP)
  )

#Compute nearest gene (TSS distance)
nearest_tbl <- lead_loci %>%
  rowwise() %>%
  mutate(
    nearest_gene = {
      chr_i <- CHR
      bp_i  <- BP
      g <- gene_annot %>% filter(CHR == chr_i)
      g$GENE[ which.min(abs(g$TSS - bp_i)) ]
    },
    nearest_dist_bp = {
      chr_i <- CHR
      bp_i  <- BP
      g <- gene_annot %>% filter(CHR == chr_i)
      min(abs(g$TSS - bp_i))
    }
  ) %>%
  ungroup()

#Compare nearest gene vs top PoPs gene
comparison <- nearest_tbl %>%
  left_join(top1_gene %>% select(locus_label, pops_top_gene = prioritised_gene),
            by = "locus_label") %>%
  mutate(same_gene = nearest_gene == pops_top_gene)

write.csv(
  comparison,
  file = file.path(data_dir, "nearest_gene_vs_PoPS_top_gene.csv"),
  row.names = FALSE
)





#Add nearest gene into 50 loci summary table
locus_table_80k <- pops_lead_80k %>%
  transmute(
    locus_label = paste0("Locus ", row_number()),
    chr = as.integer(CHR),
    rsid = SNP,
    effect_allele = effect_allele,
    other_allele  = other_allele,
    beta = effect_size,
    p_value = P_raw
  ) %>%
  left_join(top1_gene %>% select(locus_label, prioritised_gene), by = "locus_label") %>%
  left_join(nearest_tbl %>% select(locus_label, nearest_gene), by = "locus_label") %>%
  select(chr, rsid, effect_allele, other_allele, beta, p_value, prioritised_gene, nearest_gene)

write.csv(
  locus_table_80k,
  file = file.path(data_dir, "80k_loci_summary_table_with_PoPS_and_nearest_gene.csv"),
  row.names = FALSE
)
