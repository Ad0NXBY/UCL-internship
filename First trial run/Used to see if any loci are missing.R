lead_loci <- pops_lead_80k %>%
  transmute(locus_id = SNP, CHR = as.integer(CHR), BP = as.integer(BP)) %>%
  distinct()

top1_loci <- top1 %>%
  transmute(locus_id = locus_id) %>%
  distinct()

missing <- anti_join(lead_loci, top1_loci, by = "locus_id")
missing

missing_check <- missing %>%
  rowwise() %>%
  mutate(n_genes_in_window = nrow(
    pops_gene_annot %>%
      filter(CHR == CHR, abs(TSS - BP) <= window_bp)
  )) %>%
  ungroup()

missing_check

# genes near the missing locus
nearby_genes <- pops_gene_annot %>%
  filter(CHR == missing$CHR[1], abs(TSS - missing$BP[1]) <= window_bp) %>%
  select(ENSGID, GENE) %>%
  distinct()

# how many of these have PoPS scores?
nearby_genes %>%
  left_join(pops_preds %>% select(ENSGID, PoPS_score), by = "ENSGID") %>%
  summarize(
    genes_in_window = n(),
    genes_with_pops = sum(!is.na(PoPS_score))
  )
