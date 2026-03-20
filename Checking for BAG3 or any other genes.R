# check BAG3 annotation row(s)
bag3_annot <- pops_gene_annot %>% filter(GENE == "BAG3")
bag3_annot

# check distance to the lead SNP (use the top SNP you showed)
lead_bp <- 121415685
bag3_annot %>%
  mutate(dist = abs(TSS - lead_bp)) %>%
  arrange(dist) %>%
  head(10)

# BAG3 in annotation?
pops_gene_annot %>% filter(GENE == "BAG3") %>% head()

# its ENSGID
bag3_ids <- pops_gene_annot %>% filter(GENE == "BAG3") %>% pull(ENSGID) %>% unique()
bag3_ids

# does PoPS preds contain that ENSGID?
pops_preds_public %>% filter(ENSGID %in% bag3_ids) %>% head()

top_genes_per_locus %>% filter(GENE == "BAG3")


pops_preds_public %>%
  filter(ENSGID == "ENSG00000151929") %>%  # BAG3
  select(ENSGID, PoPS_score = PoPS_score, Y)

# rank BAG3 among all genes:
pops_preds_public %>%
  mutate(rank = rank(-PoPS_score, ties.method = "min")) %>%
  filter(ENSGID == "ENSG00000151929") %>%
  select(ENSGID, PoPS_score, rank) 

pops_preds_public %>%
  mutate(rank = rank(-PoPS_score, ties.method="min")) %>%
  filter(ENSGID %in% c("ENSG00000151929","<TTN_ENSGID_HERE>")) %>%
  select(ENSGID, PoPS_score, Y, rank)

pops_lead_public %>% arrange(P_raw) %>% head(10) %>% select(SNP, CHR, BP, P_raw)

pops_gene_annot %>% filter(GENE == "TTN") %>% distinct(ENSGID)

pops_preds_public %>%
  mutate(rank = rank(-PoPS_score, ties.method="min")) %>%
  filter(ENSGID == "<TTN_ENSGID>") %>%
  select(ENSGID, PoPS_score, Y, rank)

pops_preds_public %>%
  mutate(rank = rank(-PoPS_score, ties.method="min")) %>%
  filter(ENSGID == "ENSG00000155657") %>%
  select(ENSGID, PoPS_score, Y, rank)
