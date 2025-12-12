library(data.table)
gwas <- fread("lvef_00000000.b38.gnorm.gz")

# Rename columns
setnames(gwas,
         old = c("chr_name","start_pos","var_id","pvalue"),
         new = c("CHR", "BP", "SNP", "LOG10P"))  # LOG10P = –log10(p)

# Types
gwas[, CHR    := as.numeric(CHR)]
gwas[, BP     := as.numeric(BP)]
gwas[, LOG10P := as.numeric(LOG10P)]
gwas[, SNP    := as.character(SNP)]

# Remove dodgy rows
gwas <- gwas[!is.na(CHR) & !is.na(BP) & !is.na(LOG10P)]

# If you need *raw* p-values at any point:
gwas[, P := 10^(-LOG10P)]   # this is the actual p

library(dplyr)
library(ggplot2)
library(ggrepel)

gw_thresh_log10 <- 7.3

don <- gwas %>%
  group_by(CHR) %>%
  summarise(chr_len = max(BP), .groups = "drop") %>%
  mutate(tot = cumsum(chr_len) - chr_len) %>%
  select(-chr_len) %>%
  left_join(gwas, by = "CHR") %>%
  arrange(CHR, BP) %>%
  mutate(BPcum = BP + tot)

axisdf <- don %>%
  group_by(CHR) %>%
  summarise(center = (max(BPcum) + min(BPcum)) / 2)

don <- don %>%
  mutate(
    is_highlight = LOG10P >= gw_thresh_log10,
    is_annotate  = LOG10P >= gw_thresh_log10
  )

ggplot(don, aes(x = BPcum, y = LOG10P)) +
  geom_point(aes(color = as.factor(CHR)), alpha = 0.8, size = 1.3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
  scale_x_continuous(label = axisdf$CHR, breaks = axisdf$center) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_hline(yintercept = gw_thresh_log10, linetype = "dashed", colour = "red") +
  geom_point(data = subset(don, is_highlight), color = "orange", size = 2) +
  geom_label_repel(
    data  = subset(don, is_annotate),
    aes(label = SNP),
    size  = 2
  ) +
  theme_bw() +
  theme(
    legend.position        = "none",
    panel.border           = element_blank(),
    panel.grid.major.x     = element_blank(),
    panel.grid.minor.x     = element_blank()
    
  )


library(qqman)

# ensure P exists
gwas[, P := 10^(-LOG10P)]

qq(gwas$P)


gw_thresh_log10 <- 7.3

library(dplyr)

sig_SNPs <- gwas %>%
  filter(LOG10P >= gw_thresh_log10) %>%
  arrange(desc(LOG10P))   # biggest –log10(p) = smallest p

write.table(
  sig_SNPs,
  file = "significant_SNPs_gwas.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.csv(
  sig_SNPs,
  file = "significant_SNPs_gwas.csv",
  row.names = FALSE
)

library(data.table)

sig <- as.data.table(sig_SNPs)  # from above
setorder(sig, -LOG10P)          # most significant first

window_kb <- 500
window_bp <- window_kb * 1000

lead_hits <- data.table()

while (nrow(sig) > 0) {
  # take the top SNP
  lead <- sig[1]
  lead_hits <- rbind(lead_hits, lead)
  
  # remove SNPs within ±window_bp of this SNP on same chromosome
  sig <- sig[!(
    CHR == lead$CHR &
      BP  >= (lead$BP - window_bp) &
      BP  <= (lead$BP + window_bp)
  )]
}

lead_SNPs <- lead_hits

# Save
fwrite(lead_SNPs, "lead_SNPs_approx_independent.tsv", sep = "\t")
fwrite(lead_SNPs, "lead_SNPs_approx_independent.csv")

install.packages("biomaRt")
library(biomaRt)
library(dplyr)

mart <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")

# Get gene coordinates
genes <- getBM(
  attributes = c("ensembl_gene_id",
                 "hgnc_symbol",
                 "chromosome_name",
                 "start_position",
                 "end_position",
                 "gene_biotype"),
  mart = mart
)

# Keep autosomes & X/Y/MT, match your CHR format
genes <- genes %>%
  filter(chromosome_name %in% as.character(unique(lead_SNPs$CHR))) %>%
  rename(
    CHR = chromosome_name,
    gene_start = start_position,
    gene_end   = end_position
  )

lead_dt <- as.data.table(lead_SNPs)

# Find *nearest gene* to each lead SNP
nearest_list <- lapply(1:nrow(lead_dt), function(i) {
  snp <- lead_dt[i]
  cand <- genes %>% filter(CHR == as.character(snp$CHR))
  
  cand <- cand %>%
    mutate(
      distance_bp = pmin(
        abs(gene_start - snp$BP),
        abs(gene_end   - snp$BP)
      )
    ) %>%
    arrange(distance_bp) %>%
    slice(1)
  
  cbind(as.data.frame(snp), cand)
})

lead_annotated <- bind_rows(nearest_list)

write.csv(lead_annotated, "lead_SNPs_nearest_gene.csv", row.names = FALSE)

