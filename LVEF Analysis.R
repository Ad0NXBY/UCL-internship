#Directory for plots============================================================
plot_dir <- "plots"

# Create it if it doesn't exist
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}

data_dir <- "data"

# Create it if it doesn't exist
if (!dir.exists(data_dir)) {
  dir.create(data_dir, recursive = TRUE)
}

#Install packages===============================================================
install.packages('data.table')
install.packages('R.utils')
install.packages("ieugwasr")
install.packages("usethis")
install.packages("qqman")
install.packages("coloc")
install.packages("tidyverse")
install.packages("ggrepel")

library(data.table)
library(qqman)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(ieugwasr)
library(usethis)
library(coloc)


#80k CODE (MAKE SURE TO REPLACE NUMBER IF YOU WANT OTHER DATASET)===================
gwas_80k <- fread("80k_lvef_percent_000000.b38.gnorm.gz")

#Rename the columns to match the format required by qqman
setnames(gwas_80k,
         old = c("chr_name","start_pos","var_id","pvalue"),
         new = c("CHR", "BP", "SNP", "P"))

#Ensure that the columns are the right type
gwas_80k[, CHR := as.numeric(CHR)]
gwas_80k[, BP := as.numeric(BP)]
gwas_80k[, P := as.numeric(P)]
gwas_80k[, SNP := as.character(SNP)]

#Remove rows with missing/invalid values
gwas_80k <- gwas_80k[!is.na(CHR) & !is.na(BP) & !is.na(P) & P > 0]


#Send it towards the plot directory
png(
  filename = file.path(plot_dir, "80k updated manhattan plot.png"),
  width = 2000, height = 1200, res = 150
)
#Make the manhattan plot
manhattan(gwas_80k,
          chr = "CHR",
          bp = "BP",
          snp = "SNP",
          p = "P",
          logp = FALSE,
          genomewideline = 7.3,
          annotatePval = 7.3,
          annotateTop = TRUE)

dev.off()

don_80k <- gwas_80k %>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(gwas_80k, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot)

axisdf_80k = don_80k %>%
  group_by(CHR) %>%
  summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


don_80k <- don_80k %>%
  mutate(
    is_highlight = ifelse(P > 7.3, "yes", "no"),
    is_annotate  = ifelse(P > 7.3, "yes", "no")
  )



man_plot_80k <- ggplot(don_80k, aes(x = BPcum, y = P)) +
  geom_point(aes(color = as.factor(CHR)), alpha = 0.8, size = 1.3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
  scale_x_continuous(label = axisdf_80k$CHR, breaks = axisdf_80k$center) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_hline(yintercept = 7.3, linetype = "dashed", colour = "red") +  # genome-wide line
  geom_point(data=subset(don_80k, is_highlight=="yes"), color = "orange", size = 2) +
  geom_label_repel(data=subset(don_80k, is_annotate=="yes"), aes(label=SNP), size =2) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

ggsave(
  filename = file.path(plot_dir, "80k updated manhattan plot 2.png"),
  plot     = man_plot_80k,
  width    = 20,
  height   = 10,
  dpi      = 300
)



#Extracting significant SNPs into a CSV file
sig_SNPs_80k <- gwas_80k %>% 
  filter(P > 7.3) %>%
  arrange(desc(P))

#Save as TSV and CSV file
write.table(
  sig_SNPs_80k,
  file = file.path(data_dir, "80k significant_SNPs_gwas.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.csv(
  sig_SNPs_80k,
  file = file.path(data_dir, "80k significant_SNPs_gwas.csv"),
  row.names = FALSE
)

#Create QQPlot
png(
  filename = file.path(plot_dir, "80k qqplot.png"),
  width = 1400, height = 1000, res = 150
)
gwas_80k[, P_raw := 10^(-P)]
qqman::qq(gwas_80k$P_raw)

dev.off()

#perform LD clumping

#usethis::edit_r_environ() #use this to open to add the token
#OPENGWAS_JWT=eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiJibHlsMWcyMUBzb3Rvbi5hYy51ayIsImlhdCI6MTc2NTM4MjM5NywiZXhwIjoxNzY2NTkxOTk3fQ.QI90OnB-H3hyMnCcSoNH8j9sf0BHf1HYfmYG2RTZXQFCYU6adk3GddnMDHjM_5AduCNtXVGJ4ld_dYPbivw7K4Tw5KBERHHxXYJO8Ah9g_0jCJbxoXpx8LcySt7virQ9d2WJeOXr8ARCC6r-yYvc3TQRy5GN72Jl7Fnte7XKYdqu8S87WGEkjpvWuwoY4MiVUBPgj6ad8FQ3pAP4MadvqDtV_deb75sjo5pcUInPX7o82Gj_xc_Gg81oYb_hroaamP_0qbejpPx_0ZHOjLBQdCYCE5nOMj5m2TEvfF6MOef6wRo1pzkQKbRNkD-qqPG6GDoVybLIY8vaoJM3BGHSHw

ieugwasr::get_opengwas_jwt() #this is to double check things are working
ieugwasr::user() 

gw_thresh_log10 <- 7.3

clump_dat_80k <- gwas_80k %>%
  dplyr::filter(P > gw_thresh_log10) %>%   # P = -log10(p)
  dplyr::transmute(
    rsid = SNP,
    pval = 10^(-P)
  )

clumped_80k <- ieugwasr::ld_clump(
  dat      = clump_dat_80k,
  clump_kb = 10000,
  clump_r2 = 0.001,
  clump_p  = 5e-8,
  pop      = "EUR"
)

lead_SNPs_ld_80k <- gwas_80k %>%
  dplyr::inner_join(clumped_80k, by = c("SNP" = "rsid"))

# Save LD-clumped lead SNPs
write.table(
  lead_SNPs_ld_80k,
  file = file.path(data_dir, "80k_lead_SNPs_LD_clumped.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.csv(
  lead_SNPs_ld_80k,
  file = file.path(data_dir, "80k_lead_SNPs_LD_clumped.csv"),
  row.names = FALSE
)





#40k CODE (MAKE SURE TO REPLACE NUMBER IF YOU WANT OTHER DATASET)===================
gwas_40k <- fread("40k_lvef_00000000.b38.gnorm.gz")

#Rename the columns to match the format required by qqman
setnames(gwas_40k,
         old = c("chr_name","start_pos","var_id","pvalue"),
         new = c("CHR", "BP", "SNP", "P"))

#Ensure that the columns are the right type
gwas_40k[, CHR := as.numeric(CHR)]
gwas_40k[, BP := as.numeric(BP)]
gwas_40k[, P := as.numeric(P)]
gwas_40k[, SNP := as.character(SNP)]

#Remove rows with missing/invalid values
gwas_40k <- gwas_40k[!is.na(CHR) & !is.na(BP) & !is.na(P) & P > 0]


#Send it towards the plot directory
png(
  filename = file.path(plot_dir, "40k updated manhattan plot.png"),
  width = 2000, height = 1200, res = 150
)
#Make the manhattan plot
manhattan(gwas_40k,
          chr = "CHR",
          bp = "BP",
          snp = "SNP",
          p = "P",
          logp = FALSE,
          genomewideline = 7.3,
          annotatePval = 7.3,
          annotateTop = TRUE)

dev.off()


don_40k <- gwas_40k %>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(gwas_40k, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot)

axisdf_40k = don_40k %>%
  group_by(CHR) %>%
  summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

don_40k <- don_40k %>%
  mutate(
    is_highlight = ifelse(P > 7.3, "yes", "no"),
    is_annotate  = ifelse(P > 7.3, "yes", "no")
  )



man_plot_40k <- ggplot(don_40k, aes(x = BPcum, y = P)) +
  geom_point(aes(color = as.factor(CHR)), alpha = 0.8, size = 1.3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
  scale_x_continuous(label = axisdf_40k$CHR, breaks = axisdf_40k$center) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_hline(yintercept = 7.3, linetype = "dashed", colour = "red") +  # genome-wide line
  geom_point(data=subset(don_40k, is_highlight=="yes"), color = "orange", size = 2) +
  geom_label_repel(data=subset(don_40k, is_annotate=="yes"), aes(label=SNP), size =2) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

ggsave(
  filename = file.path(plot_dir, "40k updated manhattan plot 2.png"),
  plot     = man_plot_40k,
  width    = 20,
  height   = 10,
  dpi      = 300
)


#Extracting significant SNPs into a CSV file
sig_SNPs_40k <- gwas_40k %>% 
  filter(P > 7.3) %>%
  arrange(desc(P))

#Save as TSV and CSV file
write.table(
  sig_SNPs_40k,
  file = file.path(data_dir, "40k significant_SNPs_gwas.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.csv(
  sig_SNPs_40k,
  file = file.path(data_dir, "40k significant_SNPs_gwas.csv"),
  row.names = FALSE
)

#Create QQPlot
png(
  filename = file.path(plot_dir, "40k qqplot.png"),
  width = 1400, height = 1000, res = 150
)
gwas_40k[, P_raw := 10^(-P)]
qqman::qq(gwas_40k$P_raw)

dev.off()

#perform LD clumping

#usethis::edit_r_environ() #use this to open to add the token
#OPENGWAS_JWT=eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiJibHlsMWcyMUBzb3Rvbi5hYy51ayIsImlhdCI6MTc2ODI4ODA5MywiZXhwIjoxNzY5NDk3NjkzfQ.LEsghLU4QbsclW7oGtNG6j8tGq-5b781hGUdKblZygTim22GhHsP7AvhpJ90oBn87hyH7oFZKhbItuJYdV-oRy__J7OuKtP9ZOqhOUnd5zUXgGjGyJBZn3du7B_E_OosSsSQTw1Pa87LHbdGifYX627Eh8_ucs1zDs8pcD2E9QfzbN1_UdmicEcfh-4AXhvRZcie9heFezINYZ42WgO3D4MRmVGjAR8fyf8a3zx68xbUmCwsU0X-kwnzBRcaaGceC1kwBscL9nRaim6QqaFlOHAaQXD_3B5aLK56fxjTTrD0UBnKEQq1tmT3UPgOwP7P3E5AJat3I1I1GupJ1bSSLg
#This will change i think

ieugwasr::get_opengwas_jwt() #this is to double check things are working
ieugwasr::user() 

gw_thresh_log10 <- 7.3

clump_dat_40k <- gwas_40k %>%
  dplyr::filter(P > gw_thresh_log10) %>%   # P = -log10(p)
  dplyr::transmute(
    rsid = SNP,
    pval = 10^(-P)
  )

clumped_40k <- ieugwasr::ld_clump(
  dat      = clump_dat_40k,
  clump_kb = 10000,
  clump_r2 = 0.001,
  clump_p  = 5e-8,
  pop      = "EUR"
)

lead_SNPs_ld_40k <- gwas_40k %>%
  dplyr::inner_join(clumped_40k, by = c("SNP" = "rsid"))

# Save LD-clumped lead SNPs
write.table(
  lead_SNPs_ld_40k,
  file = file.path(data_dir, "40k_lead_SNPs_LD_clumped.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.csv(
  lead_SNPs_ld_40k,
  file = file.path(data_dir, "40k_lead_SNPs_LD_clumped.csv"),
  row.names = FALSE
)




#Colocolization Analysis========================================================
#define the genomic window
window_bp <- 500000

#extract region from the GWAS table)
get_region_dt <- function(gwas, chr, bp, window_bp) {
  # subset region using data.table (fast + low memory)
  r <- gwas[CHR == chr & BP >= (bp - window_bp) & BP <= (bp + window_bp),
            .(SNP, BP, P, effect_allele, other_allele,
              effect_size, standard_error, effect_allele_freq, number_of_samples)]
  
  # deduplicate: keep the most significant row per SNP (largest -log10(p) = P)
  if (nrow(r) == 0) return(r)
  r <- r[order(-P)]
  r <- r[, .SD[1], by = SNP]
  r
}

#match alleles from 40k to 80k
matching <- function(m) {
  same <- m$effect_allele_80 == m$effect_allele_40 & m$other_allele_80 == m$other_allele_40
  swap <- m$effect_allele_80 == m$other_allele_40 & m$other_allele_80 == m$effect_allele_40
  m <- m[same | swap, , drop = FALSE]
  
  swapped_rows <- which(swap[same | swap])
  if (length(swapped_rows) > 0) {
    m$effect_size_40[swapped_rows] <- -m$effect_size_40[swapped_rows]
    m$effect_allele_freq_40[swapped_rows] <- 1 - m$effect_allele_freq_40[swapped_rows]
  }
  
  m
}

#build coloc dataset
make_dataset <- function(df, beta, se, eaf, n, pos) {
  maf <- pmin(df[[eaf]], 1 - df[[eaf]])
  list(
    snp      = df$SNP,
    position = df[[pos]],
    beta     = df[[beta]],
    varbeta  = (df[[se]])^2,
    MAF      = maf,
    N        = as.integer(round(median(df[[n]], na.rm = TRUE))),
    type     = "quant"
  )
}

#run colocolization for each 80k lead SNP hit
coloc_results <- lapply(seq_len(nrow(lead_SNPs_ld_80k)), function(i) {
  
  lead <-lead_SNPs_ld_80k$SNP[i]
  chr <- lead_SNPs_ld_80k$CHR[i]
  bp <- lead_SNPs_ld_80k$BP[i]
  
  ##extract region per GWAS
  r80 <- get_region_dt(gwas_80k,chr, bp, window_bp)
  r40 <- get_region_dt(gwas_40k, chr, bp, window_bp)
  
  ##keep SNPs present in both
  m <- merge(r80, r40, by = "SNP", suffixes = c("_80", "_40"))
  
  if (nrow(m) < 50)  {
    return(data.frame(lead_snp = lead, CHR = chr, BP = bp, nsnps = nrow(m),
                      PP.H0 = NA, PP.H1 = NA, PP.H2 = NA, PP.H3 = NA, PP.H4 = NA))
  }
  
  ##match alleles (flip 40k betas if needed)
  m <- matching(m)
  
  if (nrow(m) < 50) {
    return(data.frame(lead_snp = lead, CHR = chr, BP = bp, nsnps = nrow(m),
                      PP.H0 = NA, PP.H1 = NA, PP.H2 = NA, PP.H3 = NA, PP.H4 = NA))
  }
  
  ##make colocolization input dataset
  d1 <- make_dataset(m,
                     beta = "effect_size_80",
                     se   = "standard_error_80",
                     eaf  = "effect_allele_freq_80",
                     n    = "number_of_samples_80",
                     pos  = "BP_80"
  )
  
  d2 <- make_dataset(m,
                     beta = "effect_size_40",
                     se   = "standard_error_40",
                     eaf  = "effect_allele_freq_40",
                     n    = "number_of_samples_40",
                     pos  = "BP_40"
  )
  
  ##run coloc
  res <- coloc.abf(d1, d2)
  s <- res$summary
  
  data.frame(
    lead_snp = lead, CHR = chr, BP = bp, nsnps = nrow(m),
    PP.H0 = unname(s["PP.H0.abf"]),
    PP.H1 = unname(s["PP.H1.abf"]),
    PP.H2 = unname(s["PP.H2.abf"]),
    PP.H3 = unname(s["PP.H3.abf"]),
    PP.H4 = unname(s["PP.H4.abf"])
  )
})

coloc_results_df <- dplyr::bind_rows(coloc_results)

write.csv(
  coloc_results_df,
  file = file.path(data_dir, "coloc_80k_vs_40k_per_80k_locus.csv"),
  row.names = FALSE
)

#Visualizing the plots
plot_df <- coloc_results_df %>%
  mutate(
    locus = paste0("chr", CHR, ":", BP),
    PP.H4 = as.numeric(PP.H4)
  ) %>%
  arrange(PP.H4)

pph4_plot <-ggplot(plot_df, aes(x = reorder(locus, PP.H4), y = PP.H4)) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
  coord_flip() +
  labs(
    x = "80k locus (lead SNP window)",
    y = "PP.H4 (shared causal variant)",
    title = "Colocalisation posterior per 80k locus (80k vs 40k)"
  ) +
  theme_bw()

ggsave(
  filename = file.path(plot_dir, "coloc_PP.H4_per_80k_locus_posterior_probability plot.png"),
  plot = pph4_plot,
  width = 10,
  height = 12,
  dpi = 300
)



#Working on public dataset======================================================
gwas_public <- fread("MRI_lvef_filtered/MRI_lvef_filtered.tsv")
#check the output of the dataset to change the names and make sure they match with what we want
names(gwas_public)
head(gwas_public)

#Change the names of the headers to match what we want
setnames(gwas_public,
         old = c("SNP", "CHR", "BP", "P_LINREG", "ALLELE1", "ALLELE0", "A1FREQ", "BETA", "SE", "N"),
         new = c("SNP", "CHR", "BP", "P_raw", "effect_allele", "other_allele", "effect_allele_freq", "effect_size", "standard_error", "number_of_samples"))

#Ensure that the columns are the right type
gwas_public[, CHR := as.numeric(CHR)]
gwas_public[, BP  := as.numeric(BP)]
gwas_public[, SNP := as.character(SNP)]
gwas_public[, P_raw := as.numeric(P_raw)]
gwas_public[, effect_size := as.numeric(effect_size)]
gwas_public[, standard_error := as.numeric(standard_error)]
gwas_public[, effect_allele_freq := as.numeric(effect_allele_freq)]
gwas_public[, number_of_samples := as.integer(number_of_samples)]

#convert raw p-value to match with column "P"
gwas_public[, P := -log10(P_raw)]

#Remove missing/invalid rows
gwas_public <- gwas_public[
  !is.na(CHR) & !is.na(BP) & !is.na(SNP) &
    !is.na(P_raw) & P_raw > 0
]


#LD clumping
clump_dat_public <- gwas_public %>%
  filter(P > gw_thresh_log10) %>%
  transmute(
    rsid = SNP,
    pval = P_raw
  )

clumped_pub <- ieugwasr::ld_clump(
  dat      = clump_dat_public,
  clump_kb = 10000,
  clump_r2 = 0.001,
  clump_p  = 5e-8,
  pop      = "EUR"
)

lead_SNPs_ld_pub <- gwas_public %>%
  inner_join(clumped_pub, by = c("SNP" = "rsid"))

# Save LD-clumped lead SNPs
write.table(
  lead_SNPs_ld_pub,
  file = file.path(data_dir, "published_LVEF_lead_SNPs_LD_clumped.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.csv(
  lead_SNPs_ld_pub,
  file = file.path(data_dir, "published_LVEF_lead_SNPs_LD_clumped.csv"),
  row.names = FALSE
)

#Establishing true novelty
novelty_table <- lead_SNPs_ld_80k %>%
  transmute(SNP, CHR, BP, P_80k = P) %>%
  rowwise() %>%
  mutate(
    overlaps_40k = any(lead_SNPs_ld_40k$CHR == CHR & abs(lead_SNPs_ld_40k$BP - BP) <= window_bp),
    overlaps_published = any(lead_SNPs_ld_pub$CHR == CHR & abs(lead_SNPs_ld_pub$BP - BP) <= window_bp),
    novelty = case_when(
      overlaps_40k | overlaps_published ~ "Not novel (overlaps known)",
      TRUE ~ "Novel (no overlaps)"
    )
  ) %>%
  ungroup()

#Creating data file for the novel SNPs that has all the known information we have
novel_loci_for_80k <- novelty_table %>%
  filter(novelty == "Novel (no overlaps)")

novel_loci_80k_full <- novel_loci_for_80k %>%
  left_join(lead_SNPs_ld_80k, by = c("SNP", "CHR", "BP"))


write.csv(
  novel_loci_80k_full,
  file = file.path(data_dir, "80k_novel_LVEF_loci.csv"),
  row.names = FALSE
)

#Create Manhattan Plot with novel
novel_snps <- novel_loci_for_80k$SNP

novel_don_80k <- don_80k %>%
  mutate(
    is_novel = SNP %in% novel_snps
  )

novel_man_80k_plot <- ggplot(novel_don_80k, aes(x = BPcum, y = P)) +
  geom_point(aes(color = as.factor(CHR)), alpha = 0.8, size = 1.3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22)) +
  scale_x_continuous(label = axisdf_80k$CHR, breaks = axisdf_80k$center) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_hline(yintercept = 7.3, linetype = "dashed", colour = "red") +
  geom_point(data = subset(novel_don_80k, is_novel), color = "orange", size = 2.5) +
  geom_label_repel(
    data = subset(novel_don_80k, is_novel),
    aes(label = SNP),
    size = 2,
    max.overlaps = 50
  ) +
  labs(
    x = "Chromosomal position",
    y = expression(-log[10](italic(p))),
    title = "80k LVEF GWAS with true novel loci highlighted"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

ggsave(
  filename = file.path(plot_dir, "80k_manhattan_true_novel_labeled.png"),
  plot     = novel_man_80k_plot,
  width    = 20,
  height   = 10,
  dpi      = 300
)

#Gene Prioritization============================================================
#Create the file for MAGMA
magma_80k <- gwas_80k[, .(
  SNP,
  P = P_raw,
  N = number_of_samples
)]

fwrite(
  magma_80k,
  file = "80k_LVEF_magma.txt",
  sep = "\t"
)
#Create PoPs compatible gene annotation file
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt")
library(biomaRt)
library(data.table)

mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

gene_annot <- getBM(
  attributes = c(
    "ensembl_gene_id",
    "external_gene_name",
    "chromosome_name",
    "transcription_start_site"
  ),
  mart = mart
)

gene_annot <- as.data.table(gene_annot)
setnames(gene_annot,
         old = c("ensembl_gene_id", "external_gene_name", "chromosome_name", "transcription_start_site"),
         new = c("ENSGID", "GENE", "CHR", "TSS"))

# Keep autosomes + X, drop weird contigs
gene_annot <- gene_annot[CHR %in% c(as.character(1:22), "X")]

# Save to a file PoPS expects (tab-delimited)
fwrite(gene_annot,
       file = "pops_gene_annot_hg38.tsv",
       sep = "\t")


#Create an Entrez to ENSG mapping
library(biomaRt)
library(data.table)

mart <- useEnsembl(biomart="genes", dataset="hsapiens_gene_ensembl")

map <- getBM(
  attributes = c("entrezgene_id", "ensembl_gene_id"),
  mart = mart
)

map <- as.data.table(map)
setnames(map, c("entrezgene_id","ensembl_gene_id"), c("GENE","ENSGID"))
map <- map[!is.na(GENE) & ENSGID != ""]
map[, GENE := as.character(GENE)]

fwrite(map, file="entrez_to_ensg.tsv", sep="\t")

#PoPs gene prioritization
top_n <- 10 #This can be set to 5

#Load inputs 
pops_lead_80k <- fread(file.path(data_dir, "80k_lead_SNPs_LD_clumped.csv"))
pops_preds <- fread(file.path(data_dir, "Nikhil_lvefb37_pops.preds"))
pops_gene_annot <- fread(file.path(data_dir, "pops_gene_annot_hg38.tsv"))

#Making sure ENSG and score column are detected
ensg_col <- names(pops_preds)[grepl("ENSG", names(pops_preds), ignore.case = TRUE)][1]
if (is.na(ensg_col)) ensg_col <- names(pops_preds)[1]

score_col <- names(pops_preds)[grepl("pred|score", names(pops_preds), ignore.case = TRUE)][1]
if (is.na(score_col)) score_col <- names(pops_preds)[2]

pops_preds <- pops_preds %>%
  rename(ENSGID = all_of(ensg_col),
         PoPS_score = all_of(score_col)) %>%
  mutate(
    ENSGID = as.character(ENSGID),
    PoPS_score = as.numeric(PoPS_score)
  )

#Gene Annotation
pops_gene_annot <- pops_gene_annot %>%
  mutate(
    ENSGID = as.character(ENSGID),
    CHR = as.character(CHR),
    CHR = gsub("^chr", "", CHR),
    CHR = ifelse(CHR %in% as.character(1:22), CHR, NA),
    CHR = as.integer(CHR),
    TSS = as.integer(TSS),
    GENE = ifelse(is.na(GENE) | GENE == "", ENSGID, GENE)
  ) %>%
  filter(!is.na(CHR), !is.na(TSS))

#Join PoPs scores to gene positions
pops_annot <- pops_preds %>%
  inner_join(pops_gene_annot, by = "ENSGID") %>%
  filter(!is.na(PoPS_score))

#Map genes to each locus and keep top N
top_genes_per_locus <- pops_lead_80k %>%
  transmute(
    locus_id  = SNP,
    locus     = paste0("chr", CHR, ":", BP),
    locus_CHR = as.integer(CHR),
    locus_BP  = as.integer(BP)
  ) %>%
  distinct() %>%
  rowwise() %>%
  mutate(
    genes = list(
      pops_annot %>%
        filter(CHR == locus_CHR, abs(TSS - locus_BP) <= window_bp) %>%
        arrange(desc(PoPS_score)) %>%
        slice_head(n = top_n)
    )
  ) %>%
  ungroup() %>%
  tidyr::unnest(genes) %>%
  select(
    locus_id, locus, locus_CHR, locus_BP,
    ENSGID, GENE, CHR, TSS, PoPS_score
  ) %>%
  arrange(locus, desc(PoPS_score))


write.csv(
  top_genes_per_locus,
  file = file.path(data_dir, "80k_PoPS_top_genes_per_locus.csv"),
  row.names = FALSE
)

# Save top1 per locus
top1 <- top_genes_per_locus %>%
  group_by(locus) %>%
  slice_max(PoPS_score, n = 1, with_ties = FALSE) %>%
  ungroup()

write.csv(
  top1,
  file = file.path(data_dir, "80k_PoPS_top1_gene_per_locus.csv"),
  row.names = FALSE
)

##Gene Prioritization Heatmap (Change the top_n here if you want to get more/less genes)
library(dplyr)
library(ggplot2)
library(ggtext)

window_bp <- 500000
top_n <- 10

get_top_genes_one_locus <- function(locus_chr, locus_bp, locus_label, top_n = 10, window_bp = 500000) {
  
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
  filename = file.path(plot_dir, "PoPS_per_locus_top_genes_2.png"),
  plot = p_big,
  width = 40,   # increase width
  height = 20,   # increase height
  dpi = 300
)

ggsave(
  filename = file.path(plot_dir, "PoPS_per_locus_top_genes_2.pdf"),
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


#Finding nearest gene to each lead loci=========================================
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
