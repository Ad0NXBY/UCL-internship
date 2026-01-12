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
