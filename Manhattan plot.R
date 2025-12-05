library(data.table)
gwas <- fread("lvef_00000000.b38.gnorm.gz")


#Check the dimensions of the dataset
dim(gwas)

#Check the dataset to see what we are working with
head(gwas)

#Look at only the names of the dataset to see if they match with the manhattan plot 101 guide
names(gwas)

#Rename the columns to match the format required by qqman
setnames(gwas,
         old = c("chr_name","start_pos","var_id","pvalue"),
         new = c("CHR", "BP", "SNP", "P"))
head(gwas)

#Ensure that the columns are the right type
gwas[, CHR := as.numeric(CHR)]
gwas[, BP := as.numeric(BP)]
gwas[, P := as.numeric(P)]
gwas[, SNP := as.character(SNP)]

head(gwas)

#Remove rows with missing/invalid values
gwas <- gwas[!is.na(CHR) & !is.na(BP) & !is.na(P) & P > 0]

#install the qqman package and make the manhattan plot
install.packages("qqman")
library(qqman)
manhattan(gwas,
          chr = "CHR",
          bp = "BP",
          snp = "SNP",
          p = "P",
          logp = FALSE,
          genomewideline = 7.3,
          annotatePval = 7.3,
          annotateTop = TRUE)




library(tidyverse)
don <- gwas %>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(gwas, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot)

axisdf = don %>%
  group_by(CHR) %>%
  summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

don <- don %>%
  mutate(
    is_highlight = ifelse(P > 7.3, "yes", "no"),
    is_annotate  = ifelse(P > 7.3, "yes", "no")
  )

library(ggplot2)
library(ggrepel)

ggplot(don, aes(x = BPcum, y = P)) +
  geom_point(aes(color = as.factor(CHR)), alpha = 0.8, size = 1.3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
  scale_x_continuous(label = axisdf$CHR, breaks = axisdf$center) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_hline(yintercept = 7.3, linetype = "dashed", colour = "red") +  # genome-wide line
  geom_point(data=subset(don, is_highlight=="yes"), color = "orange", size = 2) +
  geom_label_repel(data=subset(don, is_annotate=="yes"), aes(label=SNP), size =2) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

#Extracting significant SNPs into a CSV file
sig_SNPs <- gwas %>% 
  filter(P > 7.3) %>%
  arrange(desc(P))

