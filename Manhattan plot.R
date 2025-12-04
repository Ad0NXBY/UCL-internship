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
          p = "P")

manhattan(gwas, annotatePval = 0.01)


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

ggplot(don, aes(x=BPcum, y=-log10(P))) +
  
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

