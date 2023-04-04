#Load packages
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(tidyverse))

#Import descriptive metadata
metadata.csv <-
  read.csv2(here::here("analysis", "data", "raw_data", "metadata.csv"), sep = ";", header = TRUE, na = c("", "NA", "NULL"), encoding = "UTF-8")

#Import nir data, set empty fields to NA
nir.csv <-
  read.csv2(here::here("analysis", "data", "raw_data", "asd_raw_data.csv"), sep = ";", dec = ".", header = TRUE, check.names = FALSE, na = c("","NA","NULL",NULL))

#aggregate observations by group(sample) and calculate average of wavelength measurements
nir.averaged <- 
  aggregate(nir.csv[, 4:2154], list(sample_id = nir.csv$sample_id), mean)

#merge NIR data with metadata
nir.merged <- 
  as.data.frame(merge(metadata.csv, nir.averaged, by='sample_id'))

#Filter xrf data to focus on points and preforms made from quartz/quartzite material 
Points.nir <-
  nir.merged %>%
  dplyr::filter(site_id == "Vilhelmina 1069" | site_id == "Vilhelmina 109" | site_id == "Vilhelmina 112" | site_id == "Vilhelmina 1124" | site_id == "Vilhelmina 1127" | site_id == "Vilhelmina 114" |
           site_id == "Vilhelmina 115" | site_id == "Vilhelmina 117" | site_id == "Vilhelmina 118" | site_id == "Vilhelmina 1254" | site_id == "Vilhelmina 216" | site_id == "Vilhelmina 235" |
           site_id == "Vilhelmina 240" | site_id == "Vilhelmina 245" | site_id == "Vilhelmina 252" | site_id == "Vilhelmina 263" | site_id == "Vilhelmina 335" | site_id == "Vilhelmina 356" |
           site_id == "Vilhelmina 399" | site_id == "Vilhelmina 411" | site_id == "Vilhelmina 419" | site_id == "Vilhelmina 439" | site_id == "Vilhelmina 444" | site_id == "Vilhelmina 450" |
           site_id == "Vilhelmina 458" | site_id == "Vilhelmina 539" | site_id == "Vilhelmina 542" | site_id == "Vilhelmina 611" | site_id == "Vilhelmina 619" | site_id == "Vilhelmina 636" | 
           site_id == "Vilhelmina 637" | site_id == "Vilhelmina 643" | site_id == "Vilhelmina 769" | site_id == "Vilhelmina 949" | site_id == "Vilhelmina 95" | site_id == "Åsele 101" | 
           site_id == "Åsele 107" | site_id == "Åsele 115" | site_id == "Åsele 117" | site_id == "Åsele 119" | site_id == "Åsele 129" | site_id == "Åsele 182" | site_id == "Åsele 188" |
           site_id == "Åsele 393" | site_id == "Åsele 56" | site_id == "Åsele 91" | site_id == "Åsele 92" | site_id == "Åsele 99", 
         type == "Point" | type == "Point fragment" | type == "Preform", 
         material == "Brecciated quartz" | material == "Quartz" | material == "Quartzite") %>% 
  filter(!sample_id %in% c("153","167","168","169","172","174","175","177","182","183","190","191","193","194","196","198","200","204","207","210","213","214",
                           "215","216","229","234","235","237","238","251","262","265","268","269","272","278","281","282","359","377","385","392","393","397","405",
                           "406","410","411","413","414","415","416","417","424","425","426","428","430","432","55","56")) %>% 
  replace_na(list(munsell_hue = "Colourless"))

#perform PCA with SNV normalization and mean-center
nir.pca <-
  prcomp(Points.nir[,c(681:2180)], center = TRUE, scale = FALSE)

#Prepare axis labels with variance in %
nir.pc1lab <- as.data.frame(paste0("PC 1 (",as.character(round(summary(nir.pca)$importance[2,1]*100, digits=2)),"%)"))
nir.pc2lab <- as.data.frame(paste0("PC 2 (",as.character(round(summary(nir.pca)$importance[2,2]*100, digits=2)),"%)"))
nir.pc3lab <- as.data.frame(paste0("PC 3 (",as.character(round(summary(nir.pca)$importance[2,3]*100, digits=2)),"%)"))

#extract loadings from pca
loadings <- 
  as.data.frame(nir.pca$rotation)
#wavelengths are stored as rownames, add them to new column for easy visualization
loadings$wavelengths <- 
  as.numeric(row.names(loadings))

#prepare plots
l1 <-
  ggplot(loadings,aes(x = as.numeric(wavelengths))) + 
    geom_line(aes(y = PC1), linewidth = 1, stat = "identity") +
    labs(y = nir.pc1lab) +
    scale_x_continuous(limits = c(1000, 2500), breaks = scales::pretty_breaks(n = 10)) +
    theme_classic() +
    theme(legend.position = c(.9,.95),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 12, face = "bold", colour = "black"),
          legend.title = element_text(size = 12, face = "bold", colour = "black"))

l2 <- 
  ggplot(loadings,aes(x = as.numeric(wavelengths))) + 
    geom_line(aes(y = PC2), linewidth = 1, stat = "identity") +
    labs(y = nir.pc2lab) +
    scale_x_continuous(limits = c(1000, 2500), breaks = scales::pretty_breaks(n = 10)) +
    theme_classic() +
    theme(legend.position = c(.9,.95),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 12, face = "bold", colour = "black"),
          legend.title = element_text(size = 12, face = "bold", colour = "black"))

l3 <- 
  ggplot(loadings,aes(x = as.numeric(wavelengths))) + 
    geom_line(aes(y = PC3), linewidth = 1, stat = "identity") +
    labs(y = nir.pc3lab,
         x = "Wavelength (nm)") +
    scale_x_continuous(limits = c(1000, 2500), breaks = scales::pretty_breaks(n = 10)) +
    theme_classic() +
    theme(legend.position = c(.9,.95),
          axis.title.y = element_text(size = 12, face = "bold", colour = "black"),
          legend.title = element_text(size = 12, face = "bold", colour = "black"))

#Layout the plots in one figure
fig <-
  ggpubr::ggarrange(l1, l2, l3, 
                    ncol = 1, 
                    nrow = 3,
                    hjust = -1.5)

ggsave("004-nir-pca-loadings.jpeg",
       fig,
       device = "jpeg",
       here::here("analysis", "figures"),
       scale = 1, 
       width=25, 
       height=15,
       units = "cm",
       dpi = 300)