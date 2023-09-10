#Load packages
suppressPackageStartupMessages(library(factoextra))
suppressPackageStartupMessages(library(prospectr))
suppressPackageStartupMessages(library(tidyverse))

#Import descriptive metadata
metadata.csv <-
  read.csv2(here::here("analysis", "data", "raw_data", "metadata.csv"), colClasses=c("sample_id"="character"), dec = ".", sep = ";", header = TRUE, na = c("", "NA", "NULL"), encoding = "UTF-8")

raman.baseline <-
  read.csv2(here::here("analysis", "data", "derived_data", "raman_baseline_corr.csv"), sep = ",", dec = ".", header = TRUE, check.names = FALSE, na = c("", "NA", "NULL"), encoding = "UTF-8")

#merge raman de-trended data with metadata
raman.baseline.merged <-
  as.data.frame(merge(metadata.csv, raman.baseline, by='sample_id'))

#select and filter raman data to focus on relevant material and artefact types
Points.raman <-
  raman.baseline.merged %>%
  filter(site_id == "Vilhelmina 1069" | site_id == "Vilhelmina 109" | site_id == "Vilhelmina 112" | site_id == "Vilhelmina 1124" | site_id == "Vilhelmina 1127" | site_id == "Vilhelmina 114" |
           site_id == "Vilhelmina 115" | site_id == "Vilhelmina 117" | site_id == "Vilhelmina 118" | site_id == "Vilhelmina 1254" | site_id == "Vilhelmina 216" | site_id == "Vilhelmina 235" |
           site_id == "Vilhelmina 240" | site_id == "Vilhelmina 245" | site_id == "Vilhelmina 252" | site_id == "Vilhelmina 263" | site_id == "Vilhelmina 335" | site_id == "Vilhelmina 356" |
           site_id == "Vilhelmina 399" | site_id == "Vilhelmina 411" | site_id == "Vilhelmina 419" | site_id == "Vilhelmina 439" | site_id == "Vilhelmina 444" | site_id == "Vilhelmina 450" |
           site_id == "Vilhelmina 458" | site_id == "Vilhelmina 539" | site_id == "Vilhelmina 542" | site_id == "Vilhelmina 611" | site_id == "Vilhelmina 619" | site_id == "Vilhelmina 636" |
           site_id == "Vilhelmina 637" | site_id == "Vilhelmina 643" | site_id == "Vilhelmina 769" | site_id == "Vilhelmina 949" | site_id == "Vilhelmina 95" | site_id == "Åsele 101" |
           site_id == "Åsele 107" | site_id == "Åsele 115" | site_id == "Åsele 117" | site_id == "Åsele 119" | site_id == "Åsele 129" | site_id == "Åsele 182" | site_id == "Åsele 188" |
           site_id == "Åsele 393" | site_id == "Åsele 56" | site_id == "Åsele 91" | site_id == "Åsele 92" | site_id == "Åsele 99",
         type == "Point" | type == "Point fragment" | type == "Preform",
         material == "Brecciated quartz" | material == "Quartz" | material == "Quartzite") %>%
  replace_na(list(munsell_hue = "Colourless"))

#perform PCA with SNV normalization and mean-center
raman.pca <-
  prcomp(Points.raman[,c(30:322)], center = TRUE, scale = FALSE)

#prepare labels for PCs, doing it in 2 steps allows variance percentage to be called on in the Rmd file
raman.pc1var <- round(summary(raman.pca)$importance[2,1]*100, digits=1)
raman.pc2var <- round(summary(raman.pca)$importance[2,2]*100, digits=1)
raman.pc3var <- round(summary(raman.pca)$importance[2,3]*100, digits=1)
raman.pc1lab <- paste0("PC1 (",as.character(raman.pc1var),"%)")
raman.pc2lab <- paste0("PC2 (",as.character(raman.pc2var),"%)")
raman.pc3lab <- paste0("PC3 (",as.character(raman.pc3var),"%)")

pca.colors <- c("#56B4E9", "#36454F", "#E69F00", "#F9F6EE")
pca.hue <- c("Colourless", "Dark", "Light", "White")

#prepare a basic score plot with fviz_pca_ind using PC1 and PC2
basic_plot1 <-
  fviz_pca_ind(raman.pca, axes = c(1,2), label="none")


#bind the basic fviz plot for PC 1 and 2, and use as basis for a more customizeable plot in ggpplot
fig <-
  ggplot(cbind(basic_plot1$data, Points.raman[, c(9,12)]),
         aes(x=x, y=y, shape = material, fill = hue)) +
  #to add sample id as text
  #geom_text(aes(label=Points.raman$sample_id, hjust=0.5,vjust=-1.0)) +
  geom_point(size=3.5) +
  theme_bw() +
  ggtitle(bquote(bold("Raman 100 - 1 800 " ~cm^-1))) +
  labs(x = raman.pc1lab,
       y = raman.pc2lab) +
  scale_shape_manual(name = "Material",
                     values=c(25,21)) +
  scale_fill_manual(name = "Hue",
                    values=pca.colors) +
  guides(fill = guide_legend(override.aes = list(shape = 22,
                                                 fill = pca.colors,
                                                 ncol = 4,
                                                 color = "black",
                                                 size=3,
                                                 order = 1),
                             title.position="top", title.hjust = 0.5),
         shape = guide_legend(override.aes = list(fill = "black"),
                              ncol = 3,
                              title.position="top",
                              title.hjust = 0.5,
                              order = 2)) +
  theme(plot.title = element_text(size = 12, colour = "black", vjust = 1, hjust = 0.02),
        axis.title.x = element_text(size = 12, face = "bold", colour = "black"),
        axis.title.y = element_text(size = 12, face = "bold", colour = "black"),
        legend.title = element_text(size = 12, face = "bold", colour = "black"),
        legend.background = element_rect(linetype = "solid", color = "black"),
        legend.position = "bottom")

ggsave("007-raman-pca.jpeg",
       fig,
       device = "jpeg",
       here::here("analysis", "figures"),
       scale = 1,
       width=25,
       height=20,
       units = "cm",
       dpi = 300)
