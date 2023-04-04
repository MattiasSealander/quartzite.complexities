#Load packages
suppressPackageStartupMessages(library(factoextra))
suppressPackageStartupMessages(library(tidyverse))

#Import descriptive metadata
metadata.csv <-
  read.csv2(here::here("analysis", "data", "raw_data", "metadata.csv"), sep = ";", header = TRUE, na = c("", "NA", "NULL"), encoding = "UTF-8")

#Import nir data, set empty fields to NA
nir.csv <-
  read.csv2(here::here("analysis", "data", "raw_data", "asd_raw_data.csv"), sep = ";", dec = ".", header = TRUE, check.names = FALSE, na = c("","NA","NULL",NULL))

#merge NIR data with metadata
nir.merged <- 
  as.data.frame(merge(metadata.csv, nir.csv, by='sample_id'))

#Filter NIR data to focus on points and preforms made from quartz/quartzite material 
Points.nir <-
  nir.merged %>%
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
  filter(!sample_id %in% c("153","167","168","169","172","174","175","176","177","182","190","191","193","194","196","198","200","204","214",
                           "215","216","229","238","262","265","272","282","359","391","392","393","397","405",
                           "408","410","411","413","414","415","416","417","424","425","428","430","432","55","56")) %>% 
  replace_na(list(munsell_hue = "Colourless")) %>% 
  group_by(across(sample_id:weight_g)) %>% 
  summarise(across(`350.0`:`2500.0`, mean), .groups = "drop")

#perform PCA with SNV normalization and mean-center
nir.pca <-
  prcomp(Points.nir[,c(681:2180)], center = TRUE, scale = FALSE)

#prepare variance for call in paper.rmd
nir.pc1var <- round(summary(nir.pca)$importance[2,1]*100, digits=2)
nir.pc2var <- round(summary(nir.pca)$importance[2,2]*100, digits=2)
nir.pc3var <- round(summary(nir.pca)$importance[2,3]*100, digits=2)

#Prepare axis labels with variance in %
nir.pc1lab <- as.data.frame(paste0("PC1 (",as.character(round(summary(nir.pca)$importance[2,1]*100, digits=2)),"%)"))
nir.pc2lab <- as.data.frame(paste0("PC2 (",as.character(round(summary(nir.pca)$importance[2,2]*100, digits=2)),"%)"))
nir.pc3lab <- as.data.frame(paste0("PC3 (",as.character(round(summary(nir.pca)$importance[2,3]*100, digits=2)),"%)"))

pca.colors <- c("#56B4E9", "#36454F", "#E69F00", "#F9F6EE")
pca.hue <- c("Colourless", "Dark", "Light", "White")

#prepare a basic score plot with fviz_pca_ind using PC1 and PC2
basic_plot1 <-
  fviz_pca_ind(nir.pca, axes = c(1,2), label="none")
basic_plot2 <-
  fviz_pca_ind(nir.pca, axes = c(2,3), label="none")

#bind the basic fviz plot for PC 1 and 2, and use as basis for a more customizeable plot in ggpplot
fig.1 <- 
  ggplot(cbind(basic_plot1$data, Points.nir[, c(9,12)]),
         aes(x=x, y=y, shape = material, fill = hue)) +
  #to add sample id as text
  #geom_text(aes(label=Points.nir$sample_id, hjust=0.5,vjust=-1.0)) +
  geom_point(size=3.5) +
  theme_bw() +
  ggtitle("Near infrared 1 000 - 2 500 nm") +
  labs(x = nir.pc1lab,
       y = nir.pc2lab) +
  scale_shape_manual(name = "Material",
                     values=c(25,24,21)) +
  scale_fill_manual(name = "Hue",
                    values=pca.colors) +
  guides(fill = guide_legend(override.aes = list(shape = 22,
                                                 ncol = 1,
                                                 color = "black",
                                                 size=3),
                             title.position="top", 
                             title.hjust = 0.5,
                             order = 1),
         shape = guide_legend(override.aes = list(fill = "black"), 
                              ncol = 1,
                              title.position="top", 
                              title.hjust = 0.5,
                              order = 2)) +
  theme(plot.title = element_text(size = 12, face = "bold", colour = "black", vjust = 1, hjust = 0.02),
        axis.title.x = element_text(size = 12, face = "bold", colour = "black"),
        axis.title.y = element_text(size = 12, face = "bold", colour = "black"),
        legend.title = element_text(size = 12, face = "bold", colour = "black"),
        legend.background = element_rect(linetype = "solid", color = "black"),
        legend.position = "bottom")

fig.2 <- 
  ggplot(cbind(basic_plot2$data, Points.nir[, c(9,12)]),
         aes(x=x, y=y, shape = material, fill = hue)) + 
  #to add sample id as text
  #geom_text(aes(label=Points.nir$sample_id, hjust=0.5,vjust=-1.0)) +
  geom_point(size=3.5) +
  theme_bw() +
  ggtitle("Near infrared 1 000 - 2 500 nm") +
  labs(x = nir.pc2lab,
       y = nir.pc3lab) +
  scale_shape_manual(name = "Material",
                     values=c(25,24,21)) +
  scale_fill_manual(name = "Hue",
                    values=pca.colors) +
  guides(fill = guide_legend(override.aes = list(shape = 22,
                                                 ncol = 1,
                                                 color = "black",
                                                 size=3),
                             title.position="top", 
                             title.hjust = 0.5,
                             order = 1),
         shape = guide_legend(override.aes = list(fill = "black"), 
                              ncol = 1,
                              title.position="top", 
                              title.hjust = 0.5,
                              order = 2)) +
  theme(plot.title = element_blank(),
        axis.title.x = element_text(size = 12, face = "bold", colour = "black"),
        axis.title.y = element_text(size = 12, face = "bold", colour = "black"),
        legend.title = element_text(size = 12, face = "bold", colour = "black"),
        legend.background = element_rect(linetype = "solid", color = "black"),
        legend.position = "bottom")

#Layout the plots in one figure
fig <-
  ggpubr::ggarrange(fig.1,fig.2, 
                    ncol = 1, 
                    nrow = 2,
                    labels = c("A", "B"),
                    common.legend = TRUE,
                    legend = "right")

ggsave("004-nir-pca.jpeg",
       fig,
       device = "jpeg",
       here::here("analysis", "figures"),
       scale = 1, 
       width=25, 
       height=25,
       units = "cm",
       dpi = 300)