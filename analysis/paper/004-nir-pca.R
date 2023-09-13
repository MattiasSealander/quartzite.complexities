#Load packages
suppressPackageStartupMessages(library(factoextra))
suppressPackageStartupMessages(library(tidyverse))

#Import descriptive metadata
metadata.csv <-
  read.csv2(here::here("analysis", "data", "raw_data", "metadata.csv"), colClasses=c("sample_id"="character"), dec = ".", sep = ";", header = TRUE, na = c("", "NA", "NULL"), encoding = "UTF-8")

#Import nir data, set empty fields to NA
nir.csv <-
  read.csv2(here::here("analysis", "data", "raw_data", "asd_raw_data.csv"), sep = ";", dec = ".", header = TRUE, check.names = FALSE, na = c("","NA","NULL",NULL))

#merge NIR data with metadata
nir.merged <-
  as.data.frame(merge(metadata.csv, nir.csv, by='sample_id'))

#Filter NIR data to focus on points and preforms made from quartz/quartzite material
Points.nir <-
  nir.merged %>%
  filter(
         type == "Point" | type == "Point fragment" | type == "Preform",
         material == "Brecciated quartz" | material == "Quartzite") %>%
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
  ggplot(cbind(basic_plot1$data, Points.nir[, c(11,14)]),
         aes(x=x, y=y, shape = material, fill = hue)) +
  #remove hastag to add sample id as text
  #geom_text(aes(label=Points.nir$sample_id, hjust=0.5,vjust=-1.0)) +
  geom_point(size=3.5) +
  theme_bw() +
  ggtitle("Near infrared 1 000 - 2 500 nm") +
  labs(x = nir.pc1lab,
       y = nir.pc2lab) +
  scale_shape_manual(name = "Material",
                     values=c(25,21)) +
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
  #stat_ellipse(aes(color = hue, group = hue),
  #             linetype = 2,
  #             lwd = 1.2) +
  #scale_color_manual(name = "Group",
  #                   values = c("#36454F", "#E69F00"),
  #                   limits = c("Dark", "Light")) +
  theme(plot.title = element_text(size = 12, face = "bold", colour = "black", vjust = 1, hjust = 0.02),
        axis.title.x = element_text(size = 12, face = "bold", colour = "black"),
        axis.title.y = element_text(size = 12, face = "bold", colour = "black"),
        legend.title = element_text(size = 12, face = "bold", colour = "black"),
        legend.background = element_rect(linetype = "solid", color = "black"),
        legend.position = "bottom")

fig.2 <-
  ggplot(cbind(basic_plot2$data, Points.nir[, c(11,14)]),
         aes(x=x, y=y, shape = material, fill = hue)) +
  #to add sample id as text
  #geom_text(aes(label=Points.nir$sample_id, hjust=0.5,vjust=-1.0)) +
  geom_point(size=3.5) +
  theme_bw() +
  ggtitle("Near infrared 1 000 - 2 500 nm") +
  labs(x = nir.pc2lab,
       y = nir.pc3lab) +
  scale_shape_manual(name = "Material",
                     values=c(25,21)) +
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
  #remove hashtag and run this to add ellipses for group of choice
  #stat_ellipse(aes(color = hue, group = hue),
  #             linetype = 2,
  #             lwd = 1.2) +
  #scale_color_manual(name = "Group",
  #                   values = c("#36454F", "#E69F00"),
  #                   limits = c("Dark", "Light")) +
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
