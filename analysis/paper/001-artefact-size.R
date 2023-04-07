# Loading package
library(ggplot2)
suppressPackageStartupMessages(library(tidyverse))

#Import descriptive metadata
metadata.csv <-
  read.csv2(here::here("analysis", "data", "raw_data", "metadata.csv"), sep = ";", dec = ".", header = TRUE, na = c("", "NA", "NULL"), encoding = "UTF-8")

#Filter data to focus on points and preforms made from quartz/quartzite material in Vilhelmina and Åsele parish
Points <-
  metadata.csv %>%
  dplyr::filter(type == "Point" | type == "Point fragment" | type == "Preform",
                material == "Brecciated quartz" | material == "Quartz" | material == "Quartzite")  %>%
  replace_na(list(munsell_hue = "Colourless"))

#plot the length and width
fig <-
  Points %>%
  ggplot() +
  geom_point(shape = 21, size = 2, aes(y = max_width_mm, x = max_length_mm, fill = type)) +
  xlab("Length (mm)") +
  ylab("Width (mm)") +
  scale_fill_manual(name = "Type",
                    values=ggpubfigs::friendly_pal("contrast_three")) +
  theme_minimal() +
  theme(legend.position = c(.9,.1),
        axis.title.x = element_text(size = 12, face = "bold", colour = "black"),
        axis.title.y = element_text(size = 12, face = "bold", colour = "black"),
        legend.title = element_text(size = 12, face = "bold", colour = "black"))

#Save figure
ggsave("001-artefact-size.jpeg",
       fig,
       device = "jpeg",
       here::here("analysis", "figures"),
       width=20,
       height=15,
       units = "cm",
       dpi = 300)
