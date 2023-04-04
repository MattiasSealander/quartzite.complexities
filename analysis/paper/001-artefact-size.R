# Loading package
library(ggplot2)
suppressPackageStartupMessages(library(tidyverse))

#Import descriptive metadata
metadata.csv <-
  read.csv2(here::here("analysis", "data", "raw_data", "metadata.csv"), sep = ";", header = TRUE, na = c("", "NA", "NULL"), encoding = "UTF-8")

#Filter data to focus on points and preforms made from quartz/quartzite material in Vilhelmina and Åsele parish
Points <-
  metadata.csv %>% 
  dplyr::filter(site_id == "Vilhelmina 1069" | site_id == "Vilhelmina 109" | site_id == "Vilhelmina 112" | site_id == "Vilhelmina 1124" | site_id == "Vilhelmina 1127" | 
                  site_id == "Vilhelmina 114" | site_id == "Vilhelmina 115" | site_id == "Vilhelmina 117" | site_id == "Vilhelmina 118" | site_id == "Vilhelmina 1254" | 
                  site_id == "Vilhelmina 216" | site_id == "Vilhelmina 235" | site_id == "Vilhelmina 240" | site_id == "Vilhelmina 245" | site_id == "Vilhelmina 252" | 
                  site_id == "Vilhelmina 263" | site_id == "Vilhelmina 335" | site_id == "Vilhelmina 356" | site_id == "Vilhelmina 399" | site_id == "Vilhelmina 411" | 
                  site_id == "Vilhelmina 419" | site_id == "Vilhelmina 439" | site_id == "Vilhelmina 444" | site_id == "Vilhelmina 450" | site_id == "Vilhelmina 458" | 
                  site_id == "Vilhelmina 539" | site_id == "Vilhelmina 542" | site_id == "Vilhelmina 611" | site_id == "Vilhelmina 619" | site_id == "Vilhelmina 636" | 
                  site_id == "Vilhelmina 637" | site_id == "Vilhelmina 643" | site_id == "Vilhelmina 769" | site_id == "Vilhelmina 949" | site_id == "Vilhelmina 95" | 
                  site_id == "Åsele 101" | site_id == "Åsele 107" | site_id == "Åsele 115" | site_id == "Åsele 117" | site_id == "Åsele 119" | site_id == "Åsele 129" | 
                  site_id == "Åsele 182" | site_id == "Åsele 188" | site_id == "Åsele 393" | site_id == "Åsele 56" | site_id == "Åsele 91" | site_id == "Åsele 92" | 
                  site_id == "Åsele 99", 
                type == "Point" | type == "Point fragment" | type == "Preform", 
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
       width=15, 
       height=10,
       units = "cm",
       dpi = 300)
