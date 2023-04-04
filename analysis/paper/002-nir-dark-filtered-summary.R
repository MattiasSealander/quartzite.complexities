# Loading package
library(gridBase)
library(data.table)
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

#Filter xrf data to focus on points and preforms made from quartz/quartzite material 
Points.nir <-
  nir.merged %>%
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
                material == "Brecciated quartz" | material == "Quartz" | material == "Quartzite") %>% 
  dplyr::filter(!sample_id %in% c("153","167","168","169","172","174","175","177","182","183","190","191","193","194","196","198","200","204","207","210","213","214",
                           "215","216","229","234","235","237","238","251","262","265","268","269","272","278","281","282","359","377","385","392","393","397","405",
                           "406","410","411","413","414","415","416","417","424","425","426","428","430","432","55","56")) %>% 
  replace_na(list(munsell_hue = "Colourless")) %>% 
  group_by(across(sample_id:weight_g)) %>% 
  dplyr::summarise(across(`350.0`:`2500.0`, mean), .groups = "drop")

#Filter NIR data to focus on colourless 
#and select the NIR range 1 001 - 2 500 nm, exclude 1000 nm due to filter shift
n.colourless <- Points.nir %>%
  dplyr::filter(hue == "Colourless") %>% 
  dplyr::select(sample_id, `1001.0`:`2500.0`) 

#Melt into long format
n.colourless.long <- 
  suppressWarnings(data.table::melt(setDT(n.colourless), variable.name = "Wavelength", variable.factor = FALSE, value.name = "Absorbance"))

#Filter NIR data to focus on material with dark hues 
#and select the NIR range 1 001 - 2 500 nm, exclude 1000 nm due to filter shift
n.dark <- Points.nir %>%
  dplyr::filter(hue == "Dark") %>% 
  dplyr::select(sample_id, `1001.0`:`2500.0`) 

#Melt into long format
n.dark.long <- 
  suppressWarnings(data.table::melt(setDT(n.dark), variable.name = "Wavelength", variable.factor = FALSE, value.name = "Absorbance"))

#Filter NIR data to focus on material with light hues 
#and select the NIR range 1 001 - 2 500 nm, exclude 1000 nm due to filter shift
n.light <- Points.nir %>%
  dplyr::filter(hue == "Light") %>% 
  dplyr::select(sample_id, `1001.0`:`2500.0`) 

#Melt into long format
n.light.long <- 
  suppressWarnings(data.table::melt(setDT(n.light), variable.name = "Wavelength", variable.factor = FALSE, value.name = "Absorbance"))

#Filter NIR data to focus on material with white hues 
#and select the NIR range 1 001 - 2 500 nm, exclude 1000 nm due to filter shift
n.white <- Points.nir %>%
  dplyr::filter(hue == "White") %>% 
  dplyr::select(sample_id, `1001.0`:`2500.0`) 

#Melt into long format
n.white.long <- 
  suppressWarnings(data.table::melt(setDT(n.white), variable.name = "Wavelength", variable.factor = FALSE, value.name = "Absorbance"))

#plot the colourless spectra
p.c <-
  ggplot(n.colourless.long, aes(x = as.numeric(Wavelength))) + 
  geom_line(aes(y = Absorbance, colour = "", group = sample_id), linewidth = 1, stat = "identity") +
  geom_text(data=data.frame(), aes(x=2300,y=0.51,label="Colourless"), size=5, fontface=2) +
  xlab("Wavelength (nm)") +
  ylab("Absorbance") +
  scale_color_manual(name = "Colourless",
                     values = "#56B4E9") +
  scale_x_continuous(limits = c(1000, 2500), breaks = scales::pretty_breaks(n = 10)) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, face = "bold", colour = "black"),
        legend.title = element_text(size = 12, face = "bold", colour = "black"))

#plot the dark spectra
p.d <-
  ggplot(n.dark.long, aes(x = as.numeric(Wavelength))) + 
  geom_line(aes(y = Absorbance, colour = "", group = sample_id), linewidth = 1, stat = "identity") +
  geom_text(data=data.frame(), aes(x=1100,y=1.3,label="Dark"), size=5, fontface=2) +
  xlab("Wavelength (nm)") +
  ylab("Absorbance") +
  scale_color_manual(name = "Dark",
                     values = "#36454F") +
  scale_x_continuous(limits = c(1000, 2500), breaks = scales::pretty_breaks(n = 10)) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 12, face = "bold", colour = "black"))

#plot the light spectra
p.l <-
  ggplot(n.light.long, aes(x = as.numeric(Wavelength))) + 
  geom_line(aes(y = Absorbance, colour = "", group = sample_id), linewidth = 1, stat = "identity") +
  geom_text(data=data.frame(), aes(x=1100,y=1.24,label="Light"), size=5, fontface=2) +
  xlab("Wavelength (nm)") +
  ylab("Absorbance") +
  scale_color_manual(name = "Light",
                     values = "#E69F00") +
  scale_x_continuous(limits = c(1000, 2500), breaks = scales::pretty_breaks(n = 10)) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 12, face = "bold", colour = "black"),
        axis.title.y = element_text(size = 12, face = "bold", colour = "black"),
        legend.title = element_text(size = 12, face = "bold", colour = "black"))

#plot the white spectra
p.w <-
  ggplot(n.white.long, aes(x = as.numeric(Wavelength))) + 
  geom_line(aes(y = Absorbance, colour = "", group = sample_id), linewidth = 1, stat = "identity") +
  geom_text(data=data.frame(), aes(x=1100,y=0.386,label="White"), size=5, fontface=2) +
  xlab("Wavelength (nm)") +
  ylab("Absorbance") +
  scale_color_manual(name = "White",
                     values = "#CC79A7") +
  scale_x_continuous(limits = c(1000, 2500), breaks = scales::pretty_breaks(n = 10)) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 12, face = "bold", colour = "black"),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 12, face = "bold", colour = "black"))

#Layout the plots in one figure
fig <-
  ggpubr::ggarrange(p.c, p.d, p.l, p.w, 
                    ncol = 2, 
                    nrow = 2)

#Save figure
ggsave("002-nir-dark-filtered-summary.jpeg",
       fig,
       device = "jpeg",
       here::here("analysis", "figures"),
       width=20, 
       height=20,
       units = "cm",
       dpi = 300)