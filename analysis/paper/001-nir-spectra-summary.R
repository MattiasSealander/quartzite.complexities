# Loading package
library(data.table)
library(ggplot2)
library(gridBase)
suppressPackageStartupMessages(library(tidyverse))

#Import descriptive metadata
metadata.csv <-
  read.csv2(here::here("analysis", "data", "raw_data", "metadata.csv"), sep = ";", dec = ".", header = TRUE, na = c("", "NA", "NULL"), encoding = "UTF-8")

#Import nir data, set empty fields to NA
nir.csv <-
  read.csv2(here::here("analysis", "data", "raw_data", "asd_raw_data.csv"), sep = ";", dec = ".", header = TRUE, check.names = FALSE, na = c("","NA","NULL",NULL))

#merge NIR data with metadata
nir.merged <-
  as.data.frame(merge(metadata.csv, nir.csv, by='sample_id'))

#Filter nir data to focus on points and preforms made from quartz/quartzite material
Points.nir <-
  nir.merged %>%
  dplyr::filter(type == "Point" | type == "Point fragment" | type == "Preform",
                material == "Brecciated quartz" | material == "Quartz" | material == "Quartzite") %>%
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
  suppressWarnings(data.table::melt(setDT(n.colourless), id.vars = "sample_id", variable.name = "Wavelength", variable.factor = FALSE, value.name = "Absorbance"))

#Filter NIR data to focus on material with dark hues
#and select the NIR range 1 001 - 2 500 nm, exclude 1000 nm due to filter shift
n.dark <- Points.nir %>%
  dplyr::filter(hue == "Dark") %>%
  dplyr::select(sample_id, `1001.0`:`2500.0`)

#Melt into long format
n.dark.long <-
  suppressWarnings(data.table::melt(setDT(n.dark), id.vars = "sample_id", variable.name = "Wavelength", variable.factor = FALSE, value.name = "Absorbance"))

#Filter NIR data to focus on material with light hues
#and select the NIR range 1 001 - 2 500 nm, exclude 1000 nm due to filter shift
n.light <- Points.nir %>%
  dplyr::filter(hue == "Light") %>%
  dplyr::select(sample_id, `1001.0`:`2500.0`)

#Melt into long format
n.light.long <-
  suppressWarnings(data.table::melt(setDT(n.light), id.vars = "sample_id", variable.name = "Wavelength", variable.factor = FALSE, value.name = "Absorbance"))

#Filter NIR data to focus on material with white hues
#and select the NIR range 1 001 - 2 500 nm, exclude 1000 nm due to filter shift
n.white <- Points.nir %>%
  dplyr::filter(hue == "White") %>%
  dplyr::select(sample_id, `1001.0`:`2500.0`)

#Melt into long format
n.white.long <-
  suppressWarnings(data.table::melt(setDT(n.white), id.vars = "sample_id", variable.name = "Wavelength", variable.factor = FALSE, value.name = "Absorbance"))

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
  geom_line(aes(y = Absorbance, colour = "", group = sample_id), linewidth = 0.5, stat = "identity") +
  geom_text(data=data.frame(), aes(x=1100,y=1.55,label="Dark"), size=5, fontface=2) +
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
ggsave("001-nir-spectra-summary.jpeg",
       fig,
       device = "jpeg",
       here::here("analysis", "figures"),
       width=25,
       height=20,
       units = "cm",
       dpi = 300)
