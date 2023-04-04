#Load packages
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(prospectr))
suppressPackageStartupMessages(library(reshape2))
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
  replace_na(list(munsell_hue = "Colourless")) %>% 
  group_by(across(sample_id:weight_g)) %>% 
  dplyr::summarise(across(`350.0`:`2500.0`, mean), .groups = "drop")

#Filter NIR data to focus on sample 258 and 411, then select the NIR range 1 001 - 2 500 nm, exclude 1000 nm due to filter shift
raw.spec <- Points.nir %>%
  filter(sample_id == "258" | sample_id == "411") %>% 
  dplyr::select(sample_id, `1001.0`:`2500.0`) 

#Melt into long format, data.table needs to be specified to avoid r using reshape2 version that sets variable to factor
raw.spec.long <- 
  suppressWarnings(data.table::melt(setDT(raw.spec), variable.name = "Wavelength", variable.factor = FALSE, value.name = "Absorbance"))

#plot the raw spectra for sample 258, and 411
p.raw.spec <-
  raw.spec.long %>% 
  ggplot(aes(x = as.numeric(Wavelength))) + 
  geom_line(aes(y = Absorbance, colour = "", group = sample_id), linewidth = 1, stat = "identity") +
  geom_text(data = . %>%  filter(Wavelength == max(Wavelength)), aes(x=Inf, y = Absorbance, fontface = "bold"),
            hjust = 1.0, size = 4, label = raw.spec$sample_id) +
  coord_cartesian(clip = "off") +
  xlab("Wavelength (nm)") +
  ylab("Absorbance") +
  scale_color_manual(values = "black") + 
  scale_x_continuous(limits = c(1000, 2500), breaks = scales::pretty_breaks(n = 10)) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

#Filter spectra to focus on sample 411, which has more noise
sg.411 <- raw.spec %>%
  filter(sample_id == "411")

#Perform Savitzky-Golay 2nd derivative on sample 411
sg.411 <- as.data.table(gapDer(X = sg.411[,2:1501], m = 2, w = 11, s = 5))

#Transpose the data for ggplot
sg.411.long <- suppressWarnings(data.table::melt(setDT(sg.411), variable.name = "Wavelength", value.name = "Absorbance", variable.factor = FALSE))

#plot the transformed dark spectra for sample 411
p.sg.411 <-
  sg.411.long %>% 
  ggplot() + 
  geom_line(aes(x = as.numeric(Wavelength), y = as.numeric(Absorbance), colour = ""), linewidth = 1, stat = "identity") +
  xlab("Wavelength (nm)") +
  ylab("Absorbance") +
  geom_text(data=sg.411, aes(label = '411', x = 2500, y = 1.25e-05, fontface = "bold")) +
  scale_color_manual(values = "black") + 
  scale_x_continuous(limits = c(1000, 2500), breaks = scales::pretty_breaks(n = 10)) +
  #setting breaks for y axis in order to prevent ggarrange label covering the axis values
  scale_y_continuous(limits = c(-1.1e-05, 1.5e-05), breaks = seq(-1.0e-05, 1.5e-05, by = 0.75e-5)) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 12, face = "bold", colour = "black"),
        axis.title.y = element_blank())

#Filter spectra to focus on sample 258, which has less noise
sg.258 <- raw.spec %>%
  filter(sample_id == "258")

#Perform Savitzky-Golay 2nd derivative on sample 258
sg.258 <- as.data.table(gapDer(X = sg.258[,2:1501], m = 2, w = 11, s = 5))

#Transpose the data for ggplot
sg.258.long <-  suppressWarnings(data.table::melt(setDT(sg.258), variable.name = "Wavelength", value.name = "Absorbance", variable.factor = FALSE))

#plot the transformed dark spectra for sample 258
p.sg.258 <-
  sg.258.long %>% 
  ggplot() + 
  geom_line(aes(x = as.numeric(Wavelength), y = as.numeric(Absorbance), colour = ""), size = 1, stat = "identity") +
  xlab("Wavelength (nm)") +
  ylab("Absorbance") +
  geom_text(data=sg.258, aes(label = '258', x = 2500, y = 2.5e-05, fontface = "bold")) +
  scale_color_manual(values = "black") + 
  scale_x_continuous(limits = c(1000, 2500), breaks = scales::pretty_breaks(n = 10)) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

#Layout the plots in one figure
fig <-
  ggpubr::ggarrange(p.raw.spec, p.sg.258, p.sg.411, labels = c("A", "B", "C"), nrow = 3, align="v")

#Save figure
ggsave("003-nir-savitzky-transform.jpeg",
       fig,
       device = "jpeg",
       here::here("analysis", "figures"),
       width=25, 
       height=20,
       units = "cm",
       dpi = 300)