# Loading package
suppressPackageStartupMessages(library(data.table))
library(gridBase)
suppressPackageStartupMessages(library(tidyverse))

#Import descriptive metadata
metadata.csv <-
  read.csv2(here::here("analysis", "data", "raw_data", "metadata.csv"), dec = ".", sep = ";", header = TRUE, na = c("", "NA", "NULL"), encoding = "UTF-8")

#Import nir data, set empty fields to NA
nir.csv <-
  read.csv2(here::here("analysis", "data", "raw_data", "asd_raw_data.csv"), sep = ";", dec = ".", header = TRUE, check.names = FALSE, na = c("","NA","NULL",NULL))

#merge NIR data with metadata
nir.merged <-
  as.data.frame(merge(metadata.csv, nir.csv, by='sample_id'))

#Filter xrf data to focus on points and preforms made from quartz/quartzite material
Points.nir <-
  nir.merged %>%
  dplyr::filter(type == "Point" | type == "Point fragment" | type == "Preform",
                material == "Brecciated quartz" | material == "Quartz" | material == "Quartzite") %>%
  dplyr::filter(!sample_id %in% c("153","167","168","169","172","174","175","177","182","183","190","191","193","194","196","198","200","204","207","210","213","214",
                                  "215","216","229","234","235","237","238","251","262","265","268","269","272","278","281","282","359","377","385","392","393","397","405",
                                  "406","410","411","413","414","415","416","417","424","425","426","428","430","432","55","56")) %>%
  replace_na(list(munsell_hue = "Colourless")) %>%
  group_by(across(sample_id:weight_g)) %>%
  dplyr::summarise(across(`350.0`:`2500.0`, mean), .groups = "drop")

#Filter NIR data to focus on material with dark hues
#and select the NIR range 1 001 - 2 500 nm, exclude 1000 nm due to filter shift
Points.d <- Points.nir %>%
  dplyr::filter(hue == "Dark") %>%
  dplyr::select(`1001.0`:`2500.0`) %>%
  summarise_if(is.numeric, mean)

#Melt into long format
Points.d <-
  suppressWarnings(data.table::melt(setDT(Points.d), variable.name = "Wavelength", variable.factor = FALSE, value.name = "Absorbance"))

#Filter NIR data to focus on material with light hues
#and select the NIR range 1 001 - 2 500 nm, exclude 1000 nm due to filter shift
Points.l <- Points.nir %>%
  dplyr::filter(hue == "Light") %>%
  dplyr::select(`1001.0`:`2500.0`) %>%
  summarise_if(is.numeric, mean)

#Melt into long format
Points.l <-
  suppressWarnings(data.table::melt(setDT(Points.l), variable.name = "Wavelength", variable.factor = FALSE, value.name = "Absorbance"))

#Filter NIR data to focus on colourless material
#and select the NIR range 1 001 - 2 500 nm, exclude 1000 nm due to filter shift
Points.c <- Points.nir %>%
  dplyr::filter(hue == "Colourless") %>%
  dplyr::select(`1001.0`:`2500.0`) %>%
  summarise_if(is.numeric, mean)

#Melt into long format
Points.c <-
  suppressWarnings(data.table::melt(setDT(Points.c), variable.name = "Wavelength", variable.factor = FALSE, value.name = "Absorbance"))

#Filter NIR data to focus on material with a white hue
#and select the NIR range 1 001 - 2 500 nm, exclude 1000 nm due to filter shift
Points.w <- Points.nir %>%
  dplyr::filter(hue == "White") %>%
  dplyr::select(`1001.0`:`2500.0`) %>%
  summarise_if(is.numeric, mean)

#Melt into long format
Points.w <-
  suppressWarnings(data.table::melt(setDT(Points.w), variable.name = "Wavelength", variable.factor = FALSE, value.name = "Absorbance"))

p.d <-
  ggplot(Points.d, aes(x = as.numeric(Wavelength))) +
    geom_line(aes(y = Absorbance, colour = ""), linewidth = 1, stat = "identity") +
    geom_vline(xintercept=c(1065,1296,1413,1935,2220,2260,2350), linetype='dashed', col = 'red') +
  ggtext::geom_richtext(data=data.frame(),
                        aes(x = c(1065,1300), y = c(0.99, 0.99), label="Fe<sup>+</sup>?"), label.size = NA, fontface = 2) +
  ggtext::geom_richtext(data=data.frame(),
                        aes(x = c(1413,1935), y = c(0.995, 1.005), label=c("OH", "OH")), label.size = NA, fontface = 2) +
  ggtext::geom_richtext(data=data.frame(),
                        aes(x = c(2220,2260,2350), y = c(1.0125, 1.018, 1.018), label=c("AlOH / MgOH", "AlOH", "AlOH")), label.size = NA, fontface = 2, size = 3) +
  ylab("Dark") +
  scale_color_manual(name = "Dark",
                     values = "#36454F") +
  scale_x_continuous(limits = c(1000, 2500), breaks = scales::pretty_breaks(n = 10)) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y = element_text(size = 12, face = "bold", colour = "black"))

p.l <-
  ggplot(Points.l, aes(x = as.numeric(Wavelength))) +
  geom_line(aes(y = Absorbance, colour = ""), linewidth = 1, stat = "identity") +
  geom_vline(xintercept=c(1065,1296,1413,1935,2220,2260,2350), linetype='dashed', col = 'red') +
  ylab("Light") +
  scale_color_manual(name = "Light",
                     values = "#E69F00") +
  scale_x_continuous(limits = c(1000, 2500), breaks = scales::pretty_breaks(n = 10)) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y = element_text(size = 12, face = "bold", colour = "black"))

p.c <-
  ggplot(Points.c, aes(x = as.numeric(Wavelength))) +
  geom_line(aes(y = Absorbance, colour = ""), linewidth = 1, stat = "identity") +
  geom_vline(xintercept=c(1065,1296,1413,1935,2220,2260,2350), linetype='dashed', col = 'red') +
  ylab("Colourless") +
  scale_color_manual(name = "Colourless",
                     values = "#56B4E9") +
  scale_x_continuous(limits = c(1000, 2500), breaks = scales::pretty_breaks(n = 10)) +
  scale_y_continuous(labels = function(x) sprintf("%.2f", x)) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y = element_text(size = 12, face = "bold", colour = "black"))


p.w <-
  ggplot(Points.w, aes(x = as.numeric(Wavelength))) +
  geom_line(aes(y = Absorbance, colour = ""), linewidth = 1, stat = "identity") +
  geom_vline(xintercept=c(1065,1296,1413,1935,2220,2260,2350), linetype='dashed', col = 'red') +
  xlab("Wavelength (nm)") +
  ylab("White") +
  scale_color_manual(name = "White",
                     values = "#CC79A7") +
  scale_x_continuous(limits = c(1000, 2500), breaks = scales::pretty_breaks(n = 10)) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 12, face = "bold", colour = "black"),
        axis.title.y = element_text(size = 12, face = "bold", colour = "black"))

#Layout the plots in one figure
fig <-
  ggpubr::ggarrange(p.d, p.l, p.c, p.w,
                    ncol = 1,
                    nrow = 4)

#Save figure
ggsave("002-nir-peak-summary.jpeg",
       fig,
       device = "jpeg",
       here::here("analysis", "figures"),
       width=25,
       height=25,
       units = "cm",
       dpi = 300)
