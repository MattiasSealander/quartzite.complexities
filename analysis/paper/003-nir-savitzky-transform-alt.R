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

#aggregate observations by group(sample) and calculate average of wavelength measurements
nir.averaged <- 
  aggregate(nir.csv[, 4:2154], list(sample_id = nir.csv$sample_id), mean)

#merge NIR data with metadata
nir.merged <- 
  as.data.frame(merge(metadata.csv, nir.averaged, by='sample_id'))

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
  filter(!sample_id %in% c("153","167","168","169","172","174","175","176","177","182","183","190","191","193","194","196","198","200","204","207","210","213","214",
                           "215","216","229","234","235","237","238","251","262","265","268","272","278","281","282","359","377","392","393","397","405",
                           "406","410","411","413","414","415","416","417","424","425","426","428","430","55","56")) %>% 
  replace_na(list(munsell_hue = "Colourless"))


#Filter NIR data to focus on light and dark hue material
#and select the NIR range 1 000 - 2 500 nm, as well as the hue variable
Points <- Points.nir %>%
  filter(hue == "Light" | hue == "Dark") %>% 
  dplyr::select(hue, `1001.0`:`2500.0`)

#Reorder the hue factor level so light spectra will be drawn behind dark spectra
Points$hue <- factor(Points$hue, c("Light", "Dark"))

#Perform Savitzky-Golay 2nd derivative transformation, start at 1001 nm due to filter shift at 1000 nm
sg <- cbind(Points$hue, as.data.table(gapDer(X = Points[,2:1501], m = 2, w = 11, s = 5)))

#Transpose the data for ggplot
sg.long <-  suppressWarnings(data.table::melt(setDT(sg), id.vars = "V1", variable.name = "Wavelength", value.name = "Absorbance", variable.factor = FALSE))

#plot the transformed spectra
fig <-
  ggplot(sg.long, aes(x = as.numeric(Wavelength))) + 
  geom_line(aes(y = as.numeric(Absorbance), color = V1), linewidth = 0.1, stat = "identity") +
  xlab("Wavelength (nm)") +
  ylab("Absorbance") +
  scale_color_manual(name = "Hue",
                    values = c("#E69F00", "#36454F"),
                    breaks = c("Light","Dark")) + 
  scale_x_continuous(limits = c(1000, 2500), breaks = scales::pretty_breaks(n = 10)) +
  theme_classic() +
  theme(legend.position = "bottom",
        axis.title.x = element_text(size = 12, face = "bold", colour = "black"),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 12, face = "bold", colour = "black"))

#Save figure
ggsave("003-nir-savitzky-transform-alt.jpeg",
       fig,
       device = "jpeg",
       here::here("analysis", "figures"),
       width=25, 
       height=15,
       units = "cm",
       dpi = 300)