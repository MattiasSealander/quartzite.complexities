#Load packages
suppressPackageStartupMessages(library(spectrolab))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggtext))

#Import descriptive metadata
metadata.csv <-
  read.csv2(here::here("analysis", "data", "raw_data", "metadata.csv"), colClasses=c("sample_id"="character"), dec = ".", sep = ";", header = TRUE, na = c("", "NA", "NULL"), encoding = "UTF-8")

#Import raman data, set empty fields to NA
raman.csv <-
  read.csv2(here::here("analysis", "data", "raw_data", "raman_raw_data.csv"), sep = ";", dec = ",", header = TRUE, check.names = FALSE, na = c("","NA","NULL",NULL))

#aggregate observations by group(sample) and calculate average of wavelength measurements
raman.averaged <-
  aggregate(raman.csv[, 4:468], list(sample_id = raman.csv$sample_id), mean)

#merge raman de-trended data with metadata
raman.baseline.merged <-
  as.data.frame(merge(metadata.csv, raman.averaged, by='sample_id'))

Points.raman <-
  raman.baseline.merged %>%
  filter(type == "Point" | type == "Point fragment" | type == "Preform",
         material == "Brecciated quartz" | material == "Quartz" | material == "Quartzite") %>%
  replace_na(list(munsell_hue = "Colourless"))

#Filter raman data to focus on dark samples
r.dark <- Points.raman %>%
  filter(hue == "Dark") %>%
  select(sample_id, c(`92.88`:`2503.59`))

#Melt into long format
r.dark.long <-
  suppressWarnings(reshape2::melt(data.table::setDT(r.dark), id.vars = "sample_id", variable.name = "Wavenumber", variable.factor = FALSE, value.name = "Intensity"))

#Filter raman data to focus on light samples
r.light <- Points.raman %>%
  filter(hue == "Light") %>%
  select(sample_id, c(`92.88`:`2503.59`))

#Melt into long format
r.light.long <-
  suppressWarnings(reshape2::melt(data.table::setDT(r.light), id.vars = "sample_id", variable.name = "Wavenumber", variable.factor = FALSE, value.name = "Intensity"))

#Filter raman data to focus on colourless samples
r.colourless <- Points.raman %>%
  filter(hue == "Colourless") %>%
  select(sample_id, c(`92.88`:`2503.59`))

#Melt into long format
r.colourless.long <-
  suppressWarnings(reshape2::melt(data.table::setDT(r.colourless), id.vars = "sample_id", variable.name = "Wavenumber", variable.factor = FALSE, value.name = "Intensity"))

#Filter raman data to focus on white samples
r.white <- Points.raman %>%
  filter(hue == "White") %>%
  select(sample_id, c(`92.88`:`2503.59`))

#Melt into long format
r.white.long <-
  suppressWarnings(reshape2::melt(data.table::setDT(r.white), id.vars = "sample_id", variable.name = "Wavenumber", variable.factor = FALSE, value.name = "Intensity"))

#plot the colourless spectra
p.c <-
  ggplot(r.colourless.long,aes(x = as.numeric(as.character(Wavenumber)))) +
    geom_line(aes(y = Intensity, colour = "", group = sample_id), linewidth = 1, stat = "identity") +
    geom_text(data=data.frame(), aes(x=2200,y=38000,label="Colourless"), size=5, fontface=2) +
    xlab("Wavenumber (cm<sup>-1</sup>)") +
    ylab("Intensity") +
    scale_color_manual(name = "Colourless",
                       values = "#56B4E9") +
    theme_classic() +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 12, face = "bold", colour = "black"),
          legend.title = element_blank())

#plot the dark spectra
p.d <-
  ggplot(r.dark.long, aes(x = as.numeric(as.character(Wavenumber)))) +
    geom_line(aes(y = Intensity, colour = "", group = sample_id), linewidth = 0.5, stat = "identity") +
    geom_text(data=data.frame(), aes(x=2400,y=42000,label="Dark"), size=5, fontface=2) +
    xlab("Wavenumber (cm<sup>-1</sup>)") +
    ylab("Intensity") +
    scale_color_manual(name = "Dark",
                       values = "#36454F") +
    theme_classic() +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.title = element_blank())

#plot the light spectra
p.l <-
  ggplot(r.light.long, aes(x = as.numeric(as.character(Wavenumber)))) +
    geom_line(aes(y = Intensity, colour = "", group = sample_id), linewidth = 0.5, stat = "identity") +
    geom_text(data=data.frame(), aes(x=2400,y=39500,label="Light"), size=5, fontface=2) +
    xlab("Wavenumber (cm<sup>-1</sup>)") +
    ylab("Intensity") +
    scale_color_manual(name = "Light",
                       values = "#E69F00") +
    theme_classic() +
    theme(legend.position = "none",
          axis.title.x = element_markdown(size = 12, face = "bold", colour = "black"),
          axis.title.y = element_text(size = 12, face = "bold", colour = "black"),
          legend.title = element_blank())

#plot the white spectra
p.w <-
  ggplot(r.white.long, aes(x = as.numeric(as.character(Wavenumber)))) +
    geom_line(aes(y = Intensity, colour = "", group = sample_id), linewidth = 1, stat = "identity") +
    geom_text(data=data.frame(), aes(x=2400,y=18000,label="White"), size=5, fontface=2) +
    xlab("Wavenumber (cm<sup>-1</sup>)") +
    ylab("Intensity") +
    scale_color_manual(name = "White",
                       values = "#CC79A7") +
    theme_classic() +
    theme(legend.position = "none",
          axis.title.x = element_markdown(size = 12, face = "bold", colour = "black"),
          axis.title.y = element_blank(),
          legend.title = element_blank())

#Layout the plots in one figure
fig <-
  ggpubr::ggarrange(p.c, p.d, p.l, p.w,
                    ncol = 2,
                    nrow = 2)

#Save figure
ggsave("005-raman-spectra-summary.jpeg",
       fig,
       device = "jpeg",
       here::here("analysis", "figures"),
       width=20,
       height=20,
       units = "cm",
       dpi = 300)
