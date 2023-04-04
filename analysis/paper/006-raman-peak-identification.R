suppressPackageStartupMessages(library(ChemoSpec))
suppressPackageStartupMessages(library(spectrolab))
suppressPackageStartupMessages(library(tidyverse))

#Import descriptive metadata
metadata.csv <-
  read.csv2(here::here("analysis", "data", "raw_data", "metadata.csv"), sep = ";", header = TRUE, na = c("", "NA", "NULL"), encoding = "UTF-8")

#Import raman data, set empty fields to NA
raman.csv <-
  read.csv2(here::here("analysis", "data", "raw_data", "raman_raw_data.csv"), sep = ";", dec = ",", header = TRUE, check.names = FALSE, na = c("","NA","NULL",NULL))

#aggregate observations by group(sample) and calculate average of wavelength measurements
raman.averaged <- 
  aggregate(raman.csv[, 5:469], list(sample_id = raman.csv$sample_id), mean)

#select sample_id and data columns and transpose in prep for baseline de-trend
raman.transposed <- 
  raman.averaged %>% 
  select(sample_id, `92.88`:`2503.59`)

raman.transposed <- 
  t(raman.transposed)

#Write csv with transposed raman data for baseline de-trending in ChemoSpec package
#ChemoSpec requires data to be read from file
write.table(raman.transposed, file = here::here("analysis", "data", "derived_data", "raman_transposed.csv"), col.names = FALSE, sep = ";", dec = )

#prepare vector with sample_id for reading raman data into spectra object
sample_id <- 
  raman.averaged$sample_id

#read transposed raman data into spectra object
raman <- suppressWarnings(matrix2SpectraObject(
  gr.crit = sample_id,
  gr.cols = c("auto"),
  freq.unit ="Wavelength cm-1",
  int.unit ="Intensity",
  descrip ="Bifacial points measurements",
  in.file = here::here("analysis", "data", "derived_data", "raman_transposed.csv"),
  sep = ";",
  dec = ".",
  chk = TRUE,
  out.file = here::here("analysis", "data", "derived_data", "raman_spec_object")))

#call on relevant baseline method (modified polynomial fitting) from baseline package and 
#return corrected spectra to spectra object
raman.baseline <- 
  baselineSpectra(raman,
                  int = FALSE,
                  method = "modpolyfit",
                  retC = TRUE)

#prevent baseline from plotting a figure
dev.off()

#spectra are stored in "Spectra" object, not a data frame
#need to extract samples and frequencies from that object 
raman.tmp <- 
  cbind(raman.baseline$freq, t(raman.baseline$data))
#transpose data to get frequencies as columns
raman.tmp <- 
  as.data.frame(t(raman.tmp))
#Add sample_id column and remove the initial "X" character introduced to sample names by R
raman.tmp$sample_id <- 
  c("sample_id", sub("^X*", "", raman.baseline$names))

#Make first row headers and remove the first row with the names
names(raman.tmp) <-
  raman.tmp[1,]
raman.tmp <- 
  raman.tmp[-1,]

#merge raman de-trended data with metadata
raman.baseline.merged <- 
  as.data.frame(merge(metadata.csv, raman.tmp, by='sample_id'))

#select and filter raman data to focus on relevant material and artefact types
Points.raman <-
  raman.baseline.merged %>%
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
  replace_na(list(munsell_hue = "Colourless"))

#Filter raman data to focus on dark samples 
r.dark <- Points.raman %>%
  filter(hue == "Dark") %>% 
  select(`92.88`:`2503.59`)

#Calculate the mean of all dark spectra
r.dark.mean <-
  as.data.frame(colMeans(r.dark))

#Filter raman data to focus on light samples 
r.light <- Points.raman %>%
  filter(hue == "Light") %>% 
  select(`92.88`:`2503.59`)

#Calculate the mean of all dark spectra
r.light.mean <-
  as.data.frame(colMeans(r.light))

#plot the dark spectra
p.d <-
  ggplot(r.dark.mean, aes(x = as.numeric(rownames(r.dark.mean)))) + 
  geom_line(aes(y = r.dark.mean[,1], colour = ""), linewidth = 1, stat = "identity") +
  geom_text(data=data.frame(), aes(x = 2400,y = 18000, label = "Dark"), size=6, fontface=2) +
  geom_text(data=data.frame(), aes(x = 120, y = 2800, label = "Quartz"), size=3, fontface=2) +
  geom_text(data=data.frame(), aes(x = c(200, 460),y = c(5000, 10000), label = "Quartz"), size=5, fontface=2) +
  geom_text(data=data.frame(), aes(x = c(1285, 1600),y = c(7500, 5500), label = c("D", "G")), size=5, fontface=2) +
  xlab("Wavenumber (cm-1)") +
  ylab("Intensity") +
  scale_color_manual(name = "Mean dark",
                     values = "black") + 
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank())

#plot the light spectra
p.l <-
  ggplot(r.light.mean, aes(x = as.numeric(rownames(r.light.mean)))) + 
  geom_line(aes(y = r.light.mean[,1], colour = ""), linewidth = 1, stat = "identity") +
  geom_text(data=data.frame(), aes(x=2400,y=18000,label="Light"), size=6, fontface=2) +
  xlab("Wavenumber (cm-1)") +
  ylab("Intensity") +
  scale_color_manual(name = "Mean light",
                     values = "black") + 
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank())


#Layout the plots in one figure
fig <-
  ggpubr::ggarrange(p.d, p.l,
                    nrow = 2)

#Save figure
ggsave("006-raman-peak-identification.jpeg",
       fig,
       device = "jpeg",
       here::here("analysis", "figures"),
       width=20, 
       height=15,
       units = "cm",
       dpi = 300)