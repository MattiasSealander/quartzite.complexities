# Loading package
library(png)
library(loder)
suppressPackageStartupMessages(library(ChemoSpec))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(flextable))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(tidyverse))

#Import descriptive metadata
metadata.csv <-
  read.csv2(here::here("analysis", "data", "raw_data", "metadata.csv"), sep = ";", header = TRUE, na = c("", "NA", "NULL"), encoding = "UTF-8")

#Import nir data, set empty fields to NA
nir.csv <-
  read.csv2(here::here("analysis", "data", "raw_data", "asd_raw_data.csv"), sep = ";", dec = ".", header = TRUE, check.names = FALSE, na = c("","NA","NULL",NULL))

#Import raman data, set empty fields to NA
raman.csv <-
  read.csv2(here::here("analysis", "data", "raw_data", "raman_raw_data.csv"), sep = ";", dec = ",", header = TRUE, check.names = FALSE, na = c("","NA","NULL",NULL))

#Import xrf data, set empty fields to NA
xrf.csv <-
  read.csv2(here::here("analysis", "data", "raw_data", "xrf_quantitative_raw_data.csv"), sep = ";", dec = ".", header = TRUE, check.names = FALSE, na = c("","NA","NULL"))

#aggregate observations by group(sample) and calculate average of wavelength measurements
nir.averaged <- 
  aggregate(nir.csv[, 4:2154], list(sample_id = nir.csv$sample_id), mean)

#merge NIR data with metadata
nir.merged <- 
  as.data.frame(merge(metadata.csv, nir.averaged, by='sample_id'))

#merge XRF data with metadata
xrf.merged <- 
  as.data.frame(merge(metadata.csv, xrf.csv, by='sample_id'))

#aggregate observations by group(sample) and calculate average of wavelength measurements
raman.averaged <- 
  aggregate(raman.csv[, 5:469], list(sample_id = raman.csv$sample_id), mean)

#Filter nir data to focus on points and preforms made from quartz/quartzite material 
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
                material == "Brecciated quartz" | material == "Quartz" | material == "Quartzite") 

#select sample_id and data columns and transpose in prep for baseline de-trend
raman.transposed <- 
  raman.averaged %>% 
  dplyr::select(sample_id, `92.88`:`2503.59`)

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

#prevent baselineSpectra from plotting a figure which will show up in the knitted document
dev.off()

#spectra are stored in a "Spectra" object, not a data frame
#need to extract samples and frequencies from that object 
raman.tmp <- 
  #bind data from spectra object
  as.data.frame(cbind(raman.baseline$freq, t(raman.baseline$data))) %>% 
  #transpose data so that frequencies are column headers and observations rows
  tibble::rownames_to_column() %>% 
  pivot_longer(-rowname) %>% 
  pivot_wider(names_from=rowname, values_from=value) %>%
  #remove unnecessary column with previous column headers
  dplyr::select(-name) %>% 
  #Set frequencies as headers, as they are still first row
  purrr::set_names(as.character(slice(., 1))) %>%
  slice(-1) %>% 
  #add sample id from spectra object
  add_column(sample_id = sub("^X*", "", raman.baseline$names))

#merge raman de-trended data with metadata
raman.baseline.merged <- 
  as.data.frame(merge(metadata.csv, raman.tmp, by='sample_id'))

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
         material == "Brecciated quartz" | material == "Quartz" | material == "Quartzite")

#Filter xrf data to focus on points and preforms made from quartz/quartzite material 
Points.xrf <-
  xrf.merged %>%
  dplyr::select(sample_id, `Mg`, `Al`, `Si`, `P`, `S`, `K`, `Ca`, `Ti`, `Fe`, `Sr`, `Zr`,  `Ba`) %>% 
  dplyr::filter(sample_id == 386 | sample_id == 404| sample_id == 391) %>% 
  rename("Mg(%)" = `Mg`, "Al(%)" = `Al`, "Si(%)" = `Si`, "P(%)" = `P`, "S(%)" = `S`, "K(%)" = `K`, "Ca(%)" = `Ca`, "Ti(%)" = `Ti`, "Fe(%)" = `Fe`, 
         "Sr(%)" = `Sr`, "Zr(%)" = `Zr`,  "Ba(%)" = `Ba`)

ft <- 
  flextable(Points.xrf) %>% 
  set_header_labels(sample_id = "Sample") %>% 
  theme_vanilla() %>% 
  set_caption(caption = "XRF - Elemental content (%)") %>% 
  as_raster()

x.ft <- 
  ggplot() + 
  theme_void() + 
  annotation_custom(rasterGrob(ft), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)

#Filter NIR data to focus on colourless 
#and select the NIR range 1 000 - 2 500 nm
n.231 <- Points.nir %>%
  dplyr::filter(sample_id == 231) %>% 
  dplyr::select(sample_id, `1001.0`:`2500.0`)

#Melt into long format
n.231.long <- 
  suppressWarnings(melt(setDT(n.231), variable.name = "Wavelength", variable.factor = FALSE, value.name = "Absorbance"))

#Filter NIR data to focus on colourless 
#and select the NIR range 1 000 - 2 500 nm
n.404 <- Points.nir %>%
  dplyr::filter(sample_id == 249) %>% 
  dplyr::select(sample_id, `1001.0`:`2500.0`)

#Melt into long format
n.404.long <- 
  suppressWarnings(melt(setDT(n.404), variable.name = "Wavelength", variable.factor = FALSE, value.name = "Absorbance"))

#Filter NIR data to focus on colourless 
#and select the NIR range 1 000 - 2 500 nm
n.391 <- Points.nir %>%
  dplyr::filter(sample_id == 391) %>% 
  dplyr::select(sample_id, `1001.0`:`2500.0`)

#Melt into long format
n.391.long <- 
  suppressWarnings(melt(setDT(n.391), variable.name = "Wavelength", variable.factor = FALSE, value.name = "Absorbance"))

#plot the colourless spectra
pn.231 <-
  ggplot(n.231.long, aes(x = as.numeric(as.character(Wavelength)))) + 
  geom_line(aes(y = Absorbance, colour = ""), linewidth = 1, stat = "identity") +
  xlab("Wavelength (nm)") +
  ylab("Absorbance") +
  labs(title="386") +
  scale_color_manual(name = "386",
                     values = "black") + 
  scale_x_continuous(limits=c(1000,2500), breaks = scales::pretty_breaks(n = 5)) +
  theme_classic() +
  theme(legend.position = "none",
        plot.title = element_text(size=16, 
                                  hjust = .9),        
        axis.title.x = element_text(size = 10, colour = "black"),
        axis.title.y = element_text(size = 12, face = "bold", colour = "black"),
        legend.title = element_text(size = 12, face = "bold", colour = "black"))

#plot the colourless spectra
pn.404 <-
  ggplot(n.404.long, aes(x = as.numeric(as.character(Wavelength)))) + 
  geom_line(aes(y = Absorbance, colour = ""), linewidth = 1, stat = "identity") +
  xlab("Wavelength (nm)") +
  #ylab("Absorbance") +
  labs(title="404") +
  scale_color_manual(name = "404",
                     values = "black") +
  scale_x_continuous(limits=c(1000,2500), breaks = scales::pretty_breaks(n = 5)) +
  theme_classic() +
  theme(legend.position = "none",
        plot.title = element_text(size=16, 
                                  hjust = .9),
        axis.title.x = element_text(size = 10, colour = "black"),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 12, face = "bold", colour = "black"))


#plot the colourless spectra
pn.391 <-
  ggplot(n.391.long, aes(x = as.numeric(as.character(Wavelength)))) + 
  geom_line(aes(y = Absorbance, colour = ""), linewidth = 1, stat = "identity") +
  xlab("Wavelength (nm)") +
  #ylab("Absorbance") +
  labs(title="391") +
  scale_color_manual(name = "391",
                     values = "black") + 
  scale_x_continuous(limits=c(1000,2500), breaks = scales::pretty_breaks(n = 5)) +
  theme_classic() +
  theme(legend.position = "none",
        plot.title = element_text(size=16, 
                                  hjust = .9),
        axis.title.x = element_text(size = 10, colour = "black"),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 12, face = "bold", colour = "black"))

#Filter raman data for sample 231 
r.231 <- Points.raman %>%
  filter(sample_id == 231) %>% 
  select(sample_id, c(`92.88`:`2503.59`))

#Melt into long format
r.231.long <- 
  suppressWarnings(reshape2::melt(data.table::setDT(r.231), variable.name = "Wavenumber", variable.factor = FALSE, value.name = "Intensity"))

#Filter raman data for sample 249 
r.404 <- Points.raman %>%
  filter(sample_id == 404) %>% 
  select(sample_id, c(`92.88`:`2503.59`))

#Melt into long format
r.404.long <- 
  suppressWarnings(reshape2::melt(data.table::setDT(r.404), variable.name = "Wavenumber", variable.factor = FALSE, value.name = "Intensity"))

#Filter raman data for sample 391
r.391 <- Points.raman %>%
  filter(sample_id == 391) %>% 
  select(sample_id, c(`92.88`:`2503.59`))

#Melt into long format
r.391.long <- 
  suppressWarnings(reshape2::melt(data.table::setDT(r.391), variable.name = "Wavenumber", variable.factor = FALSE, value.name = "Intensity"))

#plot the 231 spectra
pr.231 <-
  ggplot(r.231.long, aes(x = as.numeric(as.character(Wavenumber)))) + 
  geom_line(aes(y = Intensity, colour = "", group = sample_id), linewidth = 1, stat = "identity") +
  xlab("Wavenumber (cm-1)") +
  ylab("Intensity") +
  labs(title="386") +
  scale_color_manual(name = "Dark",
                     values = "black") + 
  theme_classic() +
  theme(legend.position = "none",
        plot.title = element_text(size=16, 
                                  hjust = .9),
        axis.title.x = element_text(size = 10, colour = "black"),
        axis.title.y = element_text(size = 12, face = "bold", colour = "black"),
        legend.title = element_text(size = 12, face = "bold", colour = "black"))

#plot the 249 spectra
pr.404 <-
  ggplot(r.404.long, aes(x = as.numeric(as.character(Wavenumber)))) + 
  geom_line(aes(y = Intensity, colour = "", group = sample_id), linewidth = 1, stat = "identity") +
  xlab("Wavenumber (cm-1)") +
  ylab("Intensity") +
  labs(title="404") +
  scale_color_manual(name = "Dark",
                     values = "black") + 
  theme_classic() +
  theme(legend.position = "none",
        plot.title = element_text(size=16, 
                                  hjust = .9),
        axis.title.x = element_text(size = 10, colour = "black"),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 12, face = "bold", colour = "black"))

#plot the 391 spectra
pr.391 <-
  ggplot(r.391.long, aes(x = as.numeric(as.character(Wavenumber)))) + 
  geom_line(aes(y = Intensity, colour = "", group = sample_id), linewidth = 1, stat = "identity") +
  xlab("Wavenumber (cm-1)") +
  ylab("Intensity") +
  labs(title="391") +
  scale_color_manual(name = "Dark",
                     values = "black") + 
  theme_classic() +
  theme(legend.position = "none",
        plot.title = element_text(size=16, 
                                  hjust = .9),
        axis.title.x = element_text(size = 10, colour = "black"),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 12, face = "bold", colour = "black"))

img.231 <- readPNG(here::here("analysis", "figures", "386.png"),  native=TRUE)
g.231 <- rasterGrob(img.231, interpolate=TRUE)

img.404 <- readPNG(here::here("analysis", "figures", "404.png"),  native=TRUE)
g.404 <- rasterGrob(img.404, interpolate=TRUE)

img.391 <- readPNG(here::here("analysis", "figures", "391.png"),  native=TRUE)
g.391 <- rasterGrob(img.391, interpolate=TRUE)

top <- 
  plot_grid(x.ft)
p.nir <- 
  plot_grid(pn.231, pn.391, pn.404, labels = "NIRS", nrow=1)
p.raman <- 
  plot_grid(pr.231, pr.391, pr.404, labels = "Raman", nrow=1)
pics <- 
  plot_grid(g.231, g.391, g.404, nrow=1)


fig <- 
  cowplot::plot_grid(top,p.nir,p.raman,pics, 
                     nrow = 4, 
                     align="hv")

#Save figure
ggsave("010-nir-sample-comparison.jpeg",
       fig,
       device = "jpeg",
       here::here("analysis", "figures"),
       width=25, 
       height=25,
       units = "cm",
       dpi = 300)
