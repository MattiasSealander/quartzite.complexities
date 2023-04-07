# Loading package
library(png)
library(loder)
suppressPackageStartupMessages(library(ChemoSpec))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(flextable))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggtext))

#Import descriptive metadata
metadata.csv <-
  read.csv2(here::here("analysis", "data", "raw_data", "metadata.csv"), colClasses=c("sample_id"="character"), dec = ".", sep = ";", header = TRUE, na = c("", "NA", "NULL"), encoding = "UTF-8")

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
  aggregate(raman.csv[, 4:468], list(sample_id = raman.csv$sample_id), mean)

#Filter nir data to focus on points and preforms made from quartz/quartzite material
Points.nir <-
  nir.merged %>%
  dplyr::filter(type == "Point" | type == "Point fragment" | type == "Preform",
                material == "Brecciated quartz" | material == "Quartz" | material == "Quartzite")

raman.baseline <-
  read.csv2(here::here("analysis", "data", "derived_data", "raman_baseline_corr.csv"), sep = ",", dec = ".", header = TRUE, check.names = FALSE, na = c("", "NA", "NULL"), encoding = "UTF-8")

#merge raman de-trended data with metadata
raman.baseline.merged <-
  as.data.frame(merge(metadata.csv, raman.baseline, by='sample_id'))

Points.raman <-
  raman.baseline.merged %>%
  filter(type == "Point" | type == "Point fragment" | type == "Preform",
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
  xlab("Wavenumber (cm<sup>-1</sup>)") +
  ylab("Intensity") +
  labs(title="386") +
  scale_color_manual(name = "Dark",
                     values = "black") +
  theme_classic() +
  theme(legend.position = "none",
        plot.title = element_text(size=16,
                                  hjust = .9),
        axis.title.x = element_markdown(size = 10, colour = "black"),
        axis.title.y = element_text(size = 12, face = "bold", colour = "black"),
        legend.title = element_text(size = 12, face = "bold", colour = "black"))

#plot the 249 spectra
pr.404 <-
  ggplot(r.404.long, aes(x = as.numeric(as.character(Wavenumber)))) +
  geom_line(aes(y = Intensity, colour = "", group = sample_id), linewidth = 1, stat = "identity") +
  xlab("Wavenumber (cm<sup>-1</sup>)") +
  ylab("Intensity") +
  labs(title="404") +
  scale_color_manual(name = "Dark",
                     values = "black") +
  theme_classic() +
  theme(legend.position = "none",
        plot.title = element_text(size=16,
                                  hjust = .9),
        axis.title.x = element_markdown(size = 10, colour = "black"),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 12, face = "bold", colour = "black"))

#plot the 391 spectra
pr.391 <-
  ggplot(r.391.long, aes(x = as.numeric(as.character(Wavenumber)))) +
  geom_line(aes(y = Intensity, colour = "", group = sample_id), linewidth = 1, stat = "identity") +
  xlab("Wavenumber (cm<sup>-1</sup>)") +
  ylab("Intensity") +
  labs(title="391") +
  scale_color_manual(name = "Dark",
                     values = "black") +
  theme_classic() +
  theme(legend.position = "none",
        plot.title = element_text(size=16,
                                  hjust = .9),
        axis.title.x = element_markdown(size = 10, colour = "black"),
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
       width=3508,
       height=2481,
       units = "px",
       dpi = 300)
