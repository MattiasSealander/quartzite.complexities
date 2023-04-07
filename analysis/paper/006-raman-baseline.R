suppressPackageStartupMessages(library(ChemoSpec))
suppressPackageStartupMessages(library(tidyverse))

#Import raman data, set empty fields to NA
raman.csv <-
  read.csv2(here::here("analysis", "data", "raw_data", "raman_raw_data.csv"), sep = ";", dec = ",", header = TRUE, check.names = FALSE, na = c("","NA","NULL",NULL))

#aggregate observations by group(sample) and calculate average of wavelength measurements
raman.averaged <-
  aggregate(raman.csv[, 4:468], list(sample_id = raman.csv$sample_id), mean)

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

write.csv(raman.tmp, here::here("analysis", "data", "derived_data", "raman_baseline_corr.csv"), row.names=FALSE)
