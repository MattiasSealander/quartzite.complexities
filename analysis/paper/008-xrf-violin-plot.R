library(ggplot2)
suppressPackageStartupMessages(library(tidyverse))

#Import xrf data, set empty fields to NA
xrf.csv <-
  read.csv2(here::here("analysis", "data", "raw_data", "xrf_quantitative_raw_data.csv"), sep = ";", dec = ".", header = TRUE, check.names = FALSE, na = c("","NA","NULL"))

#Import descriptive metadata
metadata.csv <-
  read.csv2(here::here("analysis", "data", "raw_data", "metadata.csv"), colClasses=c("sample_id"="character"), dec = ".", sep = ";", header = TRUE, na = c("", "NA", "NULL"), encoding = "UTF-8")

#merge XRF data with metadata
xrf.merged <-
  as.data.frame(merge(metadata.csv, xrf.csv, by='sample_id'))

#Filter xrf data to focus on points and preforms made from quartz/quartzite material
Points.xrf <-
  xrf.merged %>%
  filter(type == "Point" | type == "Point fragment" | type == "Preform",
         material == "Brecciated quartz" | material == "Quartzite") %>%
  replace_na(list(munsell_hue = "Colourless"))

#Prepare XRF data in long format for violin plots, exclude samples with a dimension smaller than 10 mm
Points_long <-
  drop_na(Points.xrf %>%
            filter(max_length_mm >= 10 | max_width_mm >= 10) %>%
            #leave out Zn and Y as these have less than 30 values above LOD
            dplyr::select(`Mg`,`Al`, `Si`, `P`, `S`, `K`, `Ca`, `V`, `Mn`, `Fe`, `Sr`, `Zr`, `Ba`, `Ti`) %>%
            gather(key = "Element", value = "Value"))

# Add column with elemental groups for visualisation
Points_long <-
  Points_long %>%
  mutate(Group = case_when(
    endsWith(Element, "Mg") ~ "Lithophile",
    endsWith(Element, "Al") ~ "Lithophile",
    endsWith(Element, "Si") ~ "Lithophile",
    endsWith(Element, "P") ~ "Lithophile",
    endsWith(Element, "S") ~ "Chalcophile",
    endsWith(Element, "K") ~ "Lithophile",
    endsWith(Element, "Ca") ~ "Lithophile",
    endsWith(Element, "V") ~ "Lithophile",
    endsWith(Element, "Mn") ~ "Siderophile",
    endsWith(Element, "Fe") ~ "Siderophile",
    endsWith(Element, "Sr") ~ "Lithophile",
    endsWith(Element, "Zr") ~ "Lithophile",
    endsWith(Element, "Ba") ~ "Lithophile",
    endsWith(Element, "Ti") ~ "Chalcophile",
  ))

#Calculate and add value equal to 90% of the value range for each element, this is in preparation for
#text annotations in final plot where a text will be positioned at this values height in each plot
Points_long <-
  Points_long %>%
  group_by(Element) %>%
  dplyr::mutate(ypos = max(Value) * 0.95)

#Lookup for element labeller
Elements <- c(
  'Mg' = "Mg",
  'Al' = "Al",
  'Si' = "Si",
  'P' = "P",
  'S'= "S",
  'K' = "K",
  'Ca' ="Ca",
  'V' = "V",
  'Mn' = "Mn",
  'Fe' = "Fe",
  'Sr' = "Sr",
  'Zr' = "Zr",
  'Ba' = "Ba",
  'Ti' = "Ti"
)

#Prepare data frame with elements for annotating facets in final plot
na_values <-
  data.frame("Element" = unique(Points_long$Element))

#Calculate and add number of missing values for each element
na_values <-
  na_values %>%
  mutate(na_values = c(
    paste("n = ", length(which(!is.na(Points.xrf$`Mg`)))),
    paste("n = ", length(which(!is.na(Points.xrf$`Al`)))),
    paste("n = ", length(which(!is.na(Points.xrf$`Si`)))),
    paste("n = ", length(which(!is.na(Points.xrf$`P`)))),
    paste("n = ", length(which(!is.na(Points.xrf$`S`)))),
    paste("n = ", length(which(!is.na(Points.xrf$`K`)))),
    paste("n = ", length(which(!is.na(Points.xrf$`Ca`)))),
    paste("n = ", length(which(!is.na(Points.xrf$`V`)))),
    paste("n = ", length(which(!is.na(Points.xrf$`Mn`)))),
    paste("n = ", length(which(!is.na(Points.xrf$`Fe`)))),
    paste("n = ", length(which(!is.na(Points.xrf$`Sr`)))),
    paste("n = ", length(which(!is.na(Points.xrf$`Zr`)))),
    paste("n = ", length(which(!is.na(Points.xrf$`Ba`)))),
    paste("n = ", length(which(!is.na(Points.xrf$`Ti`)))))
  )

#Merge na_values and Points_long data frames in order to get the y position that the text annotation should be displayed
#in for each element facet in final plot
na_values <-
  merge(na_values, unique(Points_long[,c("Element", "ypos")]), by = "Element", all.x = TRUE)


#Prepare ylimits for each group in order to exclude extreme outliers in the violin plot graphic
ylims <- Points_long %>%
  group_by(Element) %>%
  mutate(height = max(Value) + .3 * sd(Value))

#Prepare a list with formula for setting new ylimits for each group in ggh4x package "facetted_pos_scales() which makes it possible to work with facet_wrap()
scales <- c(
  Element == 'Mg' ~ scale_y_continuous(limits = c(as.vector(t(ylims[1,2:3])))),
  Element == 'Al' ~ scale_y_continuous(limits = c(as.vector(t(ylims[2,2:3])))),
  Element == 'Si' ~ scale_y_continuous(limits = c(as.vector(t(ylims[3,2:3])))),
  Element == 'P' ~ scale_y_continuous(limits = c(as.vector(t(ylims[4,2:3])))),
  Element == 'S' ~ scale_y_continuous(limits = c(as.vector(t(ylims[5,2:3])))),
  Element == 'K' ~ scale_y_continuous(limits = c(as.vector(t(ylims[6,2:3])))),
  Element == 'Ca' ~ scale_y_continuous(limits = c(as.vector(t(ylims[7,2:3])))),
  Element == 'V' ~ scale_y_continuous(limits = c(as.vector(t(ylims[8,2:3])))),
  Element == 'Mn' ~ scale_y_continuous(limits = c(as.vector(t(ylims[9,2:3])))),
  Element == 'Fe' ~ scale_y_continuous(limits = c(as.vector(t(ylims[10,2:3])))),
  Element == 'Sr' ~ scale_y_continuous(limits = c(as.vector(t(ylims[12,2:3])))),
  Element == 'Zr' ~ scale_y_continuous(limits = c(as.vector(t(ylims[14,2:3])))),
  Element == 'Ba' ~ scale_y_continuous(limits = c(as.vector(t(ylims[15,2:3])))),
  Element == 'Ti' ~ scale_y_continuous(limits = c(as.vector(t(ylims[16,2:3]))))
)

#Violin plot of elemental content, trimmed, with fill based on elemental groups and number of values above LOD as text annotation
fig <-
  Points_long %>%
    ggplot(aes(x=Element, y=Value, fill = Group)) +
    see::geom_violinhalf(trim = TRUE) +
    geom_boxplot(width=.1, outlier.shape = NA, position= position_nudge(x=-.1)) +
    facet_wrap( ~ Element, scales = "free", labeller = labeller(Element = Elements)) +
    scale_fill_manual(name = "Group",
                      values = c("#E69F00", "#56B4E9", "#CC79A7")) +
    theme_bw() +
    geom_text(data=na_values, aes(x=0.28, y=ypos, label=na_values, fontface = "bold"), nudge_x = 0.4, inherit.aes = FALSE, size = 3) +
    theme(
      legend.title = element_text(size = 12, face = "bold", colour = "black"),
      strip.text.x = element_text(size = 12, face = "bold", colour = "black"),
      axis.title.x = element_blank(),
      plot.title = element_text(size=10, face = "bold", colour = "black", margin = margin(t = 10, b = -20), vjust=0.02, hjust=0.01),
      axis.title.y=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank())

ggsave("008-xrf-violin-plot.jpeg",
       fig,
       device = "jpeg",
       here::here("analysis", "figures"),
       width=20,
       height=20,
       units = "cm",
       dpi = 300)
