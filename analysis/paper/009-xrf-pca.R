suppressPackageStartupMessages(library(factoextra))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(robCompositions))
suppressPackageStartupMessages(library(tidyverse))

#### Load data ####

#Import xrf data, set empty fields to NA
xrf.csv <-
  read.csv2(here::here("analysis", "data", "raw_data", "xrf_quantitative_raw_data.csv"), sep = ";", dec = ".", header = TRUE, check.names = FALSE, na = c("","NA","NULL"))

#Import nir data, set empty fields to NA
nir.csv <-
  read.csv2(here::here("analysis", "data", "raw_data", "asd_raw_data.csv"), sep = ";", dec = ".", header = TRUE, check.names = FALSE, na = c("","NA","NULL",NULL))

#Import descriptive metadata
metadata.csv <-
  read.csv2(here::here("analysis", "data", "raw_data", "metadata.csv"), colClasses=c("sample_id"="character"), dec = ".", sep = ";", header = TRUE, na = c("", "NA", "NULL"), encoding = "UTF-8")

#merge XRF data with metadata
xrf.merged <-
  as.data.frame(merge(metadata.csv, xrf.csv, by='sample_id'))

#### Prepare dataset ####

#Filter xrf data to focus on points and preforms made from quartz/quartzite material
Points.xrf <-
  xrf.merged %>%
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

#filter XRF data on sample size (only include samples with smallest dimension >= 10mm),
#select elements to include in PCA, as well as columns to group by later on
xrf <-
  Points.xrf %>%
  filter(max_length_mm >= 10 & max_width_mm >= 10) %>%
  #percentage to ppm
  mutate_at(vars(Mg:Pb), ~ . * 10000) %>%
  dplyr::select(reading_no, site_id, sample_id, hue, material, `Al`, `Si`, `K`, `Ca`, `Fe`, `Zr`, `Ti`, `Ba`, `Mg`, `Sr`, `P`, `S`) %>%
  `row.names<-`(., NULL) %>%
  column_to_rownames(var = "reading_no")

#### Impute variables ####

#limit of detection elements Mg, Al, Si, P, S, Cl, K, Ca, V, Cr, Mn, Fe, Co, Cu, Zn, As, Rb, Sr, Y, Zr, Nb, Pd, Ag, Cd, Sn, Sb, Ba, Ti, Pb
#limit of detection values 2500, 487, 0, 47, 54, 38, 28, 15, 4, 4, 17, 14, 9, 5, 4, 2, 1, 2, 0, 2, 1, 2, 2, 2, 3, 4, 29, 9, 1

#Create vector with limit of detection for variables, this will work as a "censor" of highest possible value
dl <- c(487, 0, 28, 15, 14, 2, 9, 29, 2500, 2, 47, 54)

#use log-ratio Expectation - Maximisation algorithm to impute left-censored data (i.e values below LOD)
xrf <-
  cbind(xrf[,1:4],
        (zCompositions::lrEM(xrf[,5:16],
                             label=NA,
                             dl=dl,
                             rob=TRUE,
                             tolerance=0.005)
         #ppm to percentage
         /10000))

#### PCA full ####

#log-transform data in order to address issues with closure
#clr can also be used, in that case comment this log10 code, but gives similar results in this instance
xrf.log <-
  log10(xrf[,5:16])


#remove comment from code below if you want to transform data using centered log-ratio
#potential issues due to using geometric mean
#xrf.clr <-
#  compositions::clr(xrf[,5:16], ifclose = TRUE)

#perform pca with mean centering and unit variance scaling
xrf.pca <-
  prcomp(xrf.log, center = TRUE, scale.=TRUE)

#prepare labels for PCs
xrf.pc1var <- round(summary(xrf.pca)$importance[2,1]*100, digits=2)
xrf.pc2var <- round(summary(xrf.pca)$importance[2,2]*100, digits=2)
xrf.pc3var <- round(summary(xrf.pca)$importance[2,3]*100, digits=2)

xrf.pc1lab <- paste0("PC1 (",as.character(xrf.pc1var),"%)")
xrf.pc2lab <- paste0("PC2 (",as.character(xrf.pc2var),"%)")
xrf.pc3lab <- paste0("PC3 (",as.character(xrf.pc3var),"%)")


#prepare color/fill and symbols for score plots
pca.colors <- c("#56B4E9", "#36454F", "#E69F00", "#F9F6EE")
pca.hue <- c("Colourless", "Dark", "Light", "White")

#Extract loadings from xrf.pca
loadings <- as.data.frame(xrf.pca$rotation[,1:3])
#As elements are stored as rownames, set them as first column so it can be used as a variable
loadings <- setDT(loadings, keep.rownames = "Element")[]
#Melt loadings from wide to long format
loadings.melted <- melt(loadings, id.vars="Element")
#Set names of columns
loadings <- setNames(loadings.melted, c("Element", "PC", "Value"))

#Show eigenvalue scree plot
#fviz_eig(xrf.pca, choice = "eigenvalue")

#keep PCs above eigenvalue 1
xrf.transform = as.data.frame(-xrf.pca$x[,1:3])

#### K-means ####

#Determine no. of clusters with elbow graph
#fviz_nbclust(xrf.transform, kmeans, method = 'wss')

#Perform cluster analysis on the Principal Components
#In order to ensure that the cluster numbering is consistent between runs kmeans is calculated one time initially,
#center result is then ordered and passed to the kmeans in the second run. see: https://stackoverflow.com/questions/39906180/consistent-cluster-order-with-kmeans-in-r
centers <- kmeans(xrf.transform, centers = 3, nstart = 50)$centers

centers <- centers[order(-centers[,1], centers[,2]),]

kmeans.xrf <-
  kmeans(xrf.transform, centers = centers, nstart = 50)

#Prepare result of cluster analysis, so that they can be used in biplot
#In order to guard against potential mistakes, set rownames as column and use that for joining clusters to xrf dataframe
cluster <-
  as.data.frame(factor(kmeans.xrf$cluster)) %>%
  rename("cluster" = "factor(kmeans.xrf$cluster)") %>%
  rownames_to_column(var = "reading_no")

#Prepare legend with characteristic xrf variables for kmeans clusters to show in nir pca
cluster <-
  cluster %>%
  mutate(variables = case_when(
    endsWith(as.character(cluster), "1") ~ "Al, Ca, K, S, Si, Ti",
    endsWith(as.character(cluster), "2") ~ "Fe, Sr, Zr",
    endsWith(as.character(cluster), "3") ~ "Al, K"
  ))

#Set reading_no as column again join with kmeans cluster output
#Then group by kmeans cluster and calculate Si mean, concatenating the mean value with cluster for display in PCA legend
xrf <-
  xrf %>%
  rownames_to_column(var = "reading_no") %>%
  inner_join(cluster, by = "reading_no") %>%
  group_by(cluster) %>%
  mutate(Si_mean = round(mean(`Si`), 1)) %>%
  mutate(Si_mean = paste0("mean Si (", Si_mean, " %)")) %>%
  unite(legend, sep = " - ", cluster, Si_mean, remove = FALSE)

#pivot data and add to new dataframe, so a boxplot can be generated for each cluster and element
xrf_long <- xrf %>%
  dplyr::select(cluster, `Al`, `Si`, `K`, `Ca`, `Fe`, `Zr`, `Ti`, `Ba`, `Mg`, `Sr`, `P`, `S`) %>%
  pivot_longer(-cluster, names_to = "variable", values_to = "value")

#Prepare for facet labels
Elements <- c(
  'Al' = "Al",
  'Si' = "Si",
  'S' = "S",
  'K' = "K",
  'Ca' ="Ca",
  'Fe' = "Fe",
  'Zr' = "Zr",
  'Ti' = "Ti",
  'Ba' = "Ba",
  'Mg' = "Mg",
  'Sr' = "Sr",
  'P' = "P"
)

#### PCA low ####

#Filter out the cluster 1 samples, that feature those samples with >90% Si, then
#run pca and k-means one more time
xrf.low <-
  xrf %>%
  filter(cluster == "2" | cluster == "3") %>%
  ungroup() %>%
  select(-cluster, -legend, -Si_mean, -variables) %>%
  `row.names<-`(., NULL) %>%
  column_to_rownames(var = "reading_no")

#log-transform data in order to address issues with closure
xrf.log.low <-
  log10(xrf.low[,5:16])

#clr gives similar results in this case, at least in terms of exploration of the dataset
#xrf.clr <-
#  compositions::clr(xrf.low[,5:16], ifclose = TRUE)

#perform pca with mean centering and unit variance scaling
xrf.pca.low <-
  prcomp(xrf.log.low, center = TRUE, scale.=TRUE)

#prepare labels for PCs
xrf.pc1var.low <- round(summary(xrf.pca.low)$importance[2,1]*100, digits=2)
xrf.pc2var.low <- round(summary(xrf.pca.low)$importance[2,2]*100, digits=2)
xrf.pc3var.low <- round(summary(xrf.pca.low)$importance[2,3]*100, digits=2)
xrf.pc4var.low <- round(summary(xrf.pca.low)$importance[2,4]*100, digits=2)

xrf.pc1lab.low <- paste0("PC1 (",as.character(xrf.pc1var.low),"%)")
xrf.pc2lab.low <- paste0("PC2 (",as.character(xrf.pc2var.low),"%)")
xrf.pc3lab.low <- paste0("PC3 (",as.character(xrf.pc3var.low),"%)")
xrf.pc4lab.low <- paste0("PC4 (",as.character(xrf.pc4var.low),"%)")

#Extract loadings from xrf.pca
loadings.low <-
  as.data.frame(xrf.pca.low$rotation[,1:3])
#As elements are stored as rownames, set them as first column so it can be used as a variable
loadings.low <-
  setDT(loadings.low, keep.rownames = "Element")[]
#Melt loadings from wide to long format
loadings.melted.low <-
  melt(loadings.low, id.vars="Element")
#Set names of columns
loadings.low <-
  setNames(loadings.melted.low, c("Element", "PC", "Value"))

#Show eigenvalue scree plot
#fviz_eig(xrf.pca.low, choice = "eigenvalue")

#keep the 3 first PCs
xrf.transform.low = as.data.frame(-xrf.pca.low$x[,1:4])

#Determine no. of clusters with elbow graph
#fviz_nbclust(xrf.transform, kmeans, method = 'wss')

#Perform cluster analysis on the Principal Components
#In order to ensure that the cluster numbering is consistent between different script runs, k-means is calculated one time initially,
#the k-means center result is then ordered and passed to the k-means function in the second run.
#See: https://stackoverflow.com/questions/39906180/consistent-cluster-order-with-kmeans-in-r
centers.low <-
  kmeans(xrf.transform.low, centers = 3, nstart = 50)$centers

centers.low <-
  centers.low[order(-centers.low[,1], centers.low[,2]),]

kmeans.xrf.low <-
  kmeans(xrf.transform.low, centers = centers.low, nstart = 50)

#Prepare result of cluster analysis, so that they can be used in biplot
#In order to guard against potential mistakes, set rownames as column and use that for joining clusters to xrf dataframe
cluster.low <-
  as.data.frame(factor(kmeans.xrf.low$cluster)) %>%
  rename("cluster" = "factor(kmeans.xrf.low$cluster)") %>%
  rownames_to_column(var = "reading_no")

#Prepare legend with characteristic xrf variables for kmeans clusters to show in nir pca
cluster.low <-
  cluster.low %>%
  mutate(variables = case_when(
    endsWith(as.character(cluster), "1") ~ "Al, Ca, K, Mg, S, Sr, Ti",
    endsWith(as.character(cluster), "2") ~ "Fe, Ti, Zr",
    endsWith(as.character(cluster), "3") ~ "Si"
  ))

#Join the cluster result to xrf dataframe
xrf.low <-
  xrf.low %>%
  rownames_to_column(var = "reading_no") %>%
  inner_join(cluster.low, by = "reading_no") %>%
  group_by(cluster) %>%
  #add column storing mean Si value for each cluster for visualisation purposes
  mutate(Si_mean = round(mean(`Si`), 1)) %>%
  mutate(Si_mean = paste0("mean Si (", Si_mean, " %)")) %>%
  unite(legend, sep = " - ", cluster, Si_mean, remove = FALSE)

#pivot data and add to new dataframe, so a boxplot can be generated for each cluster and element
xrf_long_low <- xrf.low %>%
  dplyr::select(cluster, `Al`, `Si`, `K`, `Ca`, `Fe`, `Zr`, `Ti`, `Ba`, `Mg`, `Sr`, `P`, `S`) %>%
  pivot_longer(-cluster, names_to = "variable", values_to = "value")

#Prepare for facet labels
Elements.low <- c(
  'Al' = "Al",
  'Si' = "Si",
  'S' = "S",
  'K' = "K",
  'Ca' ="Ca",
  'Fe' = "Fe",
  'Zr' = "Zr",
  'Ti' = "Ti",
  'Ba' = "Ba",
  'Mg' = "Mg",
  'Sr' = "Sr",
  'P' = "P"
)


#### Figures ####

#Determine no. PC scree plot
load.1a <-
  fviz_eig(xrf.pca,
           choice = "eigenvalue",
           barfill = "#999999",
           barcolor = "black",
           linecolor = "black",
           ncp = 10,
           bar_width=0.5) +
  theme_classic() +
  theme(plot.title = element_blank(),
        axis.title.x = element_text(size = 10, face = "bold", colour = "black"),
        axis.title.y = element_text(size = 10, face = "bold", colour = "black"),
        legend.title = element_text(size = 10, face = "bold", colour = "black"),
        legend.background = element_rect(linetype = "solid", color = "black"),
        legend.position = "right")

#Determine no. PC scree plot
load.1b <-
  fviz_eig(xrf.pca.low,
           choice = "eigenvalue",
           barfill = "#999999",
           barcolor = "black",
           linecolor = "black",
           ncp = 10,
           bar_width=0.5) +
  theme_classic() +
  theme(plot.title = element_blank(),
        axis.title.x = element_text(size = 10, face = "bold", colour = "black"),
        axis.title.y = element_text(size = 10, face = "bold", colour = "black"),
        legend.title = element_text(size = 10, face = "bold", colour = "black"),
        legend.background = element_rect(linetype = "solid", color = "black"),
        legend.position = "right")


#bar graph of the loadings (P1-P2) for the XRF data
load.2a <-
  loadings %>%
  ggplot(aes(x=Element, y=Value)) +
  geom_bar(
    stat = "identity", position = "identity",
    color = "black", fill = "#999999",
    width = 0.5) +
  facet_wrap( ~ PC, scales = "free") +
  theme_bw() +
  theme(plot.margin = margin(10,10,10,10),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=10, face="bold", colour = "black"),
        axis.text.x = element_text(size=10, face="bold", colour = "black"),
        axis.title.x = element_blank(),
        strip.text.x = element_text(size=14, face="bold", colour="black"))

#bar graph of the loadings (P1-P2) for the XRF data
load.2b <-
  loadings.low %>%
  ggplot(aes(x=Element, y=Value)) +
  geom_bar(
    stat = "identity", position = "identity",
    color = "black", fill = "#999999",
    width = 0.5) +
  facet_wrap( ~ PC, scales = "free") +
  theme_bw() +
  theme(plot.margin = margin(10,10,10,10),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=10, face="bold", colour = "black"),
        axis.text.x = element_text(size=10, face="bold", colour = "black"),
        axis.title.x = element_blank(),
        strip.text.x = element_text(size=14, face="bold", colour="black"))


#Layout the scree and loading plots in one figure
fig1 <-
  ggpubr::ggarrange(ggpubr::ggarrange(load.1a,
                                      load.2a,
                                      ncol = 1,
                                      nrow = 2,
                                      labels = "A"),
                    ggpubr::ggarrange(load.1b,
                                      load.2b,
                                      ncol = 1,
                                      nrow = 2,
                                      labels = "B"),
                    ncol = 1,
                    nrow = 2)

#save fig3b object to render it in manuscript
saveRDS(fig1, here::here("analysis", "data", "derived_data", "009-xrf-loadings.rds"))

ggsave("009-xrf-pca-load-scree.jpeg",
       fig1,
       device = "jpeg",
       here::here("analysis", "figures"),
       scale = 1,
       width=25,
       height=20,
       units = "cm",
       dpi = 300)

#Biplot for the full dataset, with ellipses showing the results of the k-means cluster analysis
fig2a <-
  fviz_pca_biplot(xrf.pca,
                  axes = c(1,2),
                  geom = "text",
                  geom.var = c("point", "text"),
                  label = "var",
                  labelsize = 5,
                  col.ind = xrf$hue,
                  col.var = "black",
                  fill.var = "red",
                  repel = TRUE,
                  pointshape = 21,
                  pointsize = 2) +
  theme_bw() +
  labs(title = "XRF PCA Biplot",
       x = xrf.pc1lab,
       y = xrf.pc2lab) +
  #to add sample id as text
  #geom_text(aes(label=xrf$sample_id, hjust=0.5,vjust=-1.0)) +
  geom_point(aes(fill = xrf$hue,
                 shape = factor(xrf$material)),
             size = 2) +
  scale_shape_manual(name = "Material",
                     values=c(25,24,21),
                     guide = guide_legend(override.aes = list(fill = "black"),
                                          title.position="top",
                                          title.hjust = 0.5,
                                          order = 2)) +
  scale_color_manual(values = pca.colors,
                     guide = "none") +
  scale_fill_manual(name = "Hue",
                    values = pca.colors,
                    guide = guide_legend(override.aes = list(shape = 21,
                                                             fill = pca.colors,
                                                             color = "black",
                                                             size=3),
                                         title.position="top",
                                         title.hjust = 0.5,
                                         order = 1)) +
  ggnewscale::new_scale_fill() +
  stat_chull(aes(fill = xrf$legend),
             alpha = 0.3,
             geom = "polygon") +
  scale_fill_manual(name = "Cluster",
                    values = c("#0072B2", "#CC79A7","#F0E442"),
                    guide_legend(order = 1)) +
  theme(plot.title = element_text(size = 12, face = "bold", colour = "black"),
        axis.title.x = element_text(size = 12, face = "bold", colour = "black"),
        axis.title.y = element_text(size = 12, face = "bold", colour = "black"),
        legend.title = element_text(size = 12, face = "bold", colour = "black"),
        legend.background = element_rect(linetype = "solid", color = "black"),
        legend.position = "right")

#Biplot for low Si samples, with ellipses showing the results of the k-means cluster analysis
fig2b <-
  fviz_pca_biplot(xrf.pca.low,
                axes = c(1,2),
                geom = "text",
                geom.var = c("point", "text"),
                label = "var",
                labelsize = 5,
                col.ind = xrf.low$hue,
                col.var = "black",
                fill.var = "red",
                repel = TRUE,
                pointshape = 21,
                pointsize = 2) +
  theme_bw() +
  labs(title = "XRF PCA Biplot - Low Si ",
       x = xrf.pc1lab.low,
       y = xrf.pc2lab.low) +
  #to add sample id as text
  #geom_text(aes(label=xrf.low$sample_id, hjust=0.5,vjust=-1.0)) +
  geom_point(aes(fill = xrf.low$hue,
                 shape = factor(xrf.low$material)),
             size = 2) +
  scale_shape_manual(name = "Material",
                     values=c(25,24,21),
                     guide = "none") +
  scale_color_manual(values = pca.colors,
                     guide = "none") +
  scale_fill_manual(name = "Hue",
                    values = pca.colors,
                    guide = "none") +
  ggnewscale::new_scale_fill() +
  stat_chull(aes(fill = xrf.low$legend),
             alpha = 0.3,
             geom = "polygon") +
  scale_fill_manual(name = "Cluster",
                    values = c("#0072B2", "#CC79A7","#F0E442")) +
  theme(plot.title = element_text(size = 12, face = "bold", colour = "black"),
        axis.title.x = element_text(size = 12, face = "bold", colour = "black"),
        axis.title.y = element_text(size = 12, face = "bold", colour = "black"),
        legend.title = element_text(size = 12, face = "bold", colour = "black"),
        legend.background = element_rect(linetype = "solid", color = "black"),
        legend.position = "right")

#Layout the biplots in one figure
fig2 <-
  ggpubr::ggarrange(fig2a, fig2b,
                    ncol = 1,
                    nrow = 2,
                    labels = c("A", "B"),
                    font.label = list(size = 14, face = "bold"),
                    legend = "right")

#Save score plots for PC1-PC2
ggsave("009-xrf-pca.jpeg",
       fig2,
       device = "jpeg",
       here::here("analysis", "figures"),
       scale = 1,
       width=25,
       height=25,
       units = "cm",
       dpi = 300)
dev.off()

#Boxplot of XRF raw data by element and K-means cluster for the whole dataset
fig3a <-
  ggplot(xrf_long, aes(x=variable, y=value, fill = cluster)) +
  geom_boxplot() +
  facet_wrap( ~ variable, scales = "free", labeller = labeller(Element = Elements)) +
  theme_bw() +
  scale_fill_manual(name = "Cluster",
                    values = c("#0072B2", "#CC79A7","#F0E442")) +
  theme(
    legend.title = element_text(size = 12, face = "bold", colour = "black"),
    strip.text.x = element_text(size = 12, face = "bold", colour = "black"),
    axis.title.x = element_blank(),
    plot.title = element_text(size=10, face = "bold", colour = "black", margin = margin(t = 10, b = -20), vjust=0.02, hjust=0.01),
    axis.title.y=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank())


#Boxplot of XRF raw data by element and K-means cluster, with samples featuring >90% Si removed
#Due to the abundance of figures, not using this plot
fig3b <-
  ggplot(xrf_long_low, aes(x=variable, y=value, fill = cluster)) +
  geom_boxplot() +
  facet_wrap( ~ variable, scales = "free", labeller = labeller(Element = Elements.low)) +
  theme_bw() +
    scale_fill_manual(name = "Cluster",
                      values =c("#0072B2", "#CC79A7","#F0E442")) +
  theme(
    legend.title = element_text(size = 12, face = "bold", colour = "black"),
    strip.text.x = element_text(size = 12, face = "bold", colour = "black"),
    axis.title.x = element_blank(),
    plot.title = element_text(size=10, face = "bold", colour = "black", margin = margin(t = 10, b = -20), vjust=0.02, hjust=0.01),
    axis.title.y=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank())

#save fig3b object to render it in manuscript
saveRDS(fig3b, here::here("analysis", "data", "derived_data", "009-kmeans-box-modelb.rds"))

#Save K-means boxplot
ggsave("009-xrf-kmeans-boxplot.jpeg",
       fig3b,
       device = "jpeg",
       here::here("analysis", "figures"),
       scale = 1,
       width=20,
       height=15,
       units = "cm",
       dpi = 300)

#### NIR PCA ####
#Add XRF k-means cluster to NIR data to visualise in NIR PCA

#merge NIR data with metadata
nir.merged <-
  as.data.frame(merge(metadata.csv, nir.csv, by='sample_id'))

#identify the high Si samples and store as a matrix
high_si <-
  as.matrix(xrf %>%
  filter(cluster == "1") %>%
  ungroup() %>%
  select(sample_id))

#Filter NIR data to focus on points and preforms made from quartz/quartzite material, exclude observations on dark material with mainly noise
#then summarise measurements by group (sample)
Points.nir <-
  nir.merged %>%
  filter(type == "Point" | type == "Point fragment" | type == "Preform",
         material == "Brecciated quartz" | material == "Quartz" | material == "Quartzite") %>%
  dplyr::filter(!sample_id %in% c("153","167","168","169","172","174","175","177","182","183","190","191","193","194","196","198","200","204","207","210","213","214",
                                  "215","216","229","234","235","237","238","251","262","265","268","269","272","278","281","282","359","377","385","392","393","397","405",
                                  "406","410","411","413","414","415","416","417","424","425","426","428","430","432","55","56")) %>%
  #use previously stored high Si matrix as filter to select only the low Si samples
  dplyr::filter(!sample_id %in% high_si) %>%
  replace_na(list(munsell_hue = "Colourless")) %>%
  group_by(across(sample_id:weight_g)) %>%
  dplyr::summarise(across(`350.0`:`2500.0`, mean), .groups = "drop") %>%
  left_join(xrf.low[,c("sample_id", "variables")], by = "sample_id")

#perform PCA with mean-centering
nir.pca <-
  prcomp(Points.nir[,c(681:2180)], center = TRUE, scale = FALSE)

#Prepare axis labels with variance in %
nir.pc1lab <- as.data.frame(paste0("PC1 (",as.character(round(summary(nir.pca)$importance[2,1]*100, digits=2)),"%)"))
nir.pc2lab <- as.data.frame(paste0("PC2 (",as.character(round(summary(nir.pca)$importance[2,2]*100, digits=2)),"%)"))
nir.pc3lab <- as.data.frame(paste0("PC3 (",as.character(round(summary(nir.pca)$importance[2,3]*100, digits=2)),"%)"))

#prepare a basic score plot with fviz_pca_ind using PC1 and PC2
basic_plot1 <-
  fviz_pca_ind(nir.pca, axes = c(1,2), label="none")
#bind the basic plot with the columns storing legend items
nir <-
  cbind(basic_plot1$data, Points.nir[, c(2181,12,1)])

#PCA score plot with k-means cluster as shape fill
fig4 <-
  nir %>%
  ggplot(aes(x=x, y=y, shape = material, fill = variables)) +
  geom_text(data=subset(nir, sample_id == 386 | sample_id == 404 | sample_id == 391),
            aes(x,y,label=sample_id, hjust=1.5,vjust=-0.5)) +
  geom_point(size=4) +
  theme_bw() +
  ggtitle("Near infrared 1 000 - 2 500 nm") +
  labs(x = nir.pc1lab,
       y = nir.pc2lab) +
  scale_fill_manual(name = "XRF Cluster",
                    values=c("#0072B2", "#CC79A7","#F0E442"),
                    na.value="#D62728FF") +
  scale_shape_manual(name = "Material",
                     values=c(25,24,21)) +
  guides(fill = guide_legend(override.aes = list(shape = 22,
                                                 fill = c("#0072B2", "#CC79A7","#F0E442"), na.value="#D62728FF",
                                                 ncol = 4,
                                                 color = "black",
                                                 size=3),
                             title.position="top",
                             title.hjust = 0.5,
                             order = 1),
         shape = guide_legend(override.aes = list(fill = "black"),
                              ncol = 3,
                              title.position="top",
                              title.hjust = 0.5,
                              order = 2)) +
  ggnewscale::new_scale_color() +
  geom_point(data=subset(nir, is.na(variables)),
             aes(x=x, y=y, col="#D62728FF")) +
  scale_color_manual(name="Missing XRF data",
                     labels="NA",
                     values="#D62728FF") +
  guides(color = guide_legend(override.aes = list(shape = 22,
                                                 fill = "#D62728FF",
                                                 ncol = 1,
                                                 color = "black",
                                                 size=3,
                                                 order = 2),
                             title.position="top",
                             title.hjust = 0.5)) +
  theme(plot.title = element_text(size = 12, face = "bold", colour = "black", vjust = - 10, hjust = 0.02),
        axis.title.x = element_text(size = 12, face = "bold", colour = "black"),
        axis.title.y = element_text(size = 12, face = "bold", colour = "black"),
        legend.title = element_text(size = 12, face = "bold", colour = "black"),
        legend.background = element_rect(linetype = "solid", color = "black"),
        legend.position = "bottom")

#save fig4 object to render it in manuscript
saveRDS(fig4, here::here("analysis", "data", "derived_data", "009-kmeans-nir-pca.rds"))

#Save NIR plot with kmeans cluster
ggsave("009-nir-pca-xrf-kmeans.jpeg",
       fig4,
       device = "jpeg",
       here::here("analysis", "figures"),
       scale = 1,
       width=25,
       height=20,
       units = "cm",
       dpi = 300)
