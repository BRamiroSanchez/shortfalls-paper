
# This script provides the code for the bioregionalization using
# network analysis with the R package 'biogeonetworks' (Leroy, 2021)

# please check the package; the code below is an adaptation from
# the package tutorial: https://github.com/Farewe/biogeonetworks


###########################################################################################
#    NOTE: since some data cannot be shared for confidentiality, we provide directly the
#          network to reproduce the analysis (i.e., from STEP 3)
###########################################################################################


# load packages
list.packages <- c("sp", "raster", "rgdal", "ggplot2", "plyr", "sf", "rnaturalearth",
                   "dplyr", "RColorBrewer", "gridExtra", "stringr")
lapply(list.packages, require, character.only = TRUE)

#install.packages("devtools")
#devtools::install_github("Farewe/biogeonetworks")
library(biogeonetworks)


# Read in species data
taxa_list <- readRDS("./data/taxa_list.RDS") # curated taxa list of taxonomic names only
taxa_data <- droplevels(readRDS("./data/taxa_data.RDS")) #  obis + gbif + smithsonian + siofa + mnhn

taylor <- read.csv("./data/Taylor_and_Rogers_2017.csv") # literature data
taylor <-  taylor[,c(1,8,11,12,13,14,15)] # I'm selecting certain columns to merge with "taxa_data" dataset

species_df <- taxa_data[which(taxa_data$species %in% taxa_list$taxon), ]  # data only at species level (match with curated list)
species_df <- species_df[,c(29, 7,8,9, 16, 20, 15)] # select columns to merge with the "taylor" dataset
species_df <- rbind(species_df, taylor)

species_df <- SpatialPointsDataFrame(species_df[,c("decimalLongitude", "decimalLatitude")],
                                     data=species_df,
                                     proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))

# Read in environmental data
baseline <- readRDS("./data/baseline.RDS")  # environmental data at 1 degree spatial resolution to obtain cell id number (from Bio-Oracle)
                                                  # Alternatively, you can create a 1 degree spatial grid and generate cell id numbers
                                                  # for the same purpose.

# Aggregate records to 1D cells
species_df$cellid <- cellFromXY(baseline, species_df) # extract cell number

species_df <- as.data.frame(species_df) # transform into a dataframe to be able to use duplicated()
species_df <- species_df[-which(duplicated(species_df[, c("species", "cellid")])), ]  # remove duplicate records

vme.df <- species_df[, c("species", "cellid")] # basically, you only need this information for the analysis


# Read in spatial data to create plots later on, but not necessar for the network analysis
worldMap <- ne_countries(scale = "medium", type = "countries", returnclass = "sf")
iho <- read_sf("./data/spatial/iho/iho.shp") %>% st_as_sf() # Indian Ocean
limsSI <- st_buffer(iho, dist = 0.7) %>% st_bbox() # create some limites for plotting later on
siofa <- read_sf("./data/spatial/RFB_SIOFA.shp") %>% st_as_sf() # siofa polygon boundary shapefile
land <- readRDS("./data/spatial/ne_10m_land/land") # land shapefile



####################    NETWORK ANALYSIS    #################### 

# please check 'biogeonetworks' package (Leroy, 2021); the code below is an adaptation from
# the tutorial given for the package: https://github.com/Farewe/biogeonetworks

setwd("./data/")

# STEP 1: CREATE THE NETWORK
writePajek(vme.df, site.field=2, species.field=1, "./vmes.net")


# STEP 2: LAUNCH MAP EQUATION
system("./R/Infomap  vmes.net   ./  --flow-model undirected --tree -N1000") # where you have the executable for Infomap


# STEP 3: READ THE TREE
vmes.tree <- readInfomapTree("./vmes.tree",  
                             network.summary = TRUE, # Prints a summary of the clustering results
                             replace.leaf.names = TRUE)


# STEP 4: ANALYSE THE SITES AND SPECIES INDIVIDUALLY
sites <- getSiteTable(vme.df, site.field= 2, vmes.tree, drop.levels = TRUE)
sps <- getSpeciesTable(vme.df, species.field = 1, vmes.tree, drop.levels = TRUE)

# Inspect the number of occurrences per basin
count(sps, vars = "lvl1") #1 lvl1 1997 
count(sites, vars = "lvl1") #1 lvl1 613 

# Investigate "lvl1" and "lvl2"
plyr::count(sites$lvl1)
lvl1.list <- c("1", "2", "3", "4") # number of bioregions 

plyr::count(sites$lvl2)
lvl2.list  <- c("1.1",  "1.3", "1.2", "2.1", "1.5", "1.12", "3.1", "1.23", "2.3") # same order as freq table

plyr::count(sps$lvl1)

levels(vmes.tree$lvl1)
levels(vmes.tree$lvl2)


# STEP 5: COLOUR EACH REGION
vmes.tree <- attributeColors(network = vmes.tree, # Same object as input & output
                             lvl = "lvl1", # Which hierarchical level are we working on?
                             nb.max.colors = 4, # We chose 4 nb. clusters as significant for level 1
                             palette = "Set1",
                             colname = "colors.lvl1", # Name to give to the colour column
                             other.color = grey(0.5), # Minor clusters will be black
                             cluster.order = "sites", # Cluster order for colours (see below)
                             db = vme.df, # Database of step 1
                             site.field = 2, # Name of site column in your database
                             species.field = 1) # Name of species column in your database


vmes.tree <- attributeColors(network = vmes.tree, # Same object as input & output
                             lvl = "lvl2", # Which hierarchical level are we working on?
                             nb.max.colors = 9, # We chose 9 significant clusters for level 2 with n > 10 sites
                             palette = "Set3",
                             colname = "colors.lvl2", # Name to give to the colour column
                             other.color = grey(0.5), # Minor clusters will be black
                             cluster.order = "sites", # Cluster order for colours (see below)
                             db = vme.df, # Database of step 1
                             site.field = "cellid", # Name of site column in your database
                             species.field = "species") # Name of species column in your database


# STEP 6: MATCH THE SHAPEFILE WITH THE COLOURS OF THE REGIONS (here I plot the bioregion occurrences)

library(rgdal)
library(sf)

r <- baseline[[1]]
r2 <- setValues(r, 1:ncell(r))
names(r2) <- "cellid"
p <- rasterToPolygons(r2)
grid_poly <- st_as_sf(p) %>% st_cast("MULTIPOLYGON")

grid_poly$colors.lvl1 <- vmes.tree$colors.lvl1[match(grid_poly$cellid, vmes.tree$Name)]
grid_poly$colors.lvl2 <- vmes.tree$colors.lvl2[match(grid_poly$cellid, vmes.tree$Name)]
grid_poly$lvl2 <- vmes.tree$lvl2[match(grid_poly$cellid, vmes.tree$Name)]
grid_poly$lvl1 <- vmes.tree$lvl1[match(grid_poly$cellid, vmes.tree$Name)]

# make plot to visualise them
lvl1.clusters <- grid_poly[!is.na(grid_poly$colors.lvl1),]
plot(st_geometry(lvl1.clusters), col = lvl1.clusters$colors.lvl1)

lvl2.clusters <- grid_poly[!is.na(grid_poly$colors.lvl2),]
lvl2.clusters$lvl2 <- as.character(lvl2.clusters$lvl2)
lvl2.clusters$lvl2[lvl2.clusters$colors.lvl2 == "black"] <- "Other"


lvl1.clusters$colors.lvl1 <- factor(lvl1.clusters$colors.lvl1, levels=c("#E41A1C" ,"#377EB8" ,"#4DAF4A" , "#984EA3", "#808080"), labels= c("1", "2", "3", "4", "Other"))
lvl1.clusters.pt <- st_centroid(lvl1.clusters)


lvl1 <- ggplot(lvl1.clusters.pt) + 
  geom_sf(data = siofa, fill = NA, size = 0.3, color="black") +
  geom_sf(data= land, fill="lightgray", size=0.2) + #, color="black
  geom_sf(data= lvl1.clusters.pt,  aes(colour = colors.lvl1), size=0.8,  show.legend = TRUE) +
  scale_colour_manual(name= "Bioregions\nlevel 1",
                      values = c("#E69F00", "#66CC99", "#0072B2", "#CC6666","#333333")) + 
  coord_sf(
    xlim = c(limsSI["xmin"], limsSI["xmax"]),
    ylim = c(limsSI["ymin"], 20)
  ) +
  theme_bw() +
  theme(
    panel.background = element_rect(colour = "black")
  ) + labs(x= "Longitude", y="Latitude") +
  guides(colour = guide_legend(override.aes = list(size=2)))



lvl2.clusters$colors.lvl2 <- factor(lvl2.clusters$colors.lvl2, 
                                    levels=c("#8DD3C7","#FB8072", "#80B1D3", "#FFFFB3", "#B3DE69", 
                                             "#D9D9D9", "#BEBADA", "#FCCDE5", "#FDB462",  "black"),
                                    labels=c("1.1",  "1.2", "1.3",  "1.5", "1.12", "1.23",  "2.1", "2.3",  "3.1", "Other" ))
lvl2.clusters.pt <- st_centroid(lvl2.clusters)


lvl2 <- ggplot(lvl2.clusters.pt) + 
  geom_sf(data = siofa, fill = NA, size = 0.3, color="black") +
  geom_sf(data= land, fill="lightgray", size=0.2) + 
  geom_sf(data= lvl2.clusters.pt, aes(colour = colors.lvl2), size=0.8,  show.legend = TRUE) +
  scale_colour_manual(name= "Bioregions\nlevel 2",
                      values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442",  "#0072B2", 
                                 "#b2df8a", "#CC79A7", "#9999CC", "#E41A1C", "black")
  ) +
  coord_sf(
    xlim = c(limsSI["xmin"], limsSI["xmax"]),
    ylim = c(limsSI["ymin"], 20)
  ) +
  theme_bw() +
  theme(
    panel.background = element_rect(colour = "black")
  ) + labs(x= "Longitude", y="Latitude") +
  guides(colour = guide_legend(override.aes = list(size=2)))


saveRDS(lvl1.clusters.pt, "./data/lvl1_pts.RDS")
saveRDS(lvl2.clusters.pt, "./data/lvl2_pts.RDS")


# STEP 7: ENDEMISM
# IndVal, DiVal and robustness
projection.sio <- "+proj=aea +lat_1=0 +lat_2=-52 +lat_0=-30 +lon_0=80 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
grid.aea <- st_transform(grid_poly, projection.sio)

site.area <- data.frame(name = grid.aea$cellid,
                        area = st_area(grid.aea)) #[m^2]

site.area$area <- site.area$area/(1000^2) # from [m^2] to [km^2]

lvl.metrics <- clusterMetrics(db = vme.df, # Database of step 1
                              network = vmes.tree, # Clustering results
                              site.field = 2, # Name of site column in your database
                              species.field = 1, # Name of species column in your database
                              site.area = site.area,  # A data.frame containing at least two columns with site name ("name") and site surface ("surface")
                              level = "lvl1",  # change here your hierarchical level
                              stats.per.cluster = FALSE) # Clustering level at which we calculate metrics



lvl.metrics$network.stats # check the network
# lvl1:
# nb.nodes nb.species.network nb.sites.network nb.species.db nb.sites.db nb.links nb.endemics max.level
# 1     2610               1997              613          1997         613    12882        1903         4

# lvl2:
# nb.nodes nb.species.network nb.sites.network nb.species.db nb.sites.db nb.links nb.endemics max.level
# 1     2610               1997              613          1997         613    12882        1024         4

lvl.metrics$region.stats

# First, order the species metrics by decreasing IndVal values
lvl.metrics$species.stats <- lvl.metrics$species.stats[
  order(lvl.metrics$species.stats$IndVal, decreasing = TRUE), ]


# Create a summary table of results per region (do it for each hierarchical level: change line 217!)
region.stats <- lvl.metrics$region.stats

region.richness <- aggregate(cbind(region.stats$richness, region.stats$char.richness, region.stats$end.richness), by=list(Category=region.stats$cluster), FUN=sum)
colnames(region.richness) <- c("cluster", "richness", "char.richness", "end.richness")
region.richness$perc.end <- (region.richness$end.richness/region.stats$richness)*100
region.richness$perc.sps.world <- (region.richness$richness/sum(region.stats$char.richness))*100

 

# Store Fidelity values of species to each bioregion 
sps.fidelity <- lvl.metrics$species.stats # remember that you need to redo it by hierarchical level (line 217)!

sps.fidelity <- sps.fidelity[sps.fidelity$cluster %in% lvl1.list,] # use this or line 252 depending on what hierarchical level you are analysing
sps.fidelity <- sps.fidelity[sps.fidelity$cluster %in% lvl2.list,]

write.csv(sps.fidelity, "./data/sps.fidelity.lvl1.csv")
write.csv(sps.fidelity, "./data/sps.fidelity.lvl1.csv")







