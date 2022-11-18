# This script plots Figure 1, Figure 2 and Figure A4


#   1. Figure 1 : Richness map + Richness & Completeness graphs ------------------

library(raster)
library(sf)
library(dplyr)


# 1a. Data -----

# Read in occurrence data
taylor <- read.csv("./data/Taylor_and_Rogers_2017_network_analysis.csv") # records from the literature
taylor <-  taylor[,c(1,2,5,8,11,12,13:15)]  # select some columns to merge this dataset with the main one (OBIS + GBIF + MNHN + siofa + smithsonian)

taxa_list <- readRDS("./data/taxa_list.RDS") # curated taxonomic taxa list
taxa_data <- droplevels(readRDS("./data/taxa_data_proj.RDS"))  #  (OBIS + GBIF + MNHN + siofa + smithsonian)

taxon.level <- "species" # we'll be working at species level
taxa_list <- taxa_list[taxa_list$type == taxon.level, ] # select only species-level data
taxa_data <- taxa_data[which((taxa_data %>% pull(taxon.level)) %in% taxa_list$taxon), ] 
taxa_data <- taxa_data[,c(29, 7,8,9, 11, 16, 18, 20, 15)] # select same columns to merge both datasets

taxa_data <- rbind(taxa_data, taylor)


# Read in spatial data for context
siofa <- read_sf("./data/spatial/FAO_RFB_SIOFA/RFB_SIOFA.shp") %>% st_as_sf() # siofa polygon boundary
land <- readRDS("./data/spatial/ne_10m_land/land")

iho <- read_sf("./data/spatial/iho/iho.shp") %>% st_as_sf() # Indian Ocean
limsSI <- st_buffer(iho, dist = 0.7) %>% st_bbox() # creates a buffer to use as limits when plotting the maps

# Read in environmental data at 1D
baseline <- readRDS("./data/env/envData_all.RDS") # I use this to obtain cellid numbers to create the 'species-by-site' matrix, 
                                                #but this can be achieved by creating a grid and assigning numbers to each cell.



# 1b. Calculate observed richness and completeness per grid cell -----

library(fossil)
library(ggplot2)

# Aggregate records to 1D cells
taxa_data <- SpatialPointsDataFrame(taxa_data[,c("decimalLongitude", "decimalLatitude")],
                                    data= taxa_data,
                                    proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))

taxa_data$cellid <- cellFromXY(baseline, taxa_data) # extract cell number
taxa_data <- as.data.frame(taxa_data) 

# Calculate richness & completeness
taxa_data$date_richness <- as.Date(taxa_data$ObservationDate)

# keep only those records that have an Observation Date
sample_events <- taxa_data[which(!is.na(taxa_data$ObservationDate)),]
plyr::count(nchar(sample_events$ObservationDate)) # check again whether there are records with no date

index.obs.date <- which(nchar(sample_events$ObservationDate)==0)
sample_events <- sample_events[-index.obs.date,]
plyr::count(nchar(sample_events$ObservationDate)) # check again it is okay

# loop through all the cell samples
cell_completeness <- data.frame()

for(cell in unique(sample_events$cellid))
{
  cell_events <- sample_events[which(sample_events$cellid %in% cell), ]
  rich_cell <- table(cell_events$species,
                     cell_events$ObservationDate)
  rich_cell[rich_cell > 0] <- 1
  
  #- Richness per cell
  richness_per_cell <- colSums(rich_cell)
  
  #- Occurrence for the current phylum
  occur_cell <- rowSums(rich_cell)
  
  #- Cell richness
  cell_richness <- length(occur_cell[which(occur_cell > 0)])
  
  #- Cell richness estimators
  chao_ind <- chao2(rich_cell, taxa.row = TRUE)
  ice_ind <- ICE(rich_cell, taxa.row = TRUE)
  jack1_ind <- jack1(rich_cell, taxa.row = TRUE, abund = FALSE)
  
  #- Average estimator
  avg_ind <- round(mean(c(chao_ind, ice_ind, jack1_ind)))
  
  # completeness index
  completeness <- round(cell_richness/avg_ind, digits = 2)
  
  cell_completeness <- rbind.data.frame(cell_completeness,
                                        data.frame(
                                          Cell_id = cell,
                                          Observed_richness = cell_richness,
                                          chao2 = chao_ind,
                                          ice = ice_ind,
                                          jack1 = jack1_ind,
                                          Average_estimated_richness = avg_ind,
                                          Completeness_index = completeness))
}



# 1c. Plots ----
library(viridis)
library(RColorBrewer)
library(ggplot2)


# Plot richness

## first create an empty polygon to fill with richness values. We'll use it to populate it with calculated completeness/richness values
empty.grid <- raster(baseline[[1]]) # This will be our species stack, it is empty at first
vals <- 1:ncell(baseline[[1]]) # get cell id values
SIGrid <- setValues(empty.grid, vals) # assign values to your raster
SIGrid <- rasterToPolygons(SIGrid) # convert to a polygon
SIGrid <- st_as_sf(SIGrid) # sf object
names(SIGrid)[1] <- "cellid"

## filter the polygon grid with richness and completeness data
SIGrid$Observed_richness <- NA
SIGrid$Completeness_index <- NA

SIGrid$Observed_richness <- cell_completeness$Observed_richness[match(SIGrid$cellid,
                                                                      cell_completeness$Cell_id)]
SIGrid$Completeness_index <- cell_completeness$Completeness_index[match(SIGrid$cellid,
                                                                        cell_completeness$Cell_id)]


### plot with no richness threshold

SIGrid_index <- SIGrid[!is.na(SIGrid$Observed_richness),]

S_nothres <- ggplot() +
  geom_sf(data = siofa, fill = NA, aes(color=as.factor(RFB)), size = 0.5, show.legend = "line") + 
  geom_sf(data= land, color=NA, fill="lightgray") +
  geom_sf(data = SIGrid_index, aes(fill = Observed_richness), size = 0.1) +
  scale_fill_viridis( name= "Observed\nrichness",  guide = guide_colorbar(order = 1)) +
  scale_color_manual(values = 'grey0', name="", guide = guide_legend(order=2)) +
  coord_sf(
    xlim = c(limsSI["xmin"], limsSI["xmax"]),
    ylim = c(limsSI["ymin"], limsSI["ymax"])
  ) +
  theme_bw() +
  labs(x= "Longitude", y="Latitude") 




# Plot frequency of richness per cell using two cut-offs 
S <- 30 

threshold <- ifelse(cell_completeness$Observed_richness > S, paste("S >", S), paste("S <=", S))
cell_completeness$S_threshold <- threshold
cell_completeness$richness <- "all"

richness_plot <-  ggplot() +
  geom_boxplot(data=cell_completeness, aes(y = Observed_richness, x=richness),  colour= "#263656", fill=  "#263656") +
  annotate("text",
           x = 1:length(table(cell_completeness$richness)),
           y = aggregate(Observed_richness ~ richness, cell_completeness, max)[ , 2],
           label = table(cell_completeness$richness),
           col = "black", vjust = - 0.5, size=3) +
  geom_boxplot(data= cell_completeness,
               aes(y = Observed_richness, x = S_threshold), colour= "#566C7D", fill= "#566C7D") + 
  geom_boxplot(data= cell_completeness[cell_completeness$S_threshold == paste("S >", S),],
               aes(y = Observed_richness, x = S_threshold), colour= "#95B1C9", fill= "#95B1C9") + 
  annotate("text",
           x = 2:length(table(cell_completeness$S_threshold)),
           y = aggregate(Observed_richness ~ S_threshold, cell_completeness, max)[ 1, 2],
           label = table(cell_completeness$S_threshold)[1],
           col = "black", vjust = - 0.5, size=3) +
  annotate("text", 
           x=3, y=aggregate(Observed_richness ~ S_threshold, cell_completeness, max)[2 , 2],
           label = table(cell_completeness$S_threshold)[2],
           col = "black", vjust = - 0.5, size=3) +
  scale_x_discrete(labels = c("All \ncells", paste("Cells with \n<=", S, "species"), paste("Cells with \n>", S, "species"))) + 
  theme_classic()+
  theme(
    axis.text = element_text(size = 10),
    axis.line = element_line(color = "#6699FF",
                             arrow = arrow(type='closed', length = unit(5,'pt'))),
    plot.title = element_text(hjust = 0)) +
  labs(y=NULL, x=NULL) +
  ggtitle("Richness")



# Plot frequency of completeness per cell using two cut-offs 

# make filters here

cell_completeness$Completeness_index_thres <- cell_completeness$Completeness_index
cell_completeness$Completeness_index_thres[
  which(cell_completeness$Observed_richness <= S)] <- 0 


completeness_plot <- 
  ggplot() +
  geom_boxplot(data=cell_completeness, aes(y = Completeness_index , x=richness),  colour= "#263656", fill=  "#263656") + 
  annotate("text",
           x = 1:length(table(cell_completeness$richness)),
           y = aggregate(Completeness_index ~ richness, cell_completeness, max)[ , 2],
           label = table(cell_completeness$richness),
           col = "black", vjust = - 0.5, size=3) +
  geom_boxplot(data= cell_completeness,
               aes(y = Completeness_index, x = S_threshold), colour= "#566C7D", fill= "#566C7D") + 
  geom_boxplot(data= cell_completeness[cell_completeness$S_threshold == paste("S >", S),],
               aes(y = Completeness_index_thres, x = S_threshold),  colour= "#95B1C9", fill= "#95B1C9") + 
  annotate("text",
           x = 2:length(table(cell_completeness$S_threshold)),
           y = aggregate(Completeness_index ~ S_threshold, cell_completeness, max)[ 1, 2],
           label = table(cell_completeness$S_threshold)[1],
           col = "black", vjust = - 0.5, size=3) +
  annotate("text", 
           x=3, y=aggregate(Completeness_index ~ S_threshold, cell_completeness, max)[2 , 2],
           label =table(cell_completeness$S_threshold)[2],
           col = "black", vjust = - 0.5, size=3) +
  scale_x_discrete(labels = c("All \ncells", paste("Cells with \n<=", S, "species"), paste("Cells with \n>", S, "species"))) + 
  theme_classic()+
  theme(
    axis.text = element_text(size = 10),
    axis.line = element_line(color = "#6699FF",
                             arrow = arrow(type='closed', length = unit(5,'pt'))),
    plot.title = element_text(hjust = 0)) +
  labs(y=NULL, x=NULL) +
  ggtitle("Completeness")



# Arrange panels
library(gridExtra)
library(cowplot)
library(ggpubr)


# to draw at particular locations
fig1 <- ggdraw() +
  draw_plot(S_nothres, x = 0, y = .5, width = 1, height = .5) +
  draw_plot(completeness_plot, x = .5, y = 0, width = .5, height = .5) +
  draw_plot(richness_plot, x = 0, y = 0, width = .5, height = 0.5) +
  draw_plot_label(label = c("A", "B", "C"), size = 15,
                  x = c(0, 0, 0.5), y = c(1, 0.5, 0.5))





#   2. Figure 2 : Bioregions ------------------

library(dplyr)
library(sf)
library(ggplot2)


# Read in spatial data for context 
siofa <- read_sf("./data/spatial/FAO_RFB_SIOFA/RFB_SIOFA.shp") %>% st_as_sf() 
land <- readRDS("./data/spatial/ne_10m_land/land")

iho <- read_sf("./data/spatial/iho/iho.shp") %>% st_as_sf() # Indian Ocean
limsSI <- st_buffer(iho, dist = 0.7) %>% st_bbox()


# Read in occurrence data
lvl1.clusters.pt <- readRDS("./data/lvl1_pts.RDS") # obtained from 1.network_analysis.R
lvl2.clusters.pt <- readRDS("./data/lvl2_pts.RDS") # obtained from 1.network_analysis.R


# Plot

lvl1 <- ggplot() +
  geom_sf(data = siofa, fill = NA, size = 0.5,  show.legend = F) +
  geom_sf(data= land, color=NA, fill="lightgray") + 
  geom_sf(data= lvl1.clusters.pt[!(lvl1.clusters.pt$colors.lvl1) == "Other",], aes(colour= colors.lvl1), size=0.9,  show.legend = TRUE) +
  scale_colour_manual(name= "Bioregions\nlevel 1",
                      values = c("#3A4960", "#91323A", "#D7C969", "#669900"), 
                      guide=guide_legend(override.aes = list(size = 2))) +
  coord_sf(
    xlim = c(limsSI["xmin"], limsSI["xmax"]),
    ylim = c(limsSI["ymin"], limsSI["ymax"])
  ) +
  theme_bw() 


lvl2 <- ggplot() +
  geom_sf(data = siofa, fill = NA, size = 0.5) +
  geom_sf(data= land, color=NA, fill="lightgray") + 
  geom_sf(data= lvl2.clusters.pt[!is.na(lvl2.clusters.pt$colors.lvl2),], aes(colour = colors.lvl2), size=0.9,  show.legend = TRUE) +
  scale_colour_manual(name= "Bioregions\nlevel 2",
                      values = c("#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733","#999933" ,
                                 "#882255", "#AA4499",  "#DDCC77")
  ) +
  coord_sf(
    xlim = c(limsSI["xmin"], limsSI["xmax"]),
    ylim = c(limsSI["ymin"], limsSI["ymax"])
  ) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  theme_bw() 


# arrange figure
library(cowplot)

fig2 <- ggdraw() +
  draw_plot(lvl1, x = 0.1, y = .5, width = 1, height = .5, scale=0.95) +
  draw_plot(lvl2, x = 0.1, y = 0, width = 1, height = .5, scale=0.95) +
  draw_plot_label(label = c("A", "B"), size = 15,
                  x = c(0, 0), y = c(1, 0.5)) +
  draw_label("Longitude", x=0.55, y=  0, vjust=-0.4, angle= 0, size=12) +
  draw_label("Latitude", x= 0, y= 0.5, vjust= 13, angle=90, size=12)




#   3. Figure 3 : Fidelity ------------------
library(ggridges)
library(forcats)

# Read fidelity data
sps.fidelity1 <- read.csv("./data/sps.fidelity.lvl1.csv") # obtained from 1.network_analysis.R
sps.fidelity2 <- read.csv("./data/sps.fidelity.lvl2.csv") # obtained from 1.network_analysis.R


# Plot
lvl1.Fi  <- ggplot(sps.fidelity1, aes(x=Occ.Fi, y= as.factor(cluster), fill=as.factor(cluster))) +
  geom_density_ridges(
    jittered_points = TRUE,
    position = position_points_jitter(width = 0.05, height = 0),
    point_shape = '|', point_size = 3, point_alpha = 1, alpha = 0.8,
    show.legend = NA
  ) + 
  scale_fill_manual(values = c("#3A4960", "#91323A", "#D7C969", "#669900")) +
  scale_y_discrete(expand = expand_scale(mult = c(.1, .5)),
                   limits=rev)  +
  coord_cartesian(clip = "on", xlim = c(min(sps.fidelity1$Occ.Fi), 1)) +
  scale_x_continuous(expand = c(0.01, 0.2)) +
  theme_ridges() + 
  theme(legend.position = "none") + labs(x=NULL, y=NULL) 


lvl2.list <- c("1.1",  "1.2", "1.3",  "1.5", "1.12", "1.23",  "2.1", "2.3",  "3.1")
sps.fidelity22 <- sps.fidelity2[sps.fidelity2$cluster %in% lvl2.list,]
sps.fidelity22$cluster <- factor(sps.fidelity2$cluster, 
                                 levels=c("1.1", "1.12",  "1.2", "1.23", "1.3", "1.5", "2.1", "2.3", "3.1"),
                                 labels=c("1.1",  "1.2", "1.3",  "1.5", "1.12", "1.23",  "2.1", "2.3",  "3.1"))


lvl2.Fi <- ggplot(data=sps.fidelity22, aes(x=Occ.Fi, y= cluster, fill=cluster)) +
  geom_density_ridges(
    jittered_points = TRUE,
    position = position_points_jitter(width = 0.05, height = 0),
    point_shape = '|', point_size = 3, point_alpha = 1, alpha = 0.8,
    show.legend = NA
  ) + 
  scale_fill_manual(values =c("#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733","#999933" ,
                              "#882255", "#AA4499",  "#DDCC77")) +
  scale_y_discrete(expand = expand_scale(mult = c(.1, .2)), 
                   limits=rev)  + 
  coord_cartesian(clip="on", xlim = c(min(sps.fidelity22$Occ.Fi), 1)) +
  scale_x_continuous(expand = c(0.01, 0.3)) +
  theme_ridges() + 
  theme(legend.position = "none") + labs(x=NULL, y=NULL) 



# arrange figure

fig.Fi <- ggdraw() +
  draw_plot(lvl1.Fi, x = 0, y = 0, width = .5, height = 1, scale=0.95) +
  draw_plot(lvl2.Fi, x = .5, y = 0, width = .5, height = 1, scale=0.95) +
  draw_plot_label(label = c("A", "B"), size = 15,
                  x = c(0, 0.5), y = c(1, 1)) +
  draw_label("Fidelity of species", x=0.5, y=  0, vjust=-0.4, angle= 0, size=12) +
  draw_label("Bioregions", x= 0, y= 0.5, vjust=1.2, angle=90, size=12)




