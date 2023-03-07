###                                                             ###
###   Script to reproduce the main figures of the manuscript    ###
###                                                             ###


library(raster)
library(sf)
library(dplyr)


# 1a. Data -----

# Read in occurrence data
taylor <- read.csv("./data/Taylor_and_Rogers_2017.csv")
taylor <-  taylor[,c(1,2,5,8,11,12,13:15)]

taxa_list <- readRDS("./data/taxa_list.RDS") # checked list of species names 
taxa_data <- droplevels(readRDS("./data/taxa_data.RDS")) # dataset

taxon.level <- "species" # bioregionalization based on species level
taxa_list <- taxa_list[taxa_list$type == taxon.level, ]
taxa_data <- taxa_data[which((taxa_data %>% pull(taxon.level)) %in% taxa_list$taxon), ] 
taxa_data <- taxa_data[,c(29, 7,8,9, 11, 16, 18, 20, 15)]

taxa_data <- rbind(taxa_data, taylor)


# Read in spatial data for context
siofa <- read_sf("./data/spatial/FAO_RFB_SIOFA/RFB_SIOFA.shp") %>% st_as_sf() 
land <- readRDS("./data/spatial/ne_10m_land/land")

iho <- read_sf("./data/spatial/iho/iho.shp") %>% st_as_sf() # Indian Ocean
limsSI <- st_buffer(iho, dist = 0.7) %>% st_bbox()

# Read in some template environmental data at 1D to use as a grid, although
# there are ways to create a grid from scratch
baseline <- readRDS("./data/spatial/baseline.RDS")


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
cell_obs_richness <- data.frame()

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
  
  
  cell_obs_richness <- rbind.data.frame(cell_obs_richness,
                                        data.frame(
                                          Cell_id = cell,
                                          Observed_richness = cell_richness
                                        )
  )
}


# 1c. Plots ----
library(viridis)
library(RColorBrewer)
library(ggplot2)


## Panel a  ----

### First create an empty polygon to fill with richness values
empty.grid <- raster(baseline[[1]]) # This will be our species stack, it is empty at first
vals <- 1:ncell(baseline[[1]])
SIGrid <- setValues(empty.grid, vals)
SIGrid <- rasterToPolygons(SIGrid)
SIGrid <- st_as_sf(SIGrid)
names(SIGrid)[1] <- "cellid"

### Filter the polygon grid with richness

SIGrid$Observed_richness <- NA

SIGrid$Observed_richness <- cell_obs_richness$Observed_richness[match(SIGrid$cellid,
                                                                      cell_obs_richness$Cell_id)]

### Plot

SIGrid_index <- SIGrid[!is.na(SIGrid$Observed_richness),]

S_thres <- ggplot() +
  geom_sf(data = siofa, fill = NA, aes(color=as.factor(RFB)), size = 0.5, show.legend = "line") + 
  scale_color_manual(values = 'grey0', name="", guide = guide_legend(order=2)) +
  geom_sf(data= land, color=NA, fill="lightgray") +
  geom_sf(data = SIGrid_index, aes(fill = Observed_richness), size = 0.1) +
  scale_fill_gradientn(colours = c('#5749a0', '#0f7ab0', '#00bbb1',
                                   '#bef0b0', '#fdf4af', '#f9b64b',
                                   '#ec840e', '#ca443d', '#a51a49'), 
                       name= "Observed\nrichness", guide = guide_colorbar(order = 1)) +
  coord_sf(
    xlim = c(limsSI["xmin"], limsSI["xmax"]),
    ylim = c(limsSI["ymin"], limsSI["ymax"])
  ) +
  theme_bw() +
  labs(x= "Longitude", y="Latitude") 


## Panel b  ----

# With results at 1degree from completeness profiles
results <- readRDS("./data/results_completeness_cellwise_res1.RDS")

results$richness <- "all"

richness_plot <-  ggplot() +
  geom_boxplot(data=results, aes(y = D0, x=richness),  colour= "grey0", fill=  "#263656") +
  annotate("text",
           x = 1:length(table(results$richness)),
           y = aggregate(D0 ~ richness, results, max)[ , 2],
           label = table(results$richness),
           col = "black", vjust = - 0.5, size=3) +
  geom_boxplot(data= results[results$undersampled_cells==TRUE,],
               aes(y = D0, x = undersampled_cells), colour= "grey0", fill= "#95B1C9") + 
  geom_boxplot(data= results[results$undersampled_cells==FALSE,],
               aes(y = D0, x = undersampled_cells), colour= "grey0", fill= "#566C7D") + 
  annotate("text",
           x = 2:length(table(results$undersampled_cells)),
           y = aggregate(D0 ~ undersampled_cells, results, max)[ 1, 2],
           label = table(results$undersampled_cells)[2],
           col = "black", vjust = - 0.5, size=3) +
  annotate("text", 
           x=3, 
           y=aggregate(D0 ~ undersampled_cells, results, max)[2 , 2],
           label = table(results$undersampled_cells)[1],
           col = "black", vjust = - 0.5, size=3) +
  scale_x_discrete(labels = c("All \ncells", "Undersampled \ncells", "Adequately \nsampled cells")) + #c("All \ncells", "Cells with \n> 5 species")
  theme_classic()+
  theme(
    axis.text = element_text(size = 10),
    axis.line = element_line(color = "black", #"#6699FF",
                             arrow = arrow(type='closed', length = unit(5,'pt'))),
    plot.title = element_text(hjust = 0)) +
  labs(y=NULL, x=NULL) +
  ggtitle("Richness")

richness_plot



## Panel c  ----
results <- readRDS("./data/results_completeness_cellwise_res1.RDS")

ggcomp <- pivot_longer(results,
                       cols = c("completeness_q0", "completeness_q1", "completeness_q2"))
ggcomp$q <- as.numeric(as.factor(ggcomp$name)) - 1


C_all <- ggplot(ggcomp, aes(x = q, y = value, group = q)) +
  geom_boxplot(colour= "grey0", fill= "#263656") +
  scale_x_continuous(breaks = c(0, 1, 2)) +
  ylab("Sample completeness") +
  theme_classic()+
  theme(
    axis.text = element_text(size = 10),
    axis.line = element_line(color = "black",
                             arrow = arrow(type='closed', length = unit(5,'pt'))),
    plot.title = element_text(hjust = 0)) +
  ggtitle("Completeness profiles (all cells)")


C_filtered <- ggplot(ggcomp[!ggcomp$undersampled_cells, ], aes(x = q, y = value, group = q)) +
  geom_boxplot(colour= "grey0", fill= "#95B1C9") +
  scale_x_continuous(breaks = c(0, 1, 2)) +
  ylab("Sample completeness") +
  theme_classic()+
  theme(
    axis.text = element_text(size = 10),
    axis.line = element_line(color = "black",
                             arrow = arrow(type='closed', length = unit(5,'pt'))),
    plot.title = element_text(hjust = 0)) +
  ggtitle("Completeness profiles (adequately sampled cells)")


## Arrange figure ----
library(cowplot)

fig1 <- ggdraw() +
  draw_plot(S_thres, x = 0, y = .5, width = .65, height = .5, scale=0.98, hjust=.015) +
  draw_plot(richness_plot, x = .55, y = .5, width = .41, height = 0.5, scale=0.89, hjust=-.15) +
  draw_plot(C_all, x = 0, y = 0, width = .5, height = 0.5, scale=0.95, hjust=.015) +
  draw_plot(C_filtered, x = .5, y = 0, width = .5, height = 0.5, scale=0.95, hjust=.035) +
  draw_plot_label(label = c("A", "B", "C", "D"), size = 15,
                  x = c(0, 0.62, 0,  0.48), y = c(1, 1, 0.5, 0.5))


png("./figures/fig1.png", width = 550 * 5.04, 
    height = 550 * 4.2, res = 300) #* 4.2
print(fig1)
dev.off()



# Figure 2: Figure Bioregions --------------
library(dplyr)
library(sf)
library(ggplot2)


# Map of occurrences in SIOFA
siofa <- read_sf("./data/spatial/FAO_RFB_SIOFA/RFB_SIOFA.shp") %>% st_as_sf() 
land <- readRDS("./data/spatial/ne_10m_land/land")

iho <- read_sf("./data/spatial/iho/iho.shp")
limsSI <- st_buffer(iho, dist = 0.7) %>% st_bbox()


siofa_line = st_cast(siofa,"LINESTRING") # we cast polygon to line for plotting


# Download online bathymetric data for plotting contour lines
library(marmap)
library(ggplot2)
library(tidyr)
library(tibble)

Bathy <- getNOAA.bathy(lon1 = 10, lon2 = 160,
                       lat1 = -65, lat2 = 35, resolution = 30)
Bathy <- as.matrix(Bathy)
class(Bathy) <- "matrix"


Bathy_df <- Bathy %>%
  as.data.frame() %>%
  rownames_to_column(var = "lon") %>%
  gather(lat, value, -1) %>%
  mutate_all(funs(as.numeric)) 


Bathy_df2 <- Bathy_df[!(Bathy_df$value > 0),] # ensure there are only depth values


# Read in occurrence data
lvl1.clusters.pt <- readRDS("./data/lvl1_pts.RDS") # bioregions centroids
lvl2.clusters.pt <- readRDS("./data/lvl2_pts.RDS") # subregions centroids


lvl1 <- ggplot() +
  stat_contour(data= Bathy_df2, aes(x = lon, y = lat, z = value, colour = after_stat(level)),
               bins = 10, #colour = "black", 
               alpha = 0.3, show.legend = FALSE) +
  ggnewscale::new_scale_color() +
  geom_sf(data= siofa_line, aes(colour=as.factor(RFB)), size = 0.5) +
  scale_colour_manual(values = "grey0", name="", guide = guide_legend(order=3)) +
  geom_sf(data= land, color=NA, fill="lightgray") + 
  ggnewscale::new_scale_color() +
  geom_sf(data= lvl1.clusters.pt[!(lvl1.clusters.pt$colors.lvl1) == "Other",], aes(colour= colors.lvl1, shape=colors.lvl1), size=0.9,  show.legend = TRUE) +
  scale_shape_manual(name= "Bioregions\nlevel 1", 
                     labels = c("1", "2", "3", "4"),
                     values=c(16, 17, 18, 15)) + #
  scale_colour_manual(name= "Bioregions\nlevel 1",
                      labels = c("1", "2", "3", "4"),
                      values = c("royalblue4", "#CC6666", "gold2", "#009E73"), 
                      guide=guide_legend(override.aes = list(size = 4))) +
  coord_sf(
    xlim = c(limsSI["xmin"], limsSI["xmax"]),
    ylim = c(limsSI["ymin"], limsSI["ymax"])
  ) +
  theme_bw() + xlab(NULL) + ylab(NULL)


lvl2 <- ggplot() +
  stat_contour(data= Bathy_df2, aes(x = lon, y = lat, z = value, colour = after_stat(level)),
               bins = 10, inherit.aes=FALSE,
               alpha = 0.3) +
  scale_colour_gradient(name= "Depth (m)") +
  geom_sf(data= siofa_line, colour="black", size = 0.5, show.legend = FALSE) +
  geom_sf(data= land, color=NA, fill="lightgray") + 
  ggnewscale::new_scale_color() +
  geom_sf(data= lvl2.clusters.pt[!is.na(lvl2.clusters.pt$colors.lvl2),], aes(colour = colors.lvl2, shape=colors.lvl2), size=1.2,  show.legend = TRUE) +
  scale_shape_manual(name= "Bioregions\nlevel 2", 
                     labels = c("1.1", "1.2", "1.3", "1.5", "1.12", "1.23", "2.1", "2.3", "3.1"),
                     values=c(rep(16,6), 17, 17, 18)) + #
  scale_colour_manual(name= "Bioregions\nlevel 2",
                     labels = c("1.1", "1.2", "1.3", "1.5", "1.12", "1.23", "2.1", "2.3", "3.1"),
                      values = c("dodgerblue4", "royalblue2", "skyblue2", "deepskyblue1" ,"mediumpurple1", "purple3",
                                 "mediumvioletred", "orchid2",  "gold") 
  ) +
  coord_sf(
    xlim = c(limsSI["xmin"], limsSI["xmax"]),
    ylim = c(limsSI["ymin"], limsSI["ymax"])
  ) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_bw() + xlab(NULL) + ylab(NULL) + labs(colour = "Depth")


# arrange figure
library(cowplot)

fig2 <- ggdraw() +
  draw_plot(lvl1, x = 0.1, y = .5, width = 1, height = .5, scale=0.95, hjust= 0.05) +
  draw_plot(lvl2, x = 0.1, y = 0, width = 1, height = .5, scale=0.95, hjust= 0.05) +
  draw_plot_label(label = c("A", "B"), size = 15,
                  x = c(0, 0), y = c(1, 0.5)) +
  draw_label("Longitude", x=0.50, y=  0, vjust=-0.4, angle= 0, size=12) +
  draw_label("Latitude", x= -0.08, y= 0.5, vjust= 13, angle=90, size=12)


png(paste0("./figures/fig_bio.png"), 
    width = 550 * 4.2, height = 550 * 4.5, res = 300) 
print(fig2)
dev.off()



# Figure 3: Figure Fidelity --------------
library(ggridges)
library(forcats)

# Read in fidelity data
sps.fidelity1 <- read.csv("./data/sps.fidelity.lvl1.csv") # fidelity values for bioregions
sps.fidelity2 <- read.csv("./data/sps.fidelity.lvl2.csv") # fidelity values for subregions


# Plot
lvl1.Fi  <- ggplot(sps.fidelity1, aes(x=Occ.Fi, y= as.factor(cluster), fill=as.factor(cluster))) +
  geom_density_ridges(
    jittered_points = TRUE,
    position = position_points_jitter(width = 0.05, height = 0),
    point_shape = '|', point_size = 3, point_alpha = 1, alpha = 0.8,
    show.legend = NA
  ) + 
  scale_fill_manual(values = c("royalblue4", "#CC6666", "gold2", "#009E73")) +
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
  scale_fill_manual(values = c("dodgerblue4", "royalblue2", "skyblue2", "deepskyblue1" ,"mediumpurple1", "purple3",
                               "mediumvioletred", "orchid2",  "gold")) +
  scale_y_discrete(expand = expand_scale(mult = c(.1, .2)), 
                   limits=rev)  + 
  coord_cartesian(clip="on", xlim = c(min(sps.fidelity22$Occ.Fi), 1)) +
  scale_x_continuous(expand = c(0.01, 0.3)) +
  theme_ridges() + + 
  theme(legend.position = "none") + labs(x=NULL, y=NULL) 


# Arrange figure
fig.Fi <- ggdraw() +
  draw_plot(lvl1.Fi, x = 0, y = 0, width = .5, height = 1, scale=0.95) +
  draw_plot(lvl2.Fi, x = .5, y = 0, width = .5, height = 1, scale=0.95) +
  draw_plot_label(label = c("A", "B"), size = 15,
                  x = c(0, 0.5), y = c(1, 1)) +
  draw_label("Fidelity of species", x=0.5, y=  0, vjust=-0.4, angle= 0, size=12) +
  draw_label("Bioregions", x= 0, y= 0.5, vjust=1.2, angle=90, size=12)


png(paste0("./figures/new_fig.Fi.png"), width = 550 * 4.2, height = 550 * 4.2, res = 300)
print(fig.Fi)
dev.off()

