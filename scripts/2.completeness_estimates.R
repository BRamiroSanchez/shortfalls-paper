
# This script contains the calculation of the completeness of our inventory (Table 1),
# and the sensitivity analyses described in the Online Appendix (Figure A1, Figure A2, Figure A3).
  
# It consists of the same calculation of the completeness index at different resolutions
# and scales.

# The minimum information needed to estimate the completeness of the inventory is
# species, coordinates, and observation date (Soberón et al., 2007)
# Soberón, J., Jiménez, R., Golubov, J., & Koleff, P. (2007). Assessing completeness of 
# biodiversity databases at different spatial scales. Ecography, 30(1), 152-160.



  # Load libraries
  list.of.packages <- c("ggplot2", "tidyverse", "sp", "rgdal", "sf", "rnaturalearth", 
                        "fossil",  "RColorBrewer", "viridis", "reshape",
                        "dplyr")
  lapply(list.of.packages, library, character.only = TRUE) #load all the packages
  
  
  # Background spatial data polygons (for plots later on)
  
  SIpoly <- list(rbind(c(13, -65), c(13, 20), c(147, 20), c(147, -65), c(13, -65))) %>% 
    st_polygon() %>% st_sfc(crs = 4326) %>% st_sf()
  
  worldMap <- ne_countries(scale = "medium", type = "countries", returnclass = "sf")
  iho <- read_sf("./data/spatial/iho/iho.shp") %>% st_as_sf() # Indian Ocean
  limsSI <- st_buffer(iho, dist = 0.7) %>% st_bbox()
  
  siofa <- read_sf("./data/spatial/FAO_RFB_SIOFA/RFB_SIOFA.shp") %>% st_as_sf() 
  land <- readRDS("./data/spatial/ne_10m_land/land")
  
  
  
  
  #   1. Completeness for the whole area based on inventories   (Table 1) -----------
  
  
  # Read in your occurrence data
  taylor <- read.csv("./data/Taylor_and_Rogers_2017_network_analysis.csv") # literature data
  taylor <-  taylor[,c(1,5,8,11,12,13:15)]  # I'm selecting certain columns to merge with "taxa_data" dataset
  
  taxa_list <- readRDS("./data/taxa_list.RDS") # curated taxa list of taxonomic names only
  taxa_data <- droplevels(readRDS("./data/taxa_data_proj.RDS")) # obis + gbif + smithsonian + siofa + mnhn
  taxa_data <- taxa_data[which(taxa_data$species %in% taxa_list$taxon), ] # data only at species level (match with curated list)
  taxa_data <- taxa_data[,c(29, 7,8,9, 16, 18, 20, 15)]  # select columns to merge with the "taylor" dataset
  
  taxa_data <- rbind(taxa_data, taylor)
  
  # Create a grid (this will be the base for a species-by-site matrix by using 'cellid')
  SIGrid <- SIpoly %>%
    st_make_grid(cellsize = c(1,1)) %>%
    st_sf() %>%
    mutate(cellid = row_number())
  
  taxa_data <- st_as_sf(taxa_data, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
  taxa_data <- taxa_data %>% st_join(SIGrid) 
  
  # Look for ObservationDate (i.e. inventory)
  taxa_data$date_richness <- as.Date(taxa_data$ObservationDate)
  
  # Calculate completeness
  ocurr_sample_events <- taxa_data[which(!is.na(taxa_data$ObservationDate)),] # make sure we retain data with an observation date
  plyr::count(nchar(ocurr_sample_events$ObservationDate)) # there are some records with no date
  index.obs.date <- which(nchar(ocurr_sample_events$ObservationDate)==0) # identify them
  ocurr_sample_events <- ocurr_sample_events[-index.obs.date,] # and remove them
  plyr::count(nchar(ocurr_sample_events$ObservationDate)) # double check now it is all good to go
  
  
  #- Create a contingency table for the entire study area
  contin_whole_area <- table(ocurr_sample_events$species, 
                             ocurr_sample_events$cellid)
  
  contin_whole_area[1:10, 1:10]
  
  #- Currently it counts the number of records per cell - we want to switch to presence/absence per cell:
  contin_whole_area[contin_whole_area > 0] <- 1
  contin_whole_area[1:10, 1:10]
  
  #- Get species richness per cell
  cell_richness <- colSums(contin_whole_area)
  
  #- Get species occurrence at the current resolution
  species_occurrence <- rowSums(contin_whole_area)
  
  #- Have a quick look at the occurrence table to check that everything is alright
  hist(species_occurrence, breaks = 100) #- L-shaped, everything is fine / expected
  
  #- Check whether there are species that do not seem to occur in the study area. 
  #- They should be excluded from richness counts
  species_occurrence[which(species_occurrence == 0)]
  
  #- Total area richness is the number of rows of the contingency table
  total_area_richness <- length(species_occurrence[which(species_occurrence != 0)]) 
  
  
  # richness estimators {library:fossil} for the whole area
  chao_ind <- chao2(contin_whole_area, taxa.row = TRUE)
  ice_ind <- ICE(contin_whole_area, taxa.row = TRUE)
  jack1_ind <- jack1(contin_whole_area, taxa.row = TRUE, abund = FALSE)
  
  avg_ind <- round(mean(c(chao_ind, ice_ind, jack1_ind)))
  
  
  # completeness index
  completeness <- round(total_area_richness/avg_ind, digits = 2)
  
  # create a results table
  indices <- matrix(c(total_area_richness, avg_ind, completeness), ncol=3)
  
  colnames(indices) <- c("Species_richness","Average_estimated_richness", "Completeness_index")
  rownames(indices) <- "Database" 
  indices
  
  
  #- This table could be completed by completeness per phylum
  for(phy in unique(ocurr_sample_events$phylum))
  {
    if(!(phy %in% c("Arthropoda", "Hemichordata", "Tunicata"))) # Not enough data for these three phyla!
    {
      #- Restrict contingency table only to the current phylum:
      contin_whole_area_cur_phylum <- contin_whole_area[
        which(rownames(contin_whole_area) %in%
                ocurr_sample_events$species[
                  which(ocurr_sample_events$phylum == phy)
                ]),
      ] 
      #- Richness per cell
      richness_cur_phylum <- colSums(contin_whole_area_cur_phylum)
      #- Exclude empty cells
      if(any(richness_cur_phylum == 0))
      {
        contin_whole_area_cur_phylum <- 
          contin_whole_area_cur_phylum[, -which(richness_cur_phylum == 0), drop = FALSE]
        richness_cur_phylum <- richness_cur_phylum[-which(richness_cur_phylum == 0)]
      }
      
      #- Occurrence for the current phylum
      occur_cur_phylum <- rowSums(contin_whole_area_cur_phylum)
      
      total_area_richness_cur_phylum <- length(occur_cur_phylum[which(occur_cur_phylum > 0)])
      
      chao_ind <- chao2(contin_whole_area_cur_phylum, taxa.row = TRUE)
      ice_ind <- ICE(contin_whole_area_cur_phylum, taxa.row = TRUE)
      jack1_ind <- jack1(contin_whole_area_cur_phylum, taxa.row = TRUE, abund = FALSE)
      
      avg_ind <- round(mean(c(chao_ind, ice_ind, jack1_ind)))
      
      # completeness index
      completeness <- round(total_area_richness_cur_phylum/avg_ind, digits = 2)
      
      indices <- rbind.data.frame(indices,
                                  data.frame(
                                    Species_richness = total_area_richness_cur_phylum,
                                    Average_estimated_richness = avg_ind,
                                    Completeness_index = completeness))
      rownames(indices)[nrow(indices)] <- phy
    }
  }
  
  indices

  
  
  
  #   2. Completeness for each cell at varying spatial resolutions  (Figure A1) ----------
  
  
  cellXsize <- 1:7
  cellYsize <- 1:7
  
  cell_completeness <- data.frame()
  
  
  for(i in 1:length(cellXsize)){
    
    SIGrid <- SIpoly %>%
      st_make_grid(cellsize = c(cellXsize[i],cellXsize[i])) %>%
      st_sf() %>%
      mutate(cellid = row_number())
    
    # Read in occurrence data
    taylor <- read.csv("./data/Taylor_and_Rogers_2017_network_analysis.csv")
    taylor <-  taylor[,c(1,5,8,11,12,13:15)]
    
    taxa_list <- readRDS("./data/taxa_list.RDS") 
    taxa_data <- droplevels(readRDS("./data/taxa_data_proj.RDS")) 
    taxa_data <- taxa_data[which(taxa_data$species %in% taxa_list$taxon), ]
    taxa_data <- taxa_data[,c(29, 7,8,9, 16, 18, 20, 15)]# select same columns as the 'other'taylor' dataset
    
    taxa_data <- rbind(taxa_data, taylor)
    
    taxa_data <- st_as_sf(taxa_data, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
    taxa_data <- taxa_data %>% st_join(SIGrid) 
    
    # Look for ObservationDate (i.e. inventory)
    taxa_data$date_richness <- as.Date(taxa_data$ObservationDate)
    
    # Calculate completeness
    ocurr_sample_events <- taxa_data[which(!is.na(taxa_data$ObservationDate)),]
    plyr::count(nchar(ocurr_sample_events$ObservationDate))
    index.obs.date <- which(nchar(ocurr_sample_events$ObservationDate)==0)
    ocurr_sample_events <- ocurr_sample_events[-index.obs.date,]
    plyr::count(nchar(ocurr_sample_events$ObservationDate))
    
    for(cell in unique(ocurr_sample_events$cellid))
    {
      current_cell_events <- ocurr_sample_events[which(ocurr_sample_events$cellid %in% cell), ]
      contin_cur_cell <- table(current_cell_events$species,
                               current_cell_events$ObservationDate)
      contin_cur_cell[contin_cur_cell > 0] <- 1
      
      #- Richness per cell
      richness_cur_cell <- colSums(contin_cur_cell)
      
      #- Occurrence for the current phylum
      occur_cur_cell <- rowSums(contin_cur_cell)
      
      #- Cell richness
      cell_richness <- length(occur_cur_cell[which(occur_cur_cell > 0)])
      
      #- Cell richness estimators
      chao_ind <- chao2(contin_cur_cell, taxa.row = TRUE)
      ice_ind <- ICE(contin_cur_cell, taxa.row = TRUE)
      jack1_ind <- jack1(contin_cur_cell, taxa.row = TRUE, abund = FALSE)
      
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
                                              Completeness_index = completeness,
                                              cell_resolution = cellXsize[i],
                                              taxon = "species")
      )
    }
    
  }
  
  
  # Check Observed richness distribution
  hist(cell_completeness$Observed_richness, breaks=120)
  
  
  # Plots
  
  C <- ggplot(cell_completeness, aes(x=as.factor(cell_resolution), y= Completeness_index)) + 
    geom_boxplot(fill="slateblue", alpha=0.2) + 
    xlab(expression("Cell resolution ("*~degree*")")) + ylab("Completeness") +
    theme_bw() +
    ggtitle('All cells')
  
  # filter only with S
  cell_completeness1 <- cell_completeness
  
  cell_completeness0 <- cell_completeness1[cell_completeness1$Observed_richness >= 30,]
  
  C_more30 <- ggplot(cell_completeness0, aes(x=as.factor(cell_resolution), y= Completeness_index)) + 
    geom_boxplot(fill="slateblue", alpha=0.2) + 
    xlab(expression("Cell resolution ("*~degree*")")) + ylab("Completeness") +
    theme_bw() +
    ggtitle('Cells with S >= 30')
  
  
  cell_completeness1 <- cell_completeness1[cell_completeness1$Observed_richness < 30,]
  
  C30 <- ggplot(cell_completeness1, aes(x=as.factor(cell_resolution), y= Completeness_index)) + 
    geom_boxplot(fill="slateblue", alpha=0.2) + 
    xlab(expression("Cell resolution ("*~degree*")")) + ylab("Completeness") +
    theme_bw() +
    ggtitle('Cells with S < 30')
  
  
  # save plot
  library(cowplot)
  
  plot_grid(C, C_more30, C30, labels=c("A", "B", "C"), ncol = 3, nrow = 1)

  
  
  
  
  #   3. Completeness using an occurrences threshold (Figure A2)  ----------
  
  
  # Now, the same but with a threshold for occurrences
  freqs_per_resolution <- data.frame()
  
  for(i in 1:length(cellXsize)){
    
    SIGrid <- SIpoly %>%
      st_make_grid(cellsize = c(cellXsize[i],cellXsize[i])) %>%
      st_sf() %>%
      mutate(cellid = row_number())
    
    
    taylor <- read.csv("./Taylor_and_Rogers_2017_network_analysis.csv")
    taylor <-  taylor[,c(1,5,8,11,12,13:15)]
    
    taxa_list <- readRDS("./data/taxa_list.RDS") 
    taxa_data <- droplevels(readRDS("./data/taxa_data_proj.RDS")) 
    taxa_data <- taxa_data[which(taxa_data$species %in% taxa_list$taxon), ]
    taxa_data <- taxa_data[,c(29, 7,8,9, 16, 18, 20, 15)] # select same columns as the 'taylor' dataset
    
    taxa_data <- rbind(taxa_data, taylor)
    
    taxa_data <- st_as_sf(taxa_data, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
    taxa_data <- taxa_data %>% st_join(SIGrid) 
    
    
    # Filter here for number of occurrences
    freqs <- taxa_data %>% group_by(cellid) %>% summarise(n=n())
    hist(freqs$n, breaks=100)
    nb_records <- data.frame(
      occurrence_cutoff = c("All",10, 20, 30, "less30"), # some 'arbitrary' cut-off occurrence thresholds
      records_per_cell = c(
        nrow(freqs[freqs$n >= 1, ]),
        nrow(freqs[freqs$n >= 10, ]),
        nrow(freqs[freqs$n >= 20, ]),
        nrow(freqs[freqs$n >= 30, ]),
        nrow(freqs[freqs$n < 30, ])
      )
    )
    
    freqs$occ_qty <- ifelse(freqs$n >= 30, "Good",   # I categorise the cut-offs to find them more easily
                            ifelse(freqs$n < 30, "Fair",
                                   "Insufficient"))
    
    freqs_per_resolution <- rbind.data.frame(freqs_per_resolution,
                                             data.frame(
                                               Cell_id = freqs$cellid,
                                               occ_qty = freqs$occ_qty,
                                               n = freqs$n,
                                               cell_resolution = cellXsize[i])
    )
  }
  
  
  # Add some of this information to the completeness dataframe calculated in Section 2
  
  cell_completeness$occ_qty <- freqs_per_resolution$occ_qty[match(paste(cell_completeness$Cell_id, cell_completeness$cell_resolution), 
                                                                  paste(freqs_per_resolution$Cell_id, freqs_per_resolution$cell_resolution))]
  
  cell_completeness$n <- freqs_per_resolution$n[match(paste(cell_completeness$Cell_id, cell_completeness$cell_resolution), 
                                                      paste(freqs_per_resolution$Cell_id, freqs_per_resolution$cell_resolution))]
  
  
  # Plots:
  
  occ_completeness_good <- occ_completeness_fair <- occ_completeness_low <- occ_completeness_insuf <- cell_completeness
  
  
  occ_good <- ggplot(occ_completeness_good, aes(x=as.factor(cell_resolution), y= Completeness_index)) +
    geom_boxplot(fill="#69b3a2", alpha=0.5)+  
    xlab(expression("Cell resolution ("*~degree*")")) + ylab("Completeness") +
    theme_bw() +
    ggtitle('All cells')
  
  
  occ_completeness_more30 <- occ_completeness_good[occ_completeness_good$occ_qty == "Good",]
  
  occ_more30 <- ggplot(occ_completeness_more30, aes(x=as.factor(cell_resolution), y= Completeness_index)) +
    geom_boxplot(fill="#69b3a2", alpha=0.5)+ 
    xlab(expression("Cell resolution ("*~degree*")")) + ylab("Completeness") +
    theme_bw() +
    ggtitle('Cells with n >= 30')
  
  
  occ_completeness_less30 <- occ_completeness_good[occ_completeness_good$occ_qty == "Fair",]
  
  occ_less30 <- ggplot(occ_completeness_less30, aes(x=as.factor(cell_resolution), y= Completeness_index)) +
    geom_boxplot(fill="#69b3a2", alpha=0.5)+  
    xlab(expression("Cell resolution ("*~degree*")")) + ylab("Completeness") +
    theme_bw() +
    ggtitle('Cells with n < 30')
  

  # arrange plots
  library(cowplot)
  
  plot_grid(occ_good, occ_more30, occ_less30, labels=c("A", "B", "C"), ncol = 3, nrow = 1)

  
  
  
  #   4. SIOFA Completeness for the entire Convention Area (Table 1) ----------
  
  # Read in occurrence data
  taylor <- read.csv("./data/Taylor_and_Rogers_2017_network_analysis.csv")
  taylor <-  taylor[,c(1,5,8,11,12,13:15)]
  
  taxa_list <- readRDS("./data/taxa_list.RDS") 
  taxa_data <- droplevels(readRDS("./data/taxa_data_proj.RDS")) 
  taxa_data <- taxa_data[which(taxa_data$species %in% taxa_list$taxon), ] 
  taxa_data <- taxa_data[,c(29, 7,8,9, 16, 18, 20, 15)]
  
  taxa_data <- rbind(taxa_data, taylor)
  
  # Create SIOFA grid 
  SIGrid <- SIpoly %>%
    st_make_grid(cellsize = c(1,1)) %>%
    #st_intersection(SIpoly) %>%
    st_cast("MULTIPOLYGON") %>%
    st_sf() %>%
    mutate(cellid = row_number())
  
  sf_use_s2(FALSE) 
  
  siofaGrid <- sf::st_intersection(SIGrid, siofa) 
  ggplot() + 
    geom_sf(data = siofa, fill = "blue") +
    geom_sf(data = siofaGrid, color = "red", fill = NA) 
  
  taxa_data <- st_as_sf(taxa_data, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
  taxa_data <- taxa_data %>% st_join(siofaGrid) 
  
  
  # Look for ObservationDate (i.e. inventory)
  taxa_data$date_richness <- as.Date(taxa_data$ObservationDate)
  
  # Calculate completeness
  ocurr_sample_events <- taxa_data[which(!is.na(taxa_data$ObservationDate)),]
  plyr::count(nchar(ocurr_sample_events$ObservationDate))
  index.obs.date <- which(nchar(ocurr_sample_events$ObservationDate)==0)
  ocurr_sample_events <- ocurr_sample_events[-index.obs.date,]
  plyr::count(nchar(ocurr_sample_events$ObservationDate))
  
  
  #- Create a contingency table for the entire study area
  contin_whole_area <- table(ocurr_sample_events$species, #change as required
                             ocurr_sample_events$cellid)
  
  contin_whole_area[1:10, 1:10]
  dim(contin_whole_area)
  
  #- Currently it counts the number of records per cell - we want to switch to presence/absence per cell:
  contin_whole_area[contin_whole_area > 0] <- 1
  contin_whole_area[1:10, 1:10]
  
  #- Get species richness per cell
  cell_richness <- colSums(contin_whole_area)
  
  #- Get species occurrence at the current resolution
  species_occurrence <- rowSums(contin_whole_area)
  
  #- Have a quick look at the occurrence table to check that everything is alright
  hist(species_occurrence, breaks = 100) #- L-shaped, everything is fine / expected
  
  #- Note that we have species that do not seem to occur in the study area. 
  #- They should be excluded from richness counts
  species_occurrence[which(species_occurrence == 0)]
  
  #- Total area richness is the number of rows of the contingency table
  total_area_richness <- length(species_occurrence[which(species_occurrence != 0)]) # S = 118
  
  
  # richness estimators {library:fossil} for the whole area
  chao_ind <- chao2(contin_whole_area, taxa.row = TRUE)
  ice_ind <- ICE(contin_whole_area, taxa.row = TRUE)
  jack1_ind <- jack1(contin_whole_area, taxa.row = TRUE, abund = FALSE)
  
  avg_ind <- round(mean(c(chao_ind, ice_ind, jack1_ind)))
  
  
  # completeness index
  completeness <- round(total_area_richness/avg_ind, digits = 2)
  
  # create a results table
  indices <- matrix(c(total_area_richness, avg_ind, completeness), ncol=3)
  
  colnames(indices) <- c("Species_richness","Average_estimated_richness", "Completeness_index")
  rownames(indices) <- "Database" 
  indices
  
  
  #- This table could be completed by completeness per phylum
  for(phy in unique(ocurr_sample_events$phylum))
  {
    if(!(phy %in% c("Arthropoda", "Hemichordata", "Tunicata", "Annelida"))) # Not enough data for these phyla
    {
      #- Restrict contingency table only to the current phylum:
      contin_whole_area_cur_phylum <- contin_whole_area[
        which(rownames(contin_whole_area) %in%
                ocurr_sample_events$species[
                  which(ocurr_sample_events$phylum == phy)
                ]),
      ] 
      #- Richness per cell
      richness_cur_phylum <- colSums(contin_whole_area_cur_phylum)
      #- Exclude empty cells
      if(any(richness_cur_phylum == 0))
      {
        contin_whole_area_cur_phylum <- 
          contin_whole_area_cur_phylum[, -which(richness_cur_phylum == 0), drop = FALSE]
        richness_cur_phylum <- richness_cur_phylum[-which(richness_cur_phylum == 0)]
      }
      
      #- Occurrence for the current phylum
      occur_cur_phylum <- rowSums(contin_whole_area_cur_phylum)
      
      total_area_richness_cur_phylum <- length(occur_cur_phylum[which(occur_cur_phylum > 0)])
      
      chao_ind <- chao2(contin_whole_area_cur_phylum, taxa.row = TRUE)
      ice_ind <- ICE(contin_whole_area_cur_phylum, taxa.row = TRUE)
      jack1_ind <- jack1(contin_whole_area_cur_phylum, taxa.row = TRUE, abund = FALSE)
      
      avg_ind <- round(mean(c(chao_ind, ice_ind, jack1_ind)))
      
      # completeness index
      completeness <- round(total_area_richness_cur_phylum/avg_ind, digits = 2)
      
      indices <- rbind.data.frame(indices,
                                  data.frame(
                                    Species_richness = total_area_richness_cur_phylum,
                                    Average_estimated_richness = avg_ind,
                                    Completeness_index = completeness))
      rownames(indices)[nrow(indices)] <- phy
    }
  }
  
  indices
  
  
  
  #   5. SIOFA Completeness for each cell at varying spatial resolutions  (Figure A3) ----------------
  
  cellXsize <- 1:7
  cellYsize <- 1:7
  
  siofa_cell_completeness <- data.frame()
  
  
  for(i in 1:length(cellXsize)){
    
    # Read in occurrence data
    taylor <- read.csv("./data/Taylor_and_Rogers_2017_network_analysis.csv")
    taylor <-  taylor[,c(1,5,8,11,12,13:15)]
    
    taxa_list <- readRDS("./data/taxa_list.RDS")
    taxa_data <- droplevels(readRDS("./data/taxa_data_proj.RDS"))
    taxa_data <- taxa_data[which(taxa_data$species %in% taxa_list$taxon), ] 
    taxa_data <- taxa_data[,c(29, 7,8,9, 16, 18, 20, 15)]
    
    taxa_data <- rbind(taxa_data, taylor)
    
    # Create SIOFA grid 
    SIGrid <- SIpoly %>%
      st_make_grid(cellsize = c(cellXsize[i],cellXsize[i])) %>%
      #st_intersection(SIpoly) %>%
      st_cast("MULTIPOLYGON") %>%
      st_sf() %>%
      mutate(cellid = row_number())
    
    sf_use_s2(FALSE) # This is a consequence of sf now working with spherical geometry by default. 
    # Switch that off and sf will assume lat-long is planar and you get the result you expected:
    
    siofaGrid <- sf::st_intersection(SIGrid, siofa)
    
    taxa_data <- st_as_sf(taxa_data, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
    taxa_data <- taxa_data %>% st_join(siofaGrid) 
    taxa_data <- taxa_data[!is.na(taxa_data$cellid),]
    
    
    # Look for ObservationDate (i.e. inventory)
    taxa_data$date_richness <- as.Date(taxa_data$ObservationDate)
    
    # Calculate completeness
    ocurr_sample_events <- taxa_data[which(!is.na(taxa_data$ObservationDate)),]
    plyr::count(nchar(ocurr_sample_events$ObservationDate))
    index.obs.date <- which(nchar(ocurr_sample_events$ObservationDate)==0)
    if(length(index.obs.date)==0) {ocurr_sample_events=ocurr_sample_events} else
    {ocurr_sample_events <- ocurr_sample_events[-index.obs.date,]}
    plyr::count(nchar(ocurr_sample_events$ObservationDate))
    
    
    for(cell in unique(ocurr_sample_events$cellid))
    {
      current_cell_events <- ocurr_sample_events[which(ocurr_sample_events$cellid %in% cell), ]
      contin_cur_cell <- table(current_cell_events$species,
                               current_cell_events$ObservationDate)
      contin_cur_cell[contin_cur_cell > 0] <- 1
      
      #- Richness per cell
      richness_cur_cell <- colSums(contin_cur_cell)
      
      #- Occurrence for the current phylum
      occur_cur_cell <- rowSums(contin_cur_cell)
      
      #- Cell richness
      cell_richness <- length(occur_cur_cell[which(occur_cur_cell > 0)])
      
      #- Cell richness estimators
      chao_ind <- chao2(contin_cur_cell, taxa.row = TRUE)
      ice_ind <- ICE(contin_cur_cell, taxa.row = TRUE)
      jack1_ind <- jack1(contin_cur_cell, taxa.row = TRUE, abund = FALSE)
      
      #- Average estimator
      avg_ind <- round(mean(c(chao_ind, ice_ind, jack1_ind)))
      
      # completeness index
      completeness <- round(cell_richness/avg_ind, digits = 2)
      
      siofa_cell_completeness <- rbind.data.frame(siofa_cell_completeness,
                                                  data.frame(
                                                    Cell_id = cell,
                                                    Observed_richness = cell_richness,
                                                    chao2 = chao_ind,
                                                    ice = ice_ind,
                                                    jack1 = jack1_ind,
                                                    Average_estimated_richness = avg_ind,
                                                    Completeness_index = completeness,
                                                    cell_resolution = cellXsize[i],
                                                    taxon = "species")
      )
      
      # Check Observed richness distribution
      hist(siofa_cell_completeness$Observed_richness, breaks=120)
    }
  }
  
  
  # Plot
  
  C <- ggplot(siofa_cell_completeness, aes(x=as.factor(cell_resolution), y= Completeness_index)) + 
    geom_boxplot(fill="#ffffb3", alpha=0.8) + 
    xlab(expression("Cell resolution ("*~degree*")")) + ylab("Completeness") +
    theme_bw() +
    ggtitle('All cells')
  
  
  siofa_cell_completeness_more30 <- siofa_cell_completeness[siofa_cell_completeness$Observed_richness >= 30,]
  
  C_more30 <- ggplot(siofa_cell_completeness_more30, aes(x=as.factor(cell_resolution), y= Completeness_index)) + 
    geom_boxplot(fill="#ffffb3", alpha=0.8) + 
    xlab(expression("Cell resolution ("*~degree*")")) + ylab("Completeness") +
    theme_bw() +
    ggtitle('Cells with S >= 30')
  
  siofa_cell_completeness_less30 <- siofa_cell_completeness[siofa_cell_completeness$Observed_richness < 30,]
  
  C_less30 <- ggplot(siofa_cell_completeness_less30, aes(x=as.factor(cell_resolution), y= Completeness_index)) + 
    geom_boxplot(fill="#ffffb3", alpha=0.8) + 
    xlab(expression("Cell resolution ("*~degree*")")) + ylab("Completeness") +
    theme_bw() +
    ggtitle('Cells with S < 30')
  
  
  
  ## save plot
  library(cowplot)

  plot_grid(C, C_more30, C_less30, labels=c("A", "B", "C"), ncol = 3, nrow = 1)

  
  
  # Quick check of the percentage of cells with a desired 'X' completeness
  complete_freq <- siofa_cell_completeness[siofa_cell_completeness$cell_resolution== "1",]
  hist(complete_freq$Completeness_index)
  dim(complete_freq)
  p80 <- complete_freq[complete_freq$Completeness_index >= 0.8,]
  pless80 <- complete_freq[complete_freq$Completeness_index < 0.8,]
  p50 <- complete_freq[complete_freq$Completeness_index <= 0.5,]
  p20 <- complete_freq[complete_freq$Completeness_index <= 0.2,]
  
  perc.p80 <- nrow(p80)/nrow(complete_freq)*100
  perc.pless80 <- nrow(pless80)/nrow(complete_freq)*100
  perc.p50 <- nrow(p50)/nrow(complete_freq)*100
  perc.p20 <- nrow(p20)/nrow(complete_freq)*100
  
  
  