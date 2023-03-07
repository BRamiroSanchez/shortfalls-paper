#######################################################################
###                                                                 ### 
###     Script developed by Boris Leroy based on the working        ###
###     examples of the iNEXT.4steps R package (Chao & Hu, 2023)    ###
###     to calculate sample completeness profiles as a function of  ###
###     "q" diversity orders. 6 March 2023.                         ###
###                                                                 ### 
#######################################################################


library(raster) 
library(data.table)
library(ggplot2)
library(tidyr)
library(iNEXT.4steps)
library(sf)
library(dplyr)
library(egg)
library(grDevices)


# Data preparation --------------------------------------------------------

# Read in occurrence data
taylor <- read.csv("./data/Taylor_and_Rogers_2017.csv")
taylor <-  taylor[,c(1,2,5,8,11,12,13:15)]

taxa_list <- readRDS("./data/taxa_list.RDS") 
taxa_data <- droplevels(readRDS("./data/taxa_data.RDS")) 

taxon.level <- "species"
taxa_list <- taxa_list[taxa_list$type == taxon.level, ]
taxa_data <- taxa_data[which((taxa_data %>% pull(taxon.level)) %in% taxa_list$taxon), ] 
taxa_data <- taxa_data[,c(29, 7,8,9, 11, 16, 18, 20, 15)]

taxa_data <- rbind(taxa_data, taylor)


# Read in spatial data for context
siofa <- read_sf("./data/spatial/FAO_RFB_SIOFA/RFB_SIOFA.shp") %>% st_as_sf() 
land <- readRDS("./data/spatial/ne_10m_land/land")


# Graticule for maps
gr <- sf::st_graticule(lat = c(-89.9, seq(-80, 80, 20), 89.9),
                       lon = c(-179.9, seq(-160, 160, 20), 179.9))


region_limits=extent(20, 147, -65, 13) #xmin, xmax, ymin, ymax and #Lat = Y Long = X


# Grid cells at 0.01 resolution will be the basic sampling unit
raster_0.01=raster(res = 0.01, ext = region_limits,
                   crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" )
taxa_data$grid0.01 <- cellFromXY(raster_0.01,
                            taxa_data[,c("decimalLongitude", "decimalLatitude")])
rm(raster_0.01)


# Keep only those records that have an Observation Date
sample_events <- taxa_data[which(!is.na(taxa_data$ObservationDate)),]
plyr::count(nchar(sample_events$ObservationDate)) # check again whether there are records with no date

index.obs.date <- which(nchar(sample_events$ObservationDate)==0)
sample_events <- sample_events[-index.obs.date,]
plyr::count(nchar(sample_events$ObservationDate)) # check again it is okay


# Resolution loop ---------------------------------------------------------
sf_use_s2(FALSE)   ## needed because ggplot fails with geom_sf

global_occ <- list()
resolist <- 1:5

for(reso in resolist) {
  cat(paste0(Sys.time(), " --- starting resolution ", reso, "\n"))
  # Create rasters ----------------------------------------------------------

  cat(paste0("     ° Create raster\n"))
  
  raster_res <- raster(res = reso, ext = region_limits,
                       crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" )
  
  
  sample_events[, paste0("grid", reso)] <- cellFromXY(raster_res, # Recreating grid3 as well because of the swapped coordinates
                                                      sample_events[,c("decimalLongitude", "decimalLatitude")])
  
  
  
  
  # Global completeness profiles ------------------------------------
  
  contin <- table(sample_events$species, sample_events[, paste0("grid", reso)])
  contin[contin > 0] <- 1
  
  global_occ[[reso]] <- c(nT = ncol(contin), rowSums(contin))
  
  # Cellwise completeness profiles ------------------------------------
  
  # Create a list with occurrence frequencies within each cell
  occ_list <- list()
  length(occ_list) <- length(unique(sample_events[, paste0("grid", reso)]))
  names(occ_list) <- paste0("grid", unique(sample_events[, paste0("grid", reso)]))
  
  # Fill the list
  for (cell in unique(sample_events[, paste0("grid", reso)])) {
    suba <- sample_events[which(sample_events[, paste0("grid", reso)] == cell), ]
    nb_subcells <- length(unique(suba$grid0.01))
    
    if(any(duplicated(suba[, c("species", "grid0.01")]))) {
      suba <- suba[-which(duplicated(suba[, c("species", "grid0.01")])), ]
    }
    occ_suba <- plyr::count(suba$species)
    occ_suba2 <- occ_suba$freq
    names(occ_suba2) <- occ_suba$x
    
    
    occ_list[[paste0("grid", cell)]] <- c(nT = nb_subcells,
                                          occ_suba2)
  }
  
  
  cat(paste0("     ° Calculating completeness profiles...\n"))
  # Get the observed diversity at q = 0, 1 & 2 for each cell
  obs_s <- iNEXT.3D:::obsTD(occ_list, datatype = "incidence_freq", q = c(0, 1, 2),
                            nboot = 0,
                            conf = .95)
  
  # Prepare the data.frame to store important results
  results <- data.frame(cell = unique(sample_events[, paste0("grid", reso)]), # Grid cell number
                        cell_name = paste0("grid", unique(sample_events[, paste0("grid", reso)])),
                        nT = sapply(occ_list, function (x) x[1]), # Number of 0.01 cells with samplings
                        singleton_only = sapply(occ_list, function (x) all(x[-1] %in% 1)), # Are there only singletons in the samplings?
                        completeness_q0 = 0,
                        completeness_q1 = 0,
                        completeness_q2 = 0,
                        D0 = obs_s$qD[obs_s$Order.q == 0],
                        D1 = obs_s$qD[obs_s$Order.q == 1],
                        D2 = obs_s$qD[obs_s$Order.q == 2],
                        est_D0 = NA,
                        est_D1 = NA,
                        est_D2 = NA,
                        undetected_sp = NA,
                        undetected_frequent_sp = NA,
                        undetected_superfrequent_sp = NA,
                        standardized_coverage = 0)
  
  
  # Filter out inadequately sampled cells
  # Number of cells for which there are not enough samples (i.e., only 1 or 2 subcells, only singletons, or richness <= 10)
  results$undersampled_cells <- results$nT <=2 |
    results$singleton_only |
    results$D0 <= 10
  number_undersampled_cells <- length(results$undersampled_cells)
  if(number_undersampled_cells){
    occ_list_filtered <- occ_list[-which(results$undersampled_cells)]
  } else {
    occ_list_filtered <- occ_list
  }
  
  # Run the analysis
  res_4steps <- iNEXT4steps(occ_list_filtered,
                            q = c(0, 1, 2),
                            datatype = "incidence_freq",
                            details = FALSE) 

  saveRDS(res_4steps,
          paste0("./data/results_4steps_cellwise_res", reso, ".RDS"))
  
  # Fill values in the result table
  results$completeness_q0[match(res_4steps$summary$`STEP1. Sample completeness profiles`$Assemblage,
                                results$cell_name)] <-
    res_4steps$summary$`STEP1. Sample completeness profiles`$`q = 0`
  results$completeness_q1[match(res_4steps$summary$`STEP1. Sample completeness profiles`$Assemblage,
                                results$cell_name)] <-
    res_4steps$summary$`STEP1. Sample completeness profiles`$`q = 1`
  results$completeness_q2[match(res_4steps$summary$`STEP1. Sample completeness profiles`$Assemblage,
                                results$cell_name)] <-
    res_4steps$summary$`STEP1. Sample completeness profiles`$`q = 2`
  
  estimators <- res_4steps$summary$`STEP2. Asymptotic analysis`
  
  results$est_D0[match(estimators$Assemblage[which(estimators$Diversity == "Species richness")],
                       results$cell_name)] <-
    estimators$Estimator[which(estimators$Diversity == "Species richness")]
  
  results$est_D1[match(estimators$Assemblage[which(estimators$Diversity == "Shannon diversity")],
                       results$cell_name)] <-
    estimators$Estimator[which(estimators$Diversity == "Shannon diversity")]
  
  results$est_D2[match(estimators$Assemblage[which(estimators$Diversity == "Simpson diversity")],
                       results$cell_name)] <-
    estimators$Estimator[which(estimators$Diversity == "Simpson diversity")]
  
  results$undetected_sp <- results$est_D0 - results$D0
  results$undetected_frequent_sp <- results$est_D1 - results$D1
  results$undetected_superfrequent_sp <- results$est_D2 - results$D2
  
  std_cov <- sapply(occ_list_filtered, function(x) iNEXT.3D:::Coverage(x, "incidence_freq", 2*x[1]))
  results$standardized_coverage[match(names(std_cov),
                                      paste0(results$cell_name, ".nT"))] <- std_cov
  
  
  # Prepare figures -------------------------------------------------------------
  cat(paste0("     ° Prepare figures\n"))
  # Get centroid coordinates for all cells
  results[, c("X", "Y")] <- xyFromCell(raster_res,
                                       results$cell)
  
  
  
  # Completeness maps -------------------------------------------------------
  
  thememaps <- theme(legend.text = element_text(size = 14),
                     axis.text = element_text(size = 12),
                     axis.title = element_text(size = 14),
                     legend.title = element_text(size = 14),
                     panel.grid = element_blank(),
                     line = element_blank(),
                     #rect = element_blank(),
                     plot.margin = unit(c(2,0.1,1,0.1), "lines"))
  
  completeness_q0 <- ggplot() +
    geom_sf(data = gr, col = "lightgrey") +
    geom_sf(data = land) +
    geom_point(data = results,
               shape = 15, size = 0.5 * reso - 0.3, # NOTE: this formula for the size of points ensures no overlap
               # among cells at the current resolution AND for the chosen pdf size (see cairo_pdf below)
               aes(x = X, y = Y, col = completeness_q0)) +
    scale_color_viridis_b() +
    theme_bw() +
    thememaps +
    coord_sf(xlim= c(10, 155),
             ylim=c(-70, 25), 
             expand = FALSE) +
    xlab("Longitude") + ylab("Latitude")
  
  completeness_q1 <- ggplot() +
    geom_sf(data = gr, col = "lightgrey") +
    geom_sf(data = land) +
    geom_point(data = results, shape = 15, size = 0.5 * reso - 0.3,
               aes(x = X, y = Y, col = completeness_q1)) +
    scale_color_viridis_b() +
    theme_bw() +
    thememaps +
    coord_sf(xlim= c(10, 155),
             ylim=c(-70, 25), 
             expand = FALSE) +
    xlab("Longitude") + ylab("Latitude")
  
  completeness_q2 <- ggplot() +
    geom_sf(data = gr, col = "lightgrey") +
    geom_sf(data = land) +
    geom_point(data = results, shape = 15, size = 0.5 * reso - 0.3,
               aes(x = X, y = Y, col = completeness_q2)) +
    scale_color_viridis_b() +
    theme_bw() +
    thememaps +
    coord_sf(xlim= c(10, 155),
             ylim=c(-70, 25), 
             expand = FALSE) +
    xlab("Longitude") + ylab("Latitude")
  
  cairo_pdf(paste0("./outputs/completeness_maps_res", reso, ".pdf"),
            width = 9, height = 12)
  egg::ggarrange(completeness_q0, completeness_q1, completeness_q2,
                 labels = c("a. completeness at q = 0",
                            "b. completeness at q = 1",
                            "c. completeness at q = 2"),
                 newpage = FALSE,
                 nrow = 3)
  dev.off()

  
  # Coverage maps -----------------------------------------------------------
  
  
  gg_stdcov <- ggplot() +
    geom_sf(data = gr, col = "lightgrey") +
    geom_sf(data = land) +
    geom_point(data = results,
               shape = 15, size = 0.5 * reso - 0.3,
               aes(x = X, y = Y, col = standardized_coverage)) +
    scale_color_viridis_b(name = "Standardized\ncoverage") +
    theme_bw() +
    thememaps +
    coord_sf(xlim= c(10, 155),
             ylim=c(-70, 25), 
             expand = FALSE) +
    #coord_sf(expand = FALSE) +
    xlab("Longitude") + ylab("Latitude")
  
  
  cairo_pdf(paste0("./outputs/coverage_maps_res", reso, ".pdf"),
            width = 8, height = 4)
  print(gg_stdcov)
  dev.off()

  
  
  # Completeness profiles
  ggcomp <- pivot_longer(results,
                         cols = c("completeness_q0", "completeness_q1", "completeness_q2"))
  ggcomp$q <- as.numeric(as.factor(ggcomp$name)) - 1
  
  theme_bxp <- theme(legend.text = element_text(size = 14),
                     axis.text = element_text(size = 14),
                     axis.title = element_text(size = 14),
                     legend.title = element_text(size = 14),
                     panel.grid.minor = element_blank(),
                     plot.margin = unit(c(2,2.1,1,1.1), "lines"),
                     strip.text = element_text(size = 14))
  
  C_all <- ggplot(ggcomp, aes(x = q, y = value, group = q)) +
    geom_boxplot() + theme_bw() +
    scale_x_continuous(breaks = c(0, 1, 2)) +
    theme_bxp + ylab("Sample completeness")
  
  C_filtered <- ggplot(ggcomp[!ggcomp$undersampled_cells, ], aes(x = q, y = value, group = q)) +
    geom_boxplot() + theme_bw() +
    scale_x_continuous(breaks = c(0, 1, 2)) +
    theme_bxp + ylab("Sample completeness")
  
  bxp_std_cov <- ggplot(results[!results$undersampled_cells, ], aes(y = standardized_coverage)) +
    geom_boxplot()+ theme_bw()  +
    scale_x_continuous(breaks = NULL) +
    theme_bxp + ylab("Sample coverage")
  
  cairo_pdf(paste0("./outputs/boxplots_completeness_res", reso, ".pdf"),
            width = 12, height = 4)
  egg::ggarrange(C_all, C_filtered, bxp_std_cov,
                 labels = c("a. Completeness profiles (all cells)",
                            "b. Completeness profiles (evaluated cells only)",
                            "c. Standardized coverage\n(evaluated cells only)"),
                 newpage = FALSE,
                 ncol = 3, widths = c(1, 1, .45),
                 label.args = list(gp = grid::gpar(font = 4, cex = 1)))
  dev.off()
  
  # Undetected frequent species
  results$pc_und_q0 <- results$undetected_sp / results$est_D0
  results$pc_und_q1 <- results$undetected_frequent_sp / results$est_D1
  results$pc_und_q2 <- results$undetected_superfrequent_sp / results$est_D2
  
  ggund <- pivot_longer(results,
                        cols = c("pc_und_q0", "pc_und_q1", "pc_und_q2"))
  ggund$q <- as.numeric(as.factor(ggund$name)) - 1
  
  theme_bxp <- theme(legend.text = element_text(size = 14),
                     axis.text = element_text(size = 14),
                     axis.title = element_text(size = 14),
                     legend.title = element_text(size = 14),
                     panel.grid.minor = element_blank(),
                     plot.margin = unit(c(2,2.1,1,1.1), "lines"))
  
  ggund$q <- as.factor(ggund$name)
  levels(ggund$q) <- c("q = 0", "q = 1", "q = 2")
  
  # ggund$value2 <- ggund$value
  # ggund$value2[ggund$value2 > 100] <- 100
  undetected_sp <- ggplot(ggund, aes(y = value)) +
    geom_boxplot() + theme_bw() +
    facet_wrap(~ q) +
    scale_x_continuous(breaks = NULL)  +
    theme_bxp + ylab("Percentage of undetected diversity")
  

  cairo_pdf(paste0("./outputs/boxplots_undetected_res", reso, ".pdf"),
            width = 6, height = 4)
  print(undetected_sp)
  dev.off()
  
  # Save the output object as it will not be modified anymore
  saveRDS(results,
          paste0("./data/results_completeness_cellwise_res", reso, ".RDS"))
  
  
  cat(paste0(Sys.time(), " --- finished resolution ", reso, "\n"))
}


# Computing global completeness profiles ----------------------------------
saveRDS(global_occ,
        "./data/global_occ_list.RDS")

names(global_occ) <- paste0("global_resolution_", 1:5)

res <- iNEXT4steps(data = global_occ,
                   datatype = "incidence_freq")

res$summary$`STEP2. Asymptotic analysis`$Undetected <- 
  res$summary$`STEP2. Asymptotic analysis`$Estimator -
  res$summary$`STEP2. Asymptotic analysis`$Observed

saveRDS(res, paste0("./data/sample_completeness_global.RDS"))


cairo_pdf("./outputs/global_completeness_profiles.pdf", height = 8, width = 12)
print(res$figure)
dev.off()

C_all <- C_filtered <- undetected_sp <- list()
for(reso in resolist) {
  # Save the output object as it will not be modified anymore
  results <- readRDS(paste0("./data/results_completeness_cellwise_res", reso, ".RDS"))
  
  
  # Completeness profiles
  ggcomp <- pivot_longer(results,
                         cols = c("completeness_q0", "completeness_q1", "completeness_q2"))
  ggcomp$q <- as.numeric(as.factor(ggcomp$name)) - 1
  
  theme_bxp <- theme(legend.text = element_text(size = 14),
                     axis.text = element_text(size = 14),
                     axis.title = element_text(size = 14),
                     legend.title = element_text(size = 14),
                     panel.grid.minor = element_blank(),
                     plot.margin = unit(c(2,2.1,1,1.1), "lines"),
                     strip.text = element_text(size = 14))
  
  C_all[[reso]] <- ggplot(ggcomp, aes(x = q, y = value, group = q)) +
    geom_boxplot() + theme_bw() +
    scale_x_continuous(breaks = c(0, 1, 2)) +
    theme_bxp + ylab("Sample completeness")
  
  C_filtered[[reso]] <- ggplot(ggcomp[!ggcomp$undersampled_cells, ], aes(x = q, y = value, group = q)) +
    geom_boxplot() + theme_bw() +
    scale_x_continuous(breaks = c(0, 1, 2)) +
    theme_bxp + ylab("Sample completeness")
  
  
  
  
  # Undetected frequent species
  results$pc_und_q0 <- results$undetected_sp / results$est_D0
  results$pc_und_q1 <- results$undetected_frequent_sp / results$est_D1
  results$pc_und_q2 <- results$undetected_superfrequent_sp / results$est_D2
  
  ggund <- pivot_longer(results,
                        cols = c("pc_und_q0", "pc_und_q1", "pc_und_q2"))
  ggund$q <- as.numeric(as.factor(ggund$name)) - 1
  
  theme_bxp <- theme(legend.text = element_text(size = 14),
                     axis.text = element_text(size = 14),
                     axis.title = element_text(size = 14),
                     legend.title = element_text(size = 14),
                     panel.grid.minor = element_blank(),
                     plot.margin = unit(c(2,2.1,1,1.1), "lines"))
  
  ggund$q <- as.factor(ggund$name)
  levels(ggund$q) <- c("q = 0", "q = 1", "q = 2")
  
  # ggund$value2 <- ggund$value
  # ggund$value2[ggund$value2 > 100] <- 100
  undetected_sp[[reso]] <- ggplot(ggund, aes(y = value)) +
    geom_boxplot() + theme_bw() +
    facet_wrap(~ q) +
    scale_x_continuous(breaks = NULL)  +
    theme_bxp + ylab("Percentage of undetected diversity")
  
  cat("Resolution: ", reso, 
      " - Ncells = ", nrow(results), 
      " - N undersampled = ", length(which(results$undersampled_cells)),
      " - % = ", round(length(which(results$undersampled_cells)) /  nrow(results), 2),
      "\n")
  
}

# Resolution:  1  - Ncells =  481  - N undersampled =  392  - % =  0.81 
# Resolution:  2  - Ncells =  283  - N undersampled =  199  - % =  0.7 
# Resolution:  3  - Ncells =  208  - N undersampled =  133  - % =  0.64 
# Resolution:  4  - Ncells =  164  - N undersampled =  102  - % =  0.62 
# Resolution:  5  - Ncells =  130  - N undersampled =  73  - % =  0.56 

cairo_pdf(paste0("./outputs/boxplots_completeness.pdf"),
          width = 12, height = 20)
egg::ggarrange(C_all[[1]], C_filtered[[1]],
               C_all[[2]], C_filtered[[2]],
               C_all[[3]], C_filtered[[3]],
               C_all[[4]], C_filtered[[4]],
               C_all[[5]], C_filtered[[5]],
               
               labels = c("a. 1° - Completeness profiles (all cells)", 
                          "b. 1° - Completeness profiles (evaluated cells only)",
                          "c. 2° - Completeness profiles (all cells)", 
                          "d. 2° - Completeness profiles (evaluated cells only)",
                          "e. 3° - Completeness profiles (all cells)", 
                          "f. 3° - Completeness profiles (evaluated cells only)",
                          "g. 4° - Completeness profiles (all cells)", 
                          "h. 4° - Completeness profiles (evaluated cells only)",
                          "i. 5° - Completeness profiles (all cells)", 
                          "j. 5° - Completeness profiles (evaluated cells only)"),
               newpage = FALSE,
               ncol = 2, label.args = list(gp = grid::gpar(font = 4, cex = 1.5)))
dev.off()


cairo_pdf(paste0("./outputs/boxplots_undetected.pdf"),
          width = 6, height = 20)
egg::ggarrange(undetected_sp[[1]], 
               undetected_sp[[2]], 
               undetected_sp[[3]], 
               undetected_sp[[4]],
               undetected_sp[[5]], 
               
               labels = c("a. 1° - Undetected diversity (evaluated cells only)", 
                          "b. 2° - Undetected diversity (evaluated cells only)",
                          "c. 3° - Undetected diversity (evaluated cells only)", 
                          "d. 4° - Undetected diversity (evaluated cells only)",
                          "e. 5° - Undetected diversity (evaluated cells only)"),
               newpage = FALSE,
               ncol = 1, label.args = list(gp = grid::gpar(font = 4, cex = 1.5)))
dev.off()
  