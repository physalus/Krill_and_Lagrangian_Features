
#### Setup ####

## Load Packages
pkgTest <- function(x)
{
  if (!require(x,character.only = TRUE))
  {
    install.packages(x,dep=TRUE)
    if(!require(x,character.only = TRUE)) stop("Package not found")
  }
}
pkgTest("tidyverse")
pkgTest("stats")
pkgTest("viridis")
pkgTest("ncdf4")
pkgTest("marmap")
pkgTest("raster")
pkgTest("sf")
pkgTest("scales")
pkgTest("ggpubr") # for ggarrange
pkgTest("rnaturalearth") # for map base layers
pkgTest("paletteer")
pkgTest("metR") # for plotting bathy contours
pkgTest("ggspatial") # for plotting map annotations
pkgTest("ggnewscale")

## Load Data
load("./Figure_Files/dataGamm.RData") # Krill/Whale data output from analysis processing
load("./Figure_Files/dataCTD.RData") # CTD data output from analysis processing
load("./Figure_Files/Map_Data/bf_Sub.RData") # bathymetry subset for maps plots
load("./Figure_Files/Map_Data/HFRadarCoverageMCP.RData") # minimum convex polygon of HF Radar coverage
load("./Figure_Files/Map_Data/hfRadarStations.RData") # HF Radar station locations
load("./dataProcessed/acoustics_DF_FTLE_Line_w600.RData") # Full 200m X 5m Krill dataset 

## Global Variables
tzOffset <- "Etc/GMT+7" # all surveys are on the same timezone
FTLEdatapath <- "./dataRaw/FTLE/"
degreeBuffer <- .25 
quantilePlot <- seq(.3,.95,by=.07) 
plotLongMin <- -123.41
plotLongMax <- -122.45
plotLatMin <- 37.53
plotLatMax <- 38.04
markerSize <- 5.5
markerBump <- 0.021
# world data
world <- rnaturalearth::ne_countries(scale = 10, returnclass = 'sf')
bf_Sub1 <- bf_Sub %>% 
  dplyr::filter(x > plotLongMin-0.05,
                x < plotLongMax+0.05,
                y > plotLatMin-0.05,
                y < plotLatMax+0.05)

# Note that if you choose a Cruise-Line other than ACC1207-6, you need have the FTLE data downloaded.
c = "ACC1207" # Cruise
l = "6" # Line 
data_Sub <- dataGamm %>%
  dplyr::filter(Cruise == c, Line == l) 
cast_Sub <- dataCTD %>%
  dplyr::filter(Cruise == c, Line == l) 

# Subset full resolution acoustics data for this cruise/line
d_Sub_Line <- acoustics_DF_FTLE_W600 %>%
  dplyr::filter(Cruise == c, Line == l)
# Create a column summary for this line
lineSum <- d_Sub_Line %>% 
  mutate(minCID = min(as.integer(d_Sub_Line$ColumnID))) %>% 
  group_by(ColumnID) %>% 
  summarise(dttz = first(dttz),
            Lat = first(Lat),
            Long = first(Long),
            ftle24=first(ftle24_L),
            ftle48=first(ftle48_L),
            ftle120=first(ftle120_L),
            ftle240=first(ftle240_L),
            botDepth = first(BotDepth),
            KrillBiomass = sum(KrillBiomass),
            numWith = sum(KrillBiomass>0),
            density = if_else(numWith > 0, mean(KrillDensity[KrillDensity>0]),0),
            logDensity = if_else(density==0, -8,log10(density)),
            x = (as.integer(first(ColumnID))-first(minCID))*200)  

quantile(lineSum$logDensity)

# Load the FTLE data
if (file.exists(paste0("./Figure_Files/FTLE/",c,"_FTLE_24hrs.RData"))) {
  load(paste0("./Figure_Files/FTLE/",c,"_FTLE_24hrs.RData" ))
} else { 
  FTLE_stack24 <- stack(paste0(FTLEdatapath,c,"_FTLE_24Hrs.nc"), 
                      varname = "FTLE")
  save(FTLE_stack24,file=paste0("./Figure_Files/FTLE/",c,"_FTLE_24hrs.RData"))
}
if (file.exists(paste0("./Figure_Files/FTLE/",c,"_FTLE_48hrs.RData"))) {
  load(paste0("./Figure_Files/FTLE/",c,"_FTLE_48hrs.RData" ))
} else { 
  FTLE_stack48 <- stack(paste0(FTLEdatapath,c,"_FTLE_48Hrs.nc"), 
                        varname = "FTLE")
  save(FTLE_stack48,file=paste0("./Figure_Files/FTLE/",c,"_FTLE_48hrs.RData" ))
}
if (file.exists(paste0("./Figure_Files/FTLE/",c,"_FTLE_120hrs.RData"))) {
  load(paste0("./Figure_Files/FTLE/",c,"_FTLE_120hrs.RData" ))
} else { 
  FTLE_stack120 <- stack(paste0(FTLEdatapath,c,"_FTLE_120Hrs.nc"), 
                        varname = "FTLE")
  save(FTLE_stack120,file=paste0("./Figure_Files/FTLE/",c,"_FTLE_120hrs.RData" ))
}
if (file.exists(paste0("./Figure_Files/FTLE/",c,"_FTLE_240hrs.RData"))) {
  load(paste0("./Figure_Files/FTLE/",c,"_FTLE_240hrs.RData" ))
}else { 
  FTLE_stack240 <- stack(paste0(FTLEdatapath,c,"_FTLE_240Hrs.nc"), 
                        varname = "FTLE")
  save(FTLE_stack240,file=paste0("./Figure_Files/FTLE/",c,"_FTLE_240hrs.RData" ))
}

# FTLE PLOT
ftle_alpha <- function(n, flatten, scrunch) {
  c(scales::rescale(exp((1:flatten)/scrunch), to = c(0.25, 1)), rep(1, n - flatten))
}

#Determine Quantiles for FTLE Color Scaling
qFTLE24 <- quantile(FTLE_stack24, probs = quantilePlot,na.rm=TRUE)
qSumFTLE24 <- as.data.frame(summary(qFTLE24)) %>% 
  dplyr::select(-Var1) %>% 
  dplyr::filter(grepl("Mean",Freq)) %>%
  group_by(Var2) %>% 
  mutate(Means = unlist(str_split(Freq,":"))[2]) %>% 
  ungroup()
qFTLE48 <- quantile(FTLE_stack48, probs = quantilePlot,na.rm=TRUE)
qSumFTLE48 <- as.data.frame(summary(qFTLE48)) %>% 
  dplyr::select(-Var1) %>% 
  dplyr::filter(grepl("Mean",Freq)) %>%
  group_by(Var2) %>% 
  mutate(Means = unlist(str_split(Freq,":"))[2]) %>% 
  ungroup()
qFTLE120 <- quantile(FTLE_stack120, probs = quantilePlot,na.rm=TRUE)
qSumFTLE120 <- as.data.frame(summary(qFTLE120)) %>% 
  dplyr::select(-Var1) %>% 
  dplyr::filter(grepl("Mean",Freq)) %>%
  group_by(Var2) %>% 
  mutate(Means = unlist(str_split(Freq,":"))[2]) %>% 
  ungroup()
qFTLE240 <- quantile(FTLE_stack240,probs = quantilePlot,na.rm=TRUE)
qSumFTLE240 <- as.data.frame(summary(qFTLE240)) %>% 
  dplyr::select(-Var1) %>% 
  dplyr::filter(grepl("Mean",Freq)) %>%
  group_by(Var2) %>% 
  mutate(Means = unlist(str_split(Freq,":"))[2]) %>% 
  ungroup()

#### Maps ####

# create a dataframe of each integration
FTLE_spdf24 <- as(mean(FTLE_stack24[[first(d_Sub_Line$lineMinI):first(d_Sub_Line$lineMaxI)]]), "SpatialPixelsDataFrame")
FTLE_df24 <- as.data.frame(FTLE_spdf24)
colnames(FTLE_df24) <- c("FTLE", "x", "y")
FTLE_spdf48 <- as(mean(FTLE_stack48[[first(d_Sub_Line$lineMinI):first(d_Sub_Line$lineMaxI)]]), "SpatialPixelsDataFrame")
FTLE_df48 <- as.data.frame(FTLE_spdf48)
colnames(FTLE_df48) <- c("FTLE", "x", "y")
FTLE_spdf120 <- as(mean(FTLE_stack120[[first(d_Sub_Line$lineMinI):first(d_Sub_Line$lineMaxI)]]), "SpatialPixelsDataFrame")
FTLE_df120 <- as.data.frame(FTLE_spdf120)
colnames(FTLE_df120) <- c("FTLE", "x", "y")
FTLE_spdf240 <- as(mean(FTLE_stack240[[first(d_Sub_Line$lineMinI):first(d_Sub_Line$lineMaxI)]]), "SpatialPixelsDataFrame")
FTLE_df240 <- as.data.frame(FTLE_spdf240)
colnames(FTLE_df240) <- c("FTLE", "x", "y")

# 24 ####  
p24 <- ggplot() +
  geom_sf(data = world,color=NA,fill=NA) +
  geom_raster(data = FTLE_df24[complete.cases(FTLE_df24),], 
              mapping = aes(x, y, fill = FTLE)) +
  labs(x="", y="") +
  scale_fill_gradientn(colors = viridis(option = "plasma", n = 40, 
                                        alpha = ftle_alpha(40, 8, 8)),
                       limits=c(0, 1.25), 
                       breaks = c(0,0.6,1.2),
                       values = scales::rescale(
                         c(0,as.numeric(qSumFTLE24$Means[1:9]))),
                       oob=scales::squish,
                       na.value = NA,
                       name = expression(atop( bold('FTLE (24 hr)'), 
                                               bold(day^bold({bold("-1")})) )) ) +
  # add 100m and 200m contours
  ggplot2::geom_contour(data = bf_Sub1,
                        mapping = aes(x=x, y=y, z=z),
                        breaks=c(-100),
                        size=c(0.4),
                        colour="darkgray", show.legend = FALSE) +
  metR::geom_text_contour(data = bf_Sub1, mapping = aes(x=x,y=y,z=z),
                          breaks = c(-100), show.legend = FALSE, 
                          size = 2.5, alpha = 1, nudge_y = -.002,
                          color="darkgray") +
  ggplot2::geom_contour(data = bf_Sub1,
                        mapping = aes(x=x, y=y, z=z),
                        breaks=c(-200),
                        size=c(0.4),
                        colour="darkgray", show.legend = FALSE) +
  metR::geom_text_contour(data = bf_Sub1, mapping = aes(x=x,y=y,z=z),
                          breaks = c(-200), 
                          show.legend = FALSE, size = 2.5, alpha = 1,
                          nudge_y = -.002, color="darkgray") +
  geom_point(data = data_Sub, 
             mapping = aes(Long,Lat, color=log(KrillDensity_MeanNonZero) ),
             size=5,shape=15, show.legend = TRUE, 
             alpha = 1, fill="black") +
  geom_label(data = lineSum[c(1),],
             mapping = aes(Long-0.04,Lat),label=c("Start"),
             size=3.5)+
  geom_label(data = lineSum[c(length(lineSum$ColumnID) ),],
             mapping = aes(Long+0.04,Lat),label=c("End"),
             size=3.5)+
  scale_color_gradient(name = expression(atop(bold("Krill Density"),
                                              bold(log[e](g/m^3)) )), 
                       low="#FFCCFF",
                       high="magenta",
                       na.value = "white", 
                       limits = c(-0.5,3),
                       breaks=seq(0,3,1),
                       oob=scales::squish) +
  geom_point(data= data_Sub[data_Sub$HUWH == 1,], aes(y=Lat-markerBump,x=Long),color="black",
             fill="red",size=markerSize, shape=24)+
  geom_point(data= data_Sub[data_Sub$BLWH == 1,], aes(y=Lat-markerBump,x=Long),color="black", 
             fill="blue",size=markerSize,shape=24)+
   geom_point(data=cast_Sub, aes(y=Lat+markerBump,x=Long),color="black",fill="black",size=markerSize, shape=25)+
  geom_sf(data = world, size=1.25,fill="lightgrey",color="grey") +

  coord_sf(xlim = c(plotLongMin,
                    plotLongMax), 
           ylim = c(plotLatMin,
                    plotLatMax), 
           expand = FALSE) +
  theme_classic() +
  annotation_scale(location = "bl", width_hint = 0.25, text_face = "bold", text_col = "black") + 
  annotation_north_arrow(location = "bl", 
                         which_north = "true", 
                         pad_x = unit(0.25, "in"), 
                         pad_y = unit(0.25, "in"), 
                         style = north_arrow_fancy_orienteering(line_col = "black",
                                                                fill = c("white", "black"),
                                                                text_col = "black"  ) ) + 
  xlab("Longitude") + 
  ylab("Latitude") + 
  guides(color = guide_colourbar(order=1,direction = 'horizontal' ,
                                 title.position = "top", 
                                 title.hjust = .5,
                                 label.position = "bottom"), 
         fill = guide_colourbar(order=2,direction = 'horizontal', 
                                title.position = "top", 
                                title.hjust = .5,
                                label.position = "bottom")) +
  theme(panel.grid.major = element_line(color = gray(.5), 
                                        linetype = "dashed", size = 0.5),
        panel.background = element_rect(fill = "aliceblue"),
        axis.text = element_text(face="bold", size=12),
        axis.title = element_blank(),
        legend.text = element_text(face="bold", size=10),
        legend.position = c(.125,.725),
        legend.box = 'vertical',
        legend.title = element_text(face="bold", size=12),
        legend.box.background = element_rect(color="white",
                                             size=2.5, fill="white"),
        legend.spacing.y = unit(0.025, "cm"),
        # adds a buffer on the sides so the plot is centered
        plot.margin=margin(t = 0.5, r = 0.25, b = 0, l = 0.25, "cm"), 
        plot.title = element_text(face="bold", size=24, hjust = 0.5),
        plot.subtitle = element_text(face="bold", size=20, hjust = 0.5))

# 48 ####
p48 <- ggplot()+
  geom_sf(data = world,color=NA,fill=NA) +
  geom_raster(data = FTLE_df48[complete.cases(FTLE_df48),], 
              mapping = aes(x, y, fill = FTLE)) +
  labs(x="", y="") +
  scale_fill_gradientn(colors = viridis(option = "plasma", n = 40, 
                                        alpha = ftle_alpha(40, 8, 8)),
                       values = scales::rescale(
                         c(0,as.numeric(qSumFTLE48$Means[1:9]))),
                       limits=c(0, 1), # for >= 24hr
                       breaks = c(0,0.5,1),
                       oob=scales::squish,
                       na.value = NA,
                       name = expression(atop( bold('FTLE (48 hr)'), 
                                               bold(day^bold({bold("-1")})) )) ) +
  # add 100m and 200m contours
  ggplot2::geom_contour(data = bf_Sub1,
                        mapping = aes(x=x, y=y, z=z),
                        breaks=c(-100),
                        size=c(0.4),
                        colour="darkgray", show.legend = FALSE) +
  metR::geom_text_contour(data = bf_Sub1, mapping = aes(x=x,y=y,z=z),
                          breaks = c(-100), show.legend = FALSE, 
                          size = 2.5, alpha = 1, nudge_y = -.002,
                          color="darkgray") +
  ggplot2::geom_contour(data = bf_Sub1,
                        mapping = aes(x=x, y=y, z=z),
                        breaks=c(-200),
                        size=c(0.4),
                        colour="darkgray", show.legend = FALSE) +
  metR::geom_text_contour(data = bf_Sub1, mapping = aes(x=x,y=y,z=z),
                          breaks = c(-200), 
                          show.legend = FALSE, size = 2.5, alpha = 1,
                          nudge_y = -.002, color="darkgray") +
  geom_point(data = data_Sub, 
             mapping = aes(Long,Lat, color=log(KrillDensity_MeanNonZero) ),
             
             size=5,shape=15, show.legend = TRUE, 
             alpha = 1, fill="black") +
  geom_label(data = lineSum[c(1),],
             mapping = aes(Long-0.04,Lat),label=c("Start"),
             size=3.5)+
  geom_label(data = lineSum[c(length(lineSum$ColumnID) ),],
             mapping = aes(Long+0.04,Lat),label=c("End"),
             size=3.5)+
  scale_color_gradient(name = expression(atop(bold("Krill Density"),
                                              bold(log[e](g/m^3)) )),
                       low="#FFCCFF",
                       high="magenta",
                       na.value = "white", 
                       limits = c(-0.5,3),
                       breaks=seq(0,3,1),
                       oob=scales::squish) +
  
  geom_point(data= data_Sub[data_Sub$HUWH == 1,], aes(y=Lat-markerBump,x=Long),color="black",
             fill="red",size=markerSize, shape=24)+
  geom_point(data= data_Sub[data_Sub$BLWH == 1,], aes(y=Lat-markerBump,x=Long),color="black", 
             fill="blue",size=markerSize,shape=24)+
  geom_point(data=cast_Sub, aes(y=Lat+markerBump,x=Long),color="black",
             fill="black",size=markerSize, shape=25)+
  geom_sf(data = world, size=1.25,fill="lightgrey",color="grey") +
  coord_sf(xlim = c(plotLongMin,
                    plotLongMax), 
           ylim = c(plotLatMin,
                    plotLatMax), 
           expand = FALSE) +
  theme_classic() +
  annotation_scale(location = "bl", width_hint = 0.25, text_face = "bold", text_col = "black") + 
  annotation_north_arrow(location = "bl", 
                         which_north = "true", 
                         pad_x = unit(0.25, "in"), 
                         pad_y = unit(0.25, "in"), 
                         style = north_arrow_fancy_orienteering(line_col = "black",
                                                                fill = c("white", "black"),
                                                                text_col = "black"  ) ) + 
  xlab("Longitude") + 
  ylab("Latitude") + 
  guides(color = guide_colourbar(order=1,direction = 'horizontal' ,
                                 title.position = "top", 
                                 title.hjust = .5,
                                 label.position = "bottom"), 
         fill = guide_colourbar(order=2,direction = 'horizontal', 
                                title.position = "top", 
                                title.hjust = .5,
                                label.position = "bottom")) +
  theme(panel.grid.major = element_line(color = gray(.5), 
                                        linetype = "dashed", size = 0.5),
        panel.background = element_rect(fill = "aliceblue"),
        axis.text = element_text(face="bold", size=12),
        axis.title = element_blank(),
        legend.text = element_text(face="bold", size=10),
        legend.position = c(.125,.725),
        legend.box = 'vertical',
        legend.title = element_text(face="bold", size=12),
        legend.box.background = element_rect(color="white",
                                             size=2.5, fill="white"),
        legend.spacing.y = unit(0.025, "cm"),
        # adds a buffer on the sides so the plot is centered
        plot.margin=margin(t = 0.5, r = 0.25, b = 0, l = 0.25, "cm"), 
        plot.title = element_text(face="bold", size=24, hjust = 0.5),
        plot.subtitle = element_text(face="bold", size=20, hjust = 0.5))

# 120 ####
p120 <- ggplot()+
  geom_sf(data = world,color=NA,fill=NA) +
  geom_raster(data = FTLE_df120[complete.cases(FTLE_df120),], 
              mapping = aes(x, y, fill = FTLE)) +
  labs(x="", y="") +
  scale_fill_gradientn(colors = viridis(option = "plasma", n = 40, 
                                        alpha = ftle_alpha(40, 8, 8)),
                       values = scales::rescale(
                         c(0,as.numeric(qSumFTLE120$Means[1:9]))),
                       limits=c(0, .8), 
                       breaks = c(0,.4,.8),
                       oob=scales::squish,
                       na.value = NA,
                       name = expression(atop( bold('FTLE (120 hr)'), 
                                               bold(day^bold({bold("-1")})) )) ) +
  # add 100m and 200m contours
  ggplot2::geom_contour(data = bf_Sub1,
                        mapping = aes(x=x, y=y, z=z),
                        breaks=c(-100),
                        size=c(0.4),
                        colour="darkgray", show.legend = FALSE) +
  metR::geom_text_contour(data = bf_Sub1, mapping = aes(x=x,y=y,z=z),
                          breaks = c(-100), show.legend = FALSE, 
                          size = 2.5, alpha = 1, nudge_y = -.002,
                          color="darkgray") +
  ggplot2::geom_contour(data = bf_Sub1,
                        mapping = aes(x=x, y=y, z=z),
                        breaks=c(-200),
                        size=c(0.4),
                        colour="darkgray", show.legend = FALSE) +
  metR::geom_text_contour(data = bf_Sub1, mapping = aes(x=x,y=y,z=z),
                          breaks = c(-200), 
                          show.legend = FALSE, size = 2.5, alpha = 1,
                          nudge_y = -.002, color="darkgray") +
  geom_point(data = data_Sub, 
             mapping = aes(Long,Lat, color=log(KrillDensity_MeanNonZero) ),
             size=5,shape=15, show.legend = TRUE, 
             alpha = 1, fill="black") +
  geom_label(data = lineSum[c(1),],
             mapping = aes(Long-0.04,Lat),label=c("Start"),
             size=3.5)+
  geom_label(data = lineSum[c(length(lineSum$ColumnID) ),],
             mapping = aes(Long+0.04,Lat),label=c("End"),
             size=3.5)+
  scale_color_gradient(name = expression(atop(bold("Krill Density"),
                                              bold(log[e](g/m^3)) )),
                       low="#FFCCFF",
                       high="magenta",
                       na.value = "white", 
                       limits = c(-0.5,3),
                       breaks=seq(0,3,1),
                       oob=scales::squish) +
  geom_point(data= data_Sub[data_Sub$HUWH == 1,], aes(y=Lat-markerBump,x=Long),color="black",
             fill="red",size=markerSize, shape=24)+
  geom_point(data= data_Sub[data_Sub$BLWH == 1,], aes(y=Lat-markerBump,x=Long),color="black", 
             fill="blue",size=markerSize,shape=24)+
  geom_point(data=cast_Sub, aes(y=Lat+markerBump,x=Long),color="black",
             fill="black",size=markerSize, shape=25)+
  geom_sf(data = world, size=1.25,fill="lightgrey",color="grey") +
  coord_sf(xlim = c(plotLongMin,
                    plotLongMax), 
           ylim = c(plotLatMin,
                    plotLatMax),
           expand = FALSE) +
  theme_classic() +
  annotation_scale(location = "bl", width_hint = 0.25, text_face = "bold", text_col = "black") + 
  annotation_north_arrow(location = "bl", 
                         which_north = "true", 
                         pad_x = unit(0.25, "in"), 
                         pad_y = unit(0.25, "in"), 
                         style = north_arrow_fancy_orienteering(line_col = "black",
                                                                fill = c("white", "black"),
                                                                text_col = "black"  ) ) + 
  xlab("Longitude") + 
  ylab("Latitude") +  
  guides(color = guide_colourbar(order=1,direction = 'horizontal' ,
                                 title.position = "top", 
                                 title.hjust = .5,
                                 label.position = "bottom"), 
         fill = guide_colourbar(order=2,direction = 'horizontal', 
                                title.position = "top", 
                                title.hjust = .5,
                                label.position = "bottom")) +
  theme(panel.grid.major = element_line(color = gray(.5), 
                                        linetype = "dashed", size = 0.5),
        panel.background = element_rect(fill = "aliceblue"),
        axis.text = element_text(face="bold", size=12),
        axis.title = element_blank(),
        legend.text = element_text(face="bold", size=10),
        legend.position = c(.125,.725),
        legend.box = 'vertical',
        legend.title = element_text(face="bold", size=12),
        legend.box.background = element_rect(color="white",
                                             size=2.5, fill="white"),
        legend.spacing.y = unit(0.025, "cm"),
        # adds a buffer on the sides so the plot is centered
        plot.margin=margin(t = 0.5, r = 0.25, b = 0, l = 0.25, "cm"), 
        plot.title = element_text(face="bold", size=24, hjust = 0.5),
        plot.subtitle = element_text(face="bold", size=20, hjust = 0.5))
# p120

# 240 ####
p240 <- ggplot()+
  geom_sf(data = world,color=NA,fill=NA) +
  geom_raster(data = FTLE_df240[complete.cases(FTLE_df240),], 
              mapping = aes(x, y, fill = FTLE)) +
  labs(x="", y="") +
  scale_fill_gradientn(colors = viridis(option = "plasma", n = 40, 
                                        alpha = ftle_alpha(40, 8, 8)),
                       values = scales::rescale(
                         c(0,as.numeric(qSumFTLE240$Means[1:9]))),
                       limits=c(0, .6), 
                       breaks = c(0,.3,.6),
                       oob=scales::squish,
                       na.value = NA,
                       name = expression(atop( bold('FTLE (240 hr)'), 
                                               bold(day^bold({bold("-1")})) )) ) +
  # add 100m and 200m contours
  ggplot2::geom_contour(data = bf_Sub1,
                        mapping = aes(x=x, y=y, z=z),
                        breaks=c(-100),
                        size=c(0.4),
                        colour="darkgray", show.legend = FALSE) +
  metR::geom_text_contour(data = bf_Sub1, mapping = aes(x=x,y=y,z=z),
                          breaks = c(-100), show.legend = FALSE, 
                          size = 2.5, alpha = 1, nudge_y = -.002,
                          color="darkgray") +
  ggplot2::geom_contour(data = bf_Sub1,
                        mapping = aes(x=x, y=y, z=z),
                        breaks=c(-200),
                        size=c(0.4),
                        colour="darkgray", show.legend = FALSE) +
  metR::geom_text_contour(data = bf_Sub1, mapping = aes(x=x,y=y,z=z),
                          breaks = c(-200), 
                          show.legend = FALSE, size = 2.5, alpha = 1,
                          nudge_y = -.002, color="darkgray") +
  geom_point(data = data_Sub, 
             mapping = aes(Long,Lat, color=log(KrillDensity_MeanNonZero) ),
             size=5,shape=15, show.legend = TRUE, 
             alpha = 1, fill="black") +
  geom_label(data = lineSum[c(1),],
             mapping = aes(Long-0.04,Lat),label=c("Start"),
             size=3.5)+
  geom_label(data = lineSum[c(length(lineSum$ColumnID) ),],
             mapping = aes(Long+0.04,Lat),label=c("End"),
             size=3.5)+
  scale_color_gradient(name = expression(atop(bold("Krill Density"),
                                              bold(log[e](g/m^3)) )),
                       low="#FFCCFF",
                       high="magenta",
                       na.value = "white", 
                       limits = c(-0.5,3),
                       breaks=seq(0,3,1),
                       oob=scales::squish) +
  geom_point(data= data_Sub[data_Sub$HUWH == 1,], aes(y=Lat-markerBump,x=Long),color="black",
             fill="red",size=markerSize, shape=24)+
  geom_point(data= data_Sub[data_Sub$BLWH == 1,], aes(y=Lat-markerBump,x=Long),color="black", 
             fill="blue",size=markerSize,shape=24)+
  geom_point(data=cast_Sub, aes(y=Lat+markerBump,x=Long),color="black",
             fill="black",size=markerSize, shape=25)+
  geom_sf(data = world, size=1.25,fill="lightgrey",color="grey") +
  coord_sf(xlim = c(plotLongMin,
                    plotLongMax), 
           ylim = c(plotLatMin,
                    plotLatMax), 
           expand = FALSE) +
  theme_classic() +
  annotation_scale(location = "bl", width_hint = 0.25, text_face = "bold", text_col = "black") + 
  annotation_north_arrow(location = "bl", 
                         which_north = "true", 
                         pad_x = unit(0.25, "in"), 
                         pad_y = unit(0.25, "in"), 
                         style = north_arrow_fancy_orienteering(line_col = "black",
                                                                fill = c("white", "black"),
                                                                text_col = "black"  ) ) + 
  xlab("Longitude") + 
  ylab("Latitude") + 
  guides(color = guide_colourbar(order=1,direction = 'horizontal' ,
                                 title.position = "top", 
                                 title.hjust = .5,
                                 label.position = "bottom"), 
         fill = guide_colourbar(order=2,direction = 'horizontal', 
                                title.position = "top", 
                                title.hjust = .5,
                                label.position = "bottom")) +
  theme(panel.grid.major = element_line(color = gray(.5), 
                                        linetype = "dashed", size = 0.5),
        panel.background = element_rect(fill = "aliceblue"),
        axis.text = element_text(face="bold", size=12),
        axis.title = element_blank(),
        legend.text = element_text(face="bold", size=10),
        legend.position = c(.125,.725),
        legend.box = 'vertical',
        legend.title = element_text(face="bold", size=12),
        legend.box.background = element_rect(color="white",
                                             size=2.5, fill="white"),
        legend.spacing.y = unit(0.025, "cm"),
        # adds a buffer on the sides so the plot is centered
        plot.margin=margin(t = 0.5, r = 0.25, b = 0, l = 0.25, "cm"), 
        plot.title = element_text(face="bold", size=24, hjust = 0.5),
        plot.subtitle = element_text(face="bold", size=20, hjust = 0.5))
# p240

# Plot All ####
pMapH <- ggarrange(p24,p48,p120,p240, 
                   ncol=2, nrow = 2, 
                   align = "hv")
pMapH

ggsave(sprintf("./Output/%s-%s_FTLEmapWithLegend.png", c,l ), width=14, height=9.22, units = "in", dpi=600) # No gaps



# Acoustics Panel ####

# Create a matrix of the vertical 
x = matrix(nrow = (2+length(unique(as.integer(d_Sub_Line$CellID)))), 
           ncol = length(unique(as.integer(d_Sub_Line$ColumnID))),
           byrow=FALSE)
for(ii in 1:length(unique(as.integer(d_Sub_Line$ColumnID)))){
  numNA <- nrow(x)- (2+max(as.integer(d_Sub_Line$CellID[d_Sub_Line$ColumnID==unique(as.integer(d_Sub_Line$ColumnID))[ii]])))
  if(numNA>0){
    x[,ii] <- c(NA, NA, d_Sub_Line$KrillDensity[d_Sub_Line$ColumnID==unique(as.integer(d_Sub_Line$ColumnID))[ii]], rep(NA, numNA))
  } else {
    x[,ii] <- c(NA, NA, d_Sub_Line$KrillDensity[d_Sub_Line$ColumnID==unique(as.integer(d_Sub_Line$ColumnID))[ii]])
  }
}
# Create a raster from transect values
r <- raster(ncol=ncol(x),nrow=nrow(x))
raster::values(r) <- x
bb <- extent(0, (ncol(x)-1)*200,(nrow(x)-1)*-5, 0)
extent(r) <- bb
r <- setExtent(r, bb, keepres=TRUE)

rSPDF <- as(r, "SpatialPixelsDataFrame")
r_df <- as.data.frame(rSPDF)
colnames(r_df) <- c("KrillDensity", "x", "y")
r_df <- arrange(r_df, x, -y)
# add x and y back onto d_sub_Line
d_Sub_Line <- cbind(d_Sub_Line,x=r_df$x,y=r_df$y)

d_sub600 <- d_Sub_Line %>% 
  group_by(Col600) %>% 
  summarize(minX = min(x)-100,
            maxX= max(x)+100,
            minY = min(y)+5, # adjust to represent the top of the column
            maxY= 0,
            Krill = if_else(sum(KrillBiomass)>0,1,0)) 

cast_SubSum <- cast_Sub %>% 
  mutate(x = c(0,4180,7763,28661))

whaleCTD <- d_sub600 %>% 
  left_join(dplyr::select(data_Sub,Col600, HUWH, BLWH), by = "Col600" )


ylim.prim <- c(round(min(r_df$y))-10, 0)
ylim.sec <- c(1000,0)   

b <- diff(ylim.prim)/diff(ylim.sec)
a <- ylim.prim[1] - b*ylim.sec[1] 

# Inset location
colMin <- 7
colMax <- 10
minDepth <- -160
maxDepth <- -200

# Main Plot

pKrillPA<- ggplot() +  
  geom_rect(data=d_sub600, aes(xmin = minX, xmax = maxX, ymin = minY-2.5, ymax = -5,
                               alpha="600m Columns"),color="black", fill="white") +
  geom_rect(data=d_sub600, aes(xmin = minX, xmax = maxX, ymin = minY-2.5, ymax = -5,
                               fill=as.factor(Krill)),color="black") +
  scale_fill_manual(name = NULL, values = c("1"="#b3ffff","0"="white"),
                    labels = c("Krill Present in Column", "No Krill in Column"),
                    guide = guide_legend(order=3,direction = 'vertical') 
                    )+
  new_scale_fill() +   
  geom_tile(data=r_df,aes(x=x, y=y+5,fill=log(KrillDensity)),color="darkgrey", alpha=0.95) + 
  scale_fill_gradient(name = expression(atop(bold("Krill Density"),
                                             bold(log[e](g/m^3)) )), 
                      low="#FFCCFF",
                      high="magenta",
                      na.value = NA, 
                      limits = c(-0.5,3),
                      breaks=seq(0,3,1),
                      oob=scales::squish) +
  geom_rect(data=d_sub600, aes(xmin = minX, xmax = maxX, ymin = minY-2.5, ymax = -5,
                               alpha="600m Columns"),color="black", fill=NA, size = 0.75) +
  scale_alpha_manual(name=NULL, values=c("600m Columns"=0.1)) +
  geom_point(data= whaleCTD[whaleCTD$HUWH == 1,], aes(y=1.5,x=minX+300, shape = "Humpback Whale Sighting"),
             fill = "red", color="black", size=markerSize)+
  geom_point(data= whaleCTD[whaleCTD$BLWH == 1,], aes(y=1.5,x=minX+300, shape = "Blue Whale Sighting"),
             fill = "blue", color="black", size=markerSize)+
  geom_point(data=cast_SubSum, aes(y=6,x=x, shape = "CTD Profile"),
             fill = "black", color="black", size=markerSize)+
  scale_shape_manual(name = NULL, values = c("Humpback Whale Sighting" = 24,
                                          "Blue Whale Sighting" = 24,
                                          "CTD Profile" = 25) )+
  scale_y_continuous(limits=c(round(min(r_df$y))-5, 6), expand = c(0.025,0))+
  scale_x_continuous(expand = c(0.0115,0.00005),limits=c(-1,(max(r_df$x))-800 ),
                     breaks = c(10000,20000,30000,40000),
                     labels = c("10","20","30","40"))+
  xlab("Transect Distance (km)")+ ylab("Depth (m)")+
  guides(
         fill = guide_colourbar(order=1,direction = 'horizontal',title.position = "left", 
                                title.hjust = .5,label.position = "bottom"),
         alpha = guide_legend(order=2,direction = 'vertical'
                              ),
         shape = guide_legend(order=4,direction = 'vertical',
                              override.aes = list(fill = c("red", "blue", "black"))  )) +
  theme_classic()+
  theme(legend.position= c(.85,.375),
        legend.direction = "horizontal",
        legend.title = element_text(face="bold", size=18),
        legend.text = element_text(face="bold", size=16),
        legend.key.size = unit(1, "cm"),
        legend.spacing.y= unit(0.15, "cm"),
        plot.title = element_text(face="bold", size=24,hjust = 0.5),
        plot.subtitle = element_text( face="bold", size=22,hjust = 0.5),
        plot.caption = element_text( face="bold", size=18,hjust = 0),
        axis.title = element_text(face="bold", size=18),
        axis.text = element_text(face="bold", size=16),
        axis.title.x = element_text(face="bold", size=18),
        axis.text.x = element_text(face="bold", size=16))
pKrillPA

ggsave(sprintf("./Output/%s-%s_KrillColumnsPA.jpg", c,l ), width=14, height=7, units = "in", dpi=600) # No gaps



pKrillPAinset <- ggplot() +  
  geom_rect(data=d_sub600, aes(xmin = minX, xmax = maxX, ymin = maxDepth, ymax = minDepth,
                               alpha="600m Columns"),color="black", fill="white") +
  geom_rect(data=d_sub600, aes(xmin = minX, xmax = maxX,ymin = maxDepth, ymax = minDepth,
                               fill=as.factor(Krill)),color="black") +
  scale_fill_manual(name = "Krill Presence", values = c("1"="#b3ffff","0"="white") ,guide="none" )+
  new_scale_fill() +   
  geom_tile(data=r_df,aes(x=x, y=y+1.25,fill=log(KrillDensity)),color="darkgrey", alpha=1) + 
  geom_rect(data=d_sub600, aes(xmin = minX, xmax = maxX, ymin = maxDepth, ymax = minDepth,
                               alpha="600m Columns"),color="black", fill=NA, size = 0.75) +
  scale_y_continuous(limits=c(maxDepth-2, minDepth+2), expand = c(0,0))+
  scale_x_continuous(limits=c(d_sub600$minX[d_sub600$Col600==colMin],  d_sub600$maxX[d_sub600$Col600==colMax]),
                     expand = c(0.05,0))+
  scale_alpha_manual(name=NULL, values=c("600m Columns"=0.1)) +
  scale_fill_gradient(name = expression(atop(bold("Krill Density"),
                                             bold(log[e](g/m^3)) )), 
                      low="#FFCCFF",
                      high="magenta",
                      na.value = NA, 
                      limits = c(-0.5,3),
                      breaks=seq(0,3,1),
                      oob=scales::squish) +
  xlab("Transect Distance (m)")+ ylab("Depth (m)")+
  theme_classic()+
  theme(legend.position= 'none',
        legend.direction = "horizontal",
        legend.title = element_text(face="bold", size=22),
        legend.text = element_text(face="bold", size=22),
        legend.key.size = unit(2, "cm"),
        plot.title = element_text(face="bold", size=24,hjust = 0.5),
        plot.subtitle = element_text( face="bold", size=22,hjust = 0.5),
        plot.caption = element_text( face="bold", size=18,hjust = 0),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank())
pKrillPAinset
ggsave(sprintf("./Output/%s-%s_KrillColumnsPAinset.jpg", c,l ), width=8, height=8, units = "in", dpi=600) # No gaps


