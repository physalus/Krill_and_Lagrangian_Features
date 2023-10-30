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
pkgTest("marmap")
pkgTest("sf")
pkgTest("rnaturalearth") # for map base layers
pkgTest("paletteer")
pkgTest("metR") # for plotting bathy contours

## Load Data
load("./Figure_Files/Map_Data/bf.RData") # bathymetry subset for maps plots
load("./Figure_Files/dataGamm.RData") # Krill/Whale data output from analysis processing
load("./Figure_Files/dataCTD.RData") # CTD data output from analysis processing
load("./Figure_Files/Map_Data/HFRadarCoverageMCP.RData") # minimum convex polygon of HF Radar coverage
load("./Figure_Files/Map_Data/hfRadarStations.RData") # HF Radar station locations

# world data
world <- rnaturalearth::ne_countries(scale = 10, returnclass = 'sf')

## Global Variables
plotLongMin <- -123.41
plotLongMax <- -122.45
plotLatMin <- 37.53
plotLatMax <- 38.04
# bathymetry subsets are necessary here to get the labels to appear correctly
bf_Sub2 <- bf %>% 
  dplyr::filter(x > plotLongMin-0.05,
                x < plotLongMax+0.05,
                y > plotLatMin-0.05,
                y < plotLatMax+0.05)
plotLatMin2 <- 38.05
plotLatMax2 <- 38.85
bf_Sub3 <- bf %>% 
  dplyr::filter(x > plotLongMin-0.25,
                x < plotLongMax+0.25,
                y > plotLatMin2,
                y < plotLatMax2)

#### Main Plot ####
# All Lines Together
p_AllLines <-  dataGamm %>% 
dplyr::select(Line, Lat,Long) %>% 
  unique() %>%   
  ggplot() +
  geom_sf(data = world,color=NA,fill=NA) +
  # add bathymetry contours
  geom_contour(data = bf,
               aes(x=x, y=y, z=z),
               breaks=c(-100,-200,-400,-800,-1600),
               size=c(0.4),
               colour="darkgrey", show.legend = FALSE) +
  metR::geom_text_contour(data = bf_Sub2, aes(x=x, y=y,z = z),breaks=c(-100),
                          show.legend = FALSE, size = 4, alpha = .6, nudge_y = -.002, check_overlap=TRUE,
                          label.placer = label_placer_random(seed=41,rot_adjuster = isoband::angle_halfcircle_bottom())) +
  metR::geom_text_contour(data = bf_Sub3, aes(x=x, y=y,z = z),breaks=c(-100),
                          show.legend = FALSE, size = 4, alpha = .6, nudge_y = -.002, check_overlap=TRUE,
                          label.placer = label_placer_random(seed=41,rot_adjuster = isoband::angle_halfcircle_bottom())) +
  metR::geom_text_contour(data = bf_Sub3, aes(x=x, y=y,z = z),breaks=c(-200),
                          show.legend = FALSE, size = 4, alpha = .6, nudge_y = -.002, check_overlap=TRUE,
                          label.placer = label_placer_random(seed=75,rot_adjuster = isoband::angle_halfcircle_bottom())) +
  metR::geom_text_contour(data = bf_Sub3, aes(x=x, y=y,z = z),breaks=c(-400),
                          show.legend = FALSE, size = 4, alpha = .6, nudge_y = -.002, check_overlap=TRUE,
                          label.placer = label_placer_random(seed=88,rot_adjuster = isoband::angle_halfcircle_bottom())) +
  metR::geom_text_contour(data = bf_Sub3, aes(x=x, y=y,z = z),breaks=c(-800),
                          show.legend = FALSE, size = 4, alpha = .6, nudge_y = -.002, check_overlap=TRUE,
                          label.placer = label_placer_random(seed=94,rot_adjuster = isoband::angle_halfcircle_bottom())) +
  metR::geom_text_contour(data = bf_Sub3, aes(x=x, y=y,z = z),breaks=c(-1600),
                          show.legend = FALSE, size = 4, alpha = .6, nudge_y = -.002, check_overlap=TRUE,
                          label.placer = label_placer_random(seed=96,rot_adjuster = isoband::angle_halfcircle_bottom())) +
  geom_point(aes(Long,Lat,color="ACCESS Krill Column"),alpha=.8,stroke=0.9,size=2, shape=21,show.legend = TRUE) +
  geom_sf(data = world, size=.75,fill="grey") +
  geom_point(data=dataCTD, aes(x=Long, y=Lat, alpha="CTD Cast"),size=5)+
  geom_point(data=hfRadar_DF,
             aes(x=Long, y=Lat, shape = ID),color="black",size=4,stroke=1.5,
             fill = 'white') +
  scale_color_manual(values = c("ACCESS Krill Column"="#046C9A")) +
  scale_alpha_manual(values=c("CTD Cast"=1))+
  scale_shape_manual(values = c("HF Radar Station" = 23))+
  coord_sf(xlim = c(min(dataGamm$Long,na.rm = TRUE)-.05,
                    max(dataGamm$Long,na.rm = TRUE)+.25),
           ylim = c(min(dataGamm$Lat,na.rm = TRUE)-.15,
                    max(dataGamm$Lat,na.rm = TRUE)+.15), expand = c(0.001,0.001,0.001,0.001)) +
  annotation_scale(location = "bl", width_hint = 0.5, text_face = "bold", 
                   text_col = "black",text_cex = 1.2) + 
  annotation_north_arrow(location = "bl", 
                         which_north = "true", 
                         pad_x = unit(0.25, "in"), 
                         pad_y = unit(0.25, "in"), 
                         style = north_arrow_fancy_orienteering) + 
  
  xlab("Longitude") + 
  ylab("Latitude") + 
  # ggtitle("ACCESS Transect Lines (2012-2018)") +
  guides(color = guide_legend(order=1, override.aes = list(size=4, stroke=2)))+
  theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.5), 
        panel.background = element_rect(fill = "aliceblue"),
        axis.title = element_blank(),#element_text(face="bold", size=26),
        axis.text = element_text( face="bold", size=22),
        plot.title = element_text(face="bold", size=25,hjust = 0.5),
        plot.subtitle = element_text( face="bold", size=18,hjust = 0.5),
        plot.caption = element_text( face="bold", size=14,hjust = 0),
        # legend.position= c(0.8,0.9),#"none",
        legend.position= c(0.245,0.115),#"none",
        
        legend.key = element_blank(),
        legend.box.background = element_rect(color = "white",fill="white"),
        legend.background = element_rect(color = "white",fill="white"),
        legend.text = element_text(face="bold", size=20),
        legend.title = element_blank()) #element_text(face="bold", size=24))
p_AllLines
ggsave("./Output/TransectLines.jpg", width=12, height=12, units = "in", dpi=600)

#### Inset Plot #### 
# Center Point of Plot
latOrtho <- 37.125
lonOrtho <- -121.75

# world data
world <- rnaturalearth::ne_countries(scale = 10, returnclass = 'sf')
# Fix polygons so they don't get cut in ortho projection
world  <- st_cast(world, 'MULTILINESTRING') %>%
  st_cast('LINESTRING', do_split=TRUE) %>%
  mutate(npts = npts(geometry, by_feature = TRUE)) %>%
  st_cast('POLYGON')

# Get California polygon data
states <- ne_states(returnclass = "sf") #%>%
states  <- st_cast(states, 'MULTILINESTRING') %>%
  st_cast('LINESTRING', do_split=TRUE) %>%
  mutate(npts = npts(geometry, by_feature = TRUE)) %>%
  st_cast('POLYGON')

# Study Area Box
# Create a data frame with the coordinates of the vertices
df <- data.frame(x = c(-123.8,-122.463,-122.463,-123.8,-123.8), y = c(37.25,37.25,38.5,38.5,37.25)) 
# create matrix from df
vertices <- cbind(df$x, df$y)
# create sf polygon object
polygon <- st_polygon(list(vertices))
# Create an sf object
sfc <- st_sfc(polygon, crs="+proj=longlat +datum=WGS84")

# Define the CRS
crsLONGLAT <- "+proj=longlat +datum=WGS84 +no_defs"
crsLAEA <- paste0('+proj=laea +lat_0=', latOrtho, ' +lon_0=', lonOrtho,' +x_0=4321000 +y_0=3210000 +datum=WGS84 +units=m +no_defs')

# the bouding box polygon in long/lat projection, i.e. axis-aligned
bb <- st_sfc(
  st_polygon(list(cbind(
    c(-127, -115, -115, -127, -127), # x-coordinates (longitudes) of points A,B,C,D (LowerLeft, LowerRight, UpperRight, UpperLeft)
    c(31.3, 31.3, 42.7, 42.7, 31.3)     # y-coordinates (latitudes) of points A,B,C,D
  ))),
  crs = crsLONGLAT)

# now in in LAEA projection
laeabb <- st_transform(bb, crs = crsLAEA)

# the extent of the bounding box in the new projection
b <- st_bbox(laeabb)

# Plot
ggplot() +
  geom_sf(data=world, fill="gray90", color="gray80") + #, aes(fill=continent)
  geom_sf(data=states, color="black") + #, aes(fill=continent)
  geom_sf(data=sfc, color='black',fill="red",alpha=0.4) +
  coord_sf(crs = crsLAEA, xlim = c(b["xmin"], b["xmax"]), ylim = c(b["ymin"], b["ymax"]))+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank())

ggsave(sprintf("./Output/Figure1InsetCA.jpg"), width=14, height=14, units = "in", dpi=600, bg='white') # No gaps
