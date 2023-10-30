#### Surface Current Krill Overview ####
# James Fahlbusch
# Goldbogen Lab, Stanford University
# Version 20230912

# This script 
# 1) Imports and processes Acoustics Data (Krill), Cetacean Sightings and In-Situ Oceanographic Measurements (CTD data) 
#    from Access Cruises between 2012 and 2018
# 2) Downloads the appropriate HF Radar data (netCDF) for external processing
#    NOTE: using the tools at http://transport.me.berkeley.edu/ we restore missing HF
#     Radar data and use the restored data to calculate FTLE
# 3) Import and Extract FTLE data for each cruise, CTD cast and Cetacean Sighting
# 4) Exports processed data for analysis

#### Load Libraries and Functions ####
pkgTest <- function(x)
{
  if (!require(x,character.only = TRUE))
  {
    install.packages(x,dep=TRUE)
    if(!require(x,character.only = TRUE)) stop("Package not found")
  }
}
pkgTest("tidyverse")
pkgTest("ncdf4")
pkgTest("chron")
pkgTest("lattice")
pkgTest("raster")
pkgTest("adehabitatHR")
pkgTest("move")
pkgTest("stats")
pkgTest("lubridate")
pkgTest("rstudioapi")
pkgTest("marmap") # for bathy
pkgTest("rnaturalearth") # for map base layers
pkgTest("metR") # for plotting bathy contours
pkgTest("ggspatial") # for plotting map annotations
pkgTest("fuzzyjoin") # for joining based on a range of datetimes
pkgTest("terra")

# Helper functions specific to this analysis
# NOTE: This includes Focal Functions for Rasters
source('./functions/rasterStackFunctions.r')
# Finds pooled standard deviation
source("./functions/pooled_SD.R")
# Finds Point to Point metrics
source("./functions/pt2pt_fxns.R")
# Finds Sunrise and sunset times
source("./functions/find_Astronomical.R")
# Download HF Radar data
source('./functions/hfRadar_Download.R')
# Process netCDF file dates and times
source("./functions/ncdate.R")
# Note: getNOAA.bathy was not working, but an alternate version on github is used here
source('./functions/getNOAAbathy.R')

`%notin%` <- Negate(`%in%`)

#### Global Variables ####
# Timezone offset to convert to local time
tzOffset <- 'Etc/GMT+7' # all data for this study is in GMT-7
# HF Radar Grid Size (2km or 6km)
gridSize <- 6
# HF Radar extraction - number of seconds to center around timeP (+-) 
timeDifS <- 30*60 
numDays <- 14 # number of days +- of data to download
world <- ne_countries(scale = "large", returnclass = "sf")
# Import Bathymetry for plotting
b = getNOAA.bathy(lon1 = -126.75, lon2 = -118.75, lat1 = 41.5, lat2 = 32.25, resolution = 1)
# convert bathymetry to data frame
bf = fortify.bathy(b)

#### Import and Process Krill Data ####
  ## Acoustics Metadata ##
  # These data represent the processed acoustics data from ACCESS cruises, years 2012-2018.
  # Data are summarized into 200m long x 5m deep cells.
  # A stack of cells makes up a column.
  # The first 5m are not analyzed, so the top of the most shallow cell in any column is 5m.
  # 
  ## Field names and descriptions:	
  # Type: 3 character cruise type identifier (ACC, ACM or ACD)  
  # Cruise: The ACCESS cruise identifier; a 7-character name, using either a GOFyymm (years 2004-2009) or ACCyymm (years 2010 and after) format.
  # LocDate: Local date (in yyyymmdd format)
  # LocTime: Local time (in decimal time)
  # Latitude: Latitude (in decimal degrees) of the center of the 200m column
  # Longitude: Longitude (in decimal degrees) of the center of the 200m column
  # ColumnID: The number of the 200m column (starts at 1 for each cruise)
  # CellID: The number of the 200m x 5m cell in each column (starts at 1 for the shallowest cell). Note that the first 5m are not analyzed.
  # Depth: The depth (in m) of the bottom of the 200m x 5m cell. Note that the first 5m are not analyzed; therefore, the top of the most shallow cell in any column is 5m.
  # g/m^2: The krill biomass in the 200m x 5m cell (in g/m^2)
  # b_depth: The bottom depth (in m) of the column (i.e., the bottom of the deepest 200m x 5m cell in the column)
  # Line: Transect line number/name; if not on a transect, this field is blank.
  # BinLine: The name/number of the ACCESS transect line which is assigned using GIS on where the data fall spatially. 
  #   This will not always match Line but will indicate where the data are closest to (e.g., 1x is near or continuous with line 1 but not really on line 1, 2xN2 is the transit area between line 2 and line N2). 

  ### Import Krill Data ##
  acoustics_DF <- read_csv("./dataRaw/ACCESSdata/Acoustics/acoustics_data_03012022.csv", #file.choose() 
                           col_names = FALSE,
                           col_types = list(col_character(), col_character(),col_double(), col_double(), col_double(),
                                   col_double(), col_double(), col_double(), col_double(), col_double(),col_character(),col_character()))
  colnames(acoustics_DF) <- c("Cruise", "LocDate","LocTime","Lat","Long","ColumnID","CellID", "Depth","KrillBiomass","BotDepth","Line","BinLine")
  acoustics_DF$Type <- substr(acoustics_DF$Cruise,1,3)
  # remove .00 in date from import process
  acoustics_DF$LocDate <-str_remove(acoustics_DF$LocDate,"\\.00000")
  # Parse Date and Time
  acoustics_DF$Date <- ymd(acoustics_DF$LocDate, tz = 'Etc/GMT+7')
  
  dttz<-as.POSIXct(acoustics_DF$Date, tz = 'Etc/GMT+7')
  acoustics_DF$dttz <- dttz + (3600*acoustics_DF$LocTime)
  rm(dttz)
  # # make a GMT Datetime column
  acoustics_DF$dt <- acoustics_DF$dttz
  attr(acoustics_DF$dt, "tzone") <- 'GMT' 
  # remove unused columns
  acoustics_DF <- subset(acoustics_DF, select = -c(LocDate,LocTime))
  acoustics_DF$ColumnID <- as.integer(acoustics_DF$ColumnID)
  acoustics_DF$CellID <- as.integer(acoustics_DF$CellID)
  # Reorder the columns
  acoustics_DF <- acoustics_DF %>% dplyr::select(Type,Cruise,Line,Date,dttz,dt,Lat,Long, everything())
  # Create a column with the Average Density of Krill in each cell (g/m^3)
  acoustics_DF$KrillDensity <- acoustics_DF$KrillBiomass/5 # krill area (g/m^2) / cell height (meters)
  
  ## Save Processed Dataframe 
  save(acoustics_DF,file=paste0("./dataProcessed/ProcessedAcoustics.RData"))
  
  # Create a summary of line times for each cruise for use through analysis
  lineTimes <- acoustics_DF %>% 
    group_by(Type,Cruise,Line) %>% 
    summarize(start=min(dttz),
              end=max(dttz),
              duration = difftime(end, start,units = 'hours'))    
  ## Save Linetimes ##
  save(lineTimes,file=paste0("./dataProcessed/ProcessedAcousticsLineTimes.RData"))
  
  ### Create a 600m Krill Dataset ###
  # To match the spatial resolution of the FTLE data, we combine acoustic grid cells into 600m columns.
  # Partial columns at the end of a line are excluded. 
  # Note: Partial cells along the bottom are retained.
  chunk_length <- 3 # number of 200m X 5m cells to combine 
  
  # Create Column Grouping for all Columns in the acoustics dataset
  acoustics600columns <- acoustics_DF %>% 
    group_by(Cruise,Line,ColumnID, BotDepth) %>% 
    summarise(ColumnID = first(ColumnID),
              BotDepth = first(BotDepth),
              NumCells = n()) %>% 
    group_by(Cruise, Line) %>% 
          # For each Line, Normalize the ColumnID to start at 1
    mutate(ColNorm = ColumnID - min(ColumnID)+1,
           # Group into 3-column chunks for analysis
           Col600 =  ceiling(seq_along(ColNorm) / chunk_length)) %>%
    group_by(Cruise,Line,Col600) %>% 
           # Check if a column is complete (i.e. contains 3 columns)
    mutate(CompleteCols = if_else(n()<chunk_length,0,1)) %>% 
    ungroup() %>% 
    arrange(Cruise, Line, ColumnID)

  # Create a DF to identify the midpoint of each 600 meter column.
  midPoints<-acoustics_DF %>% 
    left_join(dplyr::select(acoustics600columns,Cruise:BotDepth,ColNorm:CompleteCols),by = c("Cruise", "Line", "ColumnID", "BotDepth")) %>% 
    group_by(Type,Cruise,Line,Col600,Lat,Long, dt, dttz, BotDepth) %>% 
    summarise(Lat =first(Lat),Long = first(Long),dt = first(dt),dttz = first(dttz),BotDepth = first(BotDepth)) %>%  
    group_by(Type, Cruise,Line,Col600) %>% 
    # select the middle value of Lat/Long and Time
    summarise(Lat = nth(Lat,2),
              Long = nth(Long,2),
              dt = nth(dt,2),
              dttz = nth(dttz,2),
              # This is the bottom depth of the krill data (max of ~250 meters)
              meanBotDepth = round(mean(BotDepth),0),
              BotDepth = nth(BotDepth,2)) %>%  
    dplyr::filter(!is.na(Lat))
  
  ## Create a dataset of 600 meter columns
  # Note: Incomplete 600 m columns at the end of a transect are discarded
  acoustics600cols <- acoustics_DF %>%
    left_join(dplyr::select(acoustics600columns, -NumCells, -BotDepth),by = c("Cruise", "Line", "ColumnID")) %>% 
    group_by(Type, Cruise, Line, Col600) %>%
    summarize(KrillBiomass = sum(KrillDensity*5),
              # Calculate Density Metrics for each column
              KrillDensity_Mean = mean(KrillDensity,na.rm=TRUE),
              KrillDensity_Median = median(KrillDensity,na.rm=TRUE),
              KrillDensity_Max = max(KrillDensity,na.rm=TRUE),
              KrillDensity_MeanNonZero = mean(KrillDensity[KrillDensity>0]),
              KrillDensity_geoMeanNonZero = 10^(mean(log10( KrillDensity[KrillDensity>0]  ) )),
              # Mean Depth of Krill Layer 
              KrillLayer_MeanDepth = mean(Depth[KrillDensity>0]-2.5), # -2.5 is the center of the 5m cell 
              # Weighted Mean Depth of Krill Layer (weighted by % contribution of krill biomass)
              KrillLayer_MeanDepth_Weighted = weighted.mean(x=Depth[KrillDensity>0]-2.5,
                                                            w= KrillDensity[KrillDensity>0]/sum(KrillDensity)), # -2.5 is the center of the 5m cell 
              NumWith = sum(KrillBiomass>0),
              NumCells = n()) %>% 
    left_join(midPoints, by = c("Type", "Cruise", "Line", "Col600")) %>% 
    dplyr::filter(!is.na(Lat))
  # Replace NANs with NAs in KrillDensiy_MeanNonZero 
  is.na(acoustics600cols) <- do.call(cbind,lapply(acoustics600cols, is.nan)) 
  
  # Determine Day/Night Designation based on Time and Location 
  #  Used for filtering during analysis (i.e. to exclude nightime observations)
  astro_dataD <- find_Astronomical(acoustics600cols$Long, acoustics600cols$Lat, acoustics600cols$dttz)
  acoustics600cols <- left_join(acoustics600cols,dplyr::select(astro_dataD,dttz,astronomical), by = "dttz")
  rm(astro_dataD)
  
  ## Save Processed Dataframe ##
  save(acoustics600cols,file=paste0("./dataProcessed/acoustics600cols.RData"))

#### Import and Process CTD Data ####
  ## CTD Metadata ##
  # In most cases, data were collected with a SBE 19plus (S/N 4041), which is equipped with an SBE 5T submersible pump, 
  #   WETStar fluorometer, Campbell backscatterance sensor, and SBE oxygen sensor (added in early 2010).
  #   Other CTDs used were a SBE 19plus (S/N 5032, which is equipped with a WETStar fluorometer) and a SBE 9+ 
  #   (on Shimada cruises, which is equipped with a SBE 43 oxygen sensor and a WETlabs FLNTU fluorometer).
  # Raw data are processed with the SBEDataProcessing-Win32 software, with the following processing steps:
  #   1. Data Conversion: converting files to .cnv format.
  #   2. Align CTD: defining time lags for each sensor on the CTD.
  #   3. Filter: running the pressure through a low pass filter to smooth the pressure/depth record.
  #   4. Loop Edit: marking "bad" scans with a flag value, which happens with pressure slowdowns or reversals (e.g., ship heave).
  #   5. Derive: using pressure, temperature and conductivity, it computes depth, density, salinity, dissolved oxygen, and percent saturation of oxygen.
  #   6. Bin Average: averages and bins all variables into 1-m depth intervals.
  #   7. ASCII Out: converts the .cnv files into .asc files.
  # The resulting .asc files are then run through two Excel macros to format it into the file imported below
  # 
  ## Field names and descriptions:	
  # ID_Date	A unique identifier of the cast. Includes date, line, station, and cast number.
  # Cruise	The name of the cruise. (Either GOF or ACC, followed by the 2-digit year and month of the cruise.)
  # Date	Date of CTD cast.
  # Time	Time of CTD cast (in local time).
  # Line	Line where CTD cast was completed.
  # Station	Station where CTD cast was completed.
  # CTD_Cast	Cast number within the year. (In some cases, the numbering started over with a new cruise, so letters had to be added to one duplicate sequence of numbers.)
  # Latitude	Latitude of the CTD cast (in decimal degrees).
  # Longitude	Longitude of the CTD cast (in decimal degrees).
  # Distance	Distance to the most inshore station of the line (in kilometers).
  # Bottom	The deepest depth of the cast (in meters).
  # Depth	Depth of the cast (in meters).
  # Pressure	Pressure reading (in psi).
  # Temperature	Temperature reading (in deg C).
  # Conductivity	Conductivity reading (in S/m).
  # Fluorescence	Fluorescence reading (in mg/m3).
  # Backscatterance	Backscatter reading (in NTU).
  # Oxygen raw	Raw dissolved oxygen sensor (in V).
  # Dissolved oxygen	Dissolved oxygen (in mg/L).
  # % saturation	Percent oxygen saturation.
  # Sigma-t	Density reading (in kg/m3). (An estimate of potential density of seawater, minus 1000 kg/m^3)
  # Salinity	Salinity reading (in PSU).
  
  # Sigma-t is a measure the density of seawater at a given temperature. σT is defined as ρ(S,T)-1000 kg m−3, 
  #  where ρ(S,T) is the density of a sample of seawater at temperature T and salinity S, measured in kg m−3, 
  #  at standard atmospheric pressure. For example, a water sample with a density of 1.027 g/cm3 has a σT value of 27.
  
  ## Import CTD Data ##
  ctd_DF <- read_csv("./dataRaw/ACCESSdata/CTD/CTD_data_032022.csv") 
  # Parse Date and Time
  ctd_DF$dttz <- mdy_hms(paste0(ctd_DF$Date," ",ctd_DF$Time), tz = 'Etc/GMT+7')
  ctd_DF$Date <- mdy(ctd_DF$Date, tz = 'Etc/GMT+7')
  # # make a GMT Datetime column
  ctd_DF$dt <- ctd_DF$dttz
  attr(ctd_DF$dt, "tzone") <- 'GMT' 
  # Filter to Cruises overlapping krill data
  ctd_DF <- ctd_DF %>% 
    dplyr::filter(Cruise %in% unique(acoustics_DF$Cruise)) %>% 
    rename("SigmaT"="Sigma-t",
           "Long"="Longitude",
           "Lat"="Latitude" ) %>% 
    mutate(Depth=-Depth, # change convention of depth to + for negative
           PressureDbar = Pressure/1.45037744, # Convert PSI to dbar
           # Calculate the density of the seawater parcel
           Rho = oce::swRho(Salinity, Temperature, PressureDbar, eos="unesco"),
           # Calculate the potential temperature of seawater, denoted theta in the UNESCO system
           # (i.e. What temperature would the parcel have if raised adiabatically to the surface? )
           Theta  = oce::swTheta(Salinity, Temperature, PressureDbar, eos="unesco"),
           # Calculate the potential density of seawater, denoted Sigma-theta in the UNESCO system
           # (i.e. What density would the parcel have if raised adiabatically to the surface? )
           SigmaTheta  = oce::swSigmaTheta(Salinity, Temperature, PressureDbar, eos="unesco"),
           # What density would the parcel have if raised adiabatically to the surface? 
           RhoTheta = oce::swRho(Salinity, oce::swTheta(Salinity, Temperature, PressureDbar, eos="unesco"), 0, eos="unesco")
           )

  # Extract the stations for plotting
  ctd_Stations <- ctd_DF %>% 
    group_by(Cruise,Line,ID_Date, Station, dttz,dt ) %>% 
    summarize(Lat = first(Lat),Long=first(Long)) %>% ungroup()

  # Create a summary of CTD Data
  CTD_Sum <- ctd_Stations %>% 
    group_by(Cruise,Line) %>% 
    summarize(Num_Casts = n())
  
  ## Save Processed Dataframe ##
  save(ctd_DF,file=paste0("./dataProcessed/Processed_ctd_DF.RData"))

 ## Calculate metrics on CTD ##
  ## Variables:
  # Sea Surface Temperature
  #   SST  is the Temperature at the minumum Depth value (range 2-4m)
  # Sea Surface Density
  #   SSD is the Density (i.e. SigmaT) at the minumum Depth value (range 2-4m)
  # 10m Temperature, Salinity, Density
  # 20m Temperature, Salinity, Density
  # 30m Temperature, Salinity, Density
  # 40m Temperature, Salinity, Density
  # 50m Temperature, Salinity, Density
  # Depth of SigmaT 26 

  # Variables will be added incrementally to this dataframe
  datasetCTD_Vars <- ctd_DF 
  castSummary <- ctd_Stations 
  ## Calculate SST
  sst <- ctd_DF %>% 
    arrange(ID_Date,Depth) %>% 
    group_by(ID_Date) %>% 
    summarize(SST = first(Temperature),
              SSTdepth = first(Depth),
              SSD = first(SigmaT),
              SigmaTheta0m = first(SigmaTheta),
              minDepth= min(Depth),
              maxDepth=max(Depth))
  # remake CTD dataframe to include sst
  castSummary <- left_join(castSummary,sst,by = c("ID_Date"))
  datasetCTD_Vars <- left_join(datasetCTD_Vars,sst,by = c("ID_Date"))
  
  ## Temperature, Salinity and Density at 10m
  vals10m <- ctd_DF %>% 
    arrange(ID_Date,Depth) %>% 
    group_by(ID_Date) %>%
    dplyr::filter(Depth==10) %>% 
    summarize(T10m = first(Temperature),
              S10m = first(Salinity),
              D10m = first(SigmaT),
              SigmaTheta10m = first(SigmaTheta)) 
  # remake CTD dataframe to include vals10m
  castSummary <- left_join(castSummary,vals10m,by = c("ID_Date"))
  datasetCTD_Vars <- left_join(datasetCTD_Vars,vals10m,by =  c("ID_Date"))
  
  ## Temperature, Salinity and Density at 20m
  vals20m <- ctd_DF %>% 
    arrange(ID_Date,Depth) %>% 
    group_by(ID_Date) %>%
    dplyr::filter(Depth==20) %>% 
    summarize(T20m = first(Temperature),
              S20m = first(Salinity),
              D20m = first(SigmaT),
              SigmaTheta20m = first(SigmaTheta)) 
  # remake CTD dataframe to include vals10m
  castSummary <- left_join(castSummary,vals20m,by = c("ID_Date"))
  datasetCTD_Vars <- left_join(datasetCTD_Vars,vals20m,by =  c("ID_Date"))
  
  ## Temperature, Salinity and Density at 30m
  vals30m <- ctd_DF %>% 
    arrange(ID_Date,Depth) %>% 
    group_by(ID_Date) %>%
    dplyr::filter(Depth==30) %>% 
    summarize(T30m = first(Temperature),
              S30m = first(Salinity),
              D30m = first(SigmaT),
              SigmaTheta30m = first(SigmaTheta)) 
  # remake CTD dataframe to include vals30m
  castSummary <- left_join(castSummary,vals30m,by = c("ID_Date"))
  datasetCTD_Vars <- left_join(datasetCTD_Vars,vals30m,by =  c("ID_Date"))
  
  ## Temperature, Salinity and Density at 40m
  vals40m <- ctd_DF %>% 
    arrange(ID_Date,Depth) %>% 
    group_by(ID_Date) %>%
    dplyr::filter(Depth==40) %>% 
    summarize(T40m = first(Temperature),
              S40m = first(Salinity),
              D40m = first(SigmaT),
              SigmaTheta40m = first(SigmaTheta)) 
  # remake CTD dataframe to include vals40m
  castSummary <- left_join(castSummary,vals40m,by = c("ID_Date"))
  datasetCTD_Vars <- left_join(datasetCTD_Vars,vals40m,by =  c("ID_Date"))
  
  ## Temperature, Salinity and Density at 50m
  vals50m <- ctd_DF %>% 
    arrange(ID_Date,Depth) %>% 
    group_by(ID_Date) %>%
    dplyr::filter(Depth==50) %>% 
    summarize(T50m = first(Temperature),
              S50m = first(Salinity),
              D50m = first(SigmaT),
              SigmaTheta50m = first(SigmaTheta)) 
  # remake CTD dataframe to include vals50m
  castSummary <- left_join(castSummary,vals50m,by = c("ID_Date"))
  datasetCTD_Vars <- left_join(datasetCTD_Vars,vals50m,by =  c("ID_Date"))
  
  ## Depth at Density (SigmaT) of 25
  valsSigma25 <- ctd_DF %>%
    arrange(ID_Date,Depth) %>%
    group_by(ID_Date) %>%
    # Must be deeper than SST
    dplyr::filter(SigmaTheta>=25) %>%
    summarize(Sigma25 = first(Depth))
  # remake CTD dataframe to include valsSigma26
  castSummary <- left_join(castSummary,valsSigma25,by = c("ID_Date"))
  datasetCTD_Vars <- left_join(datasetCTD_Vars,valsSigma25,by =  c("ID_Date"))
  
  ## Depth at Density (SigmaT) of 25.5
  valsSigma25_5 <- ctd_DF %>%
    arrange(ID_Date,Depth) %>%
    group_by(ID_Date) %>%
    # Must be deeper than SST
    dplyr::filter(SigmaTheta>=25.5) %>%
    summarize(Sigma25_5 = first(Depth))
  # remake CTD dataframe to include valsSigma26
  castSummary <- left_join(castSummary,valsSigma25_5,by = c("ID_Date"))
  datasetCTD_Vars <- left_join(datasetCTD_Vars,valsSigma25_5,by =  c("ID_Date"))
  
  ## Depth at Density (SigmaT) of 25.6
  valsSigma25_6 <- ctd_DF %>%
    arrange(ID_Date,Depth) %>%
    group_by(ID_Date) %>%
    # Must be deeper than SST
    dplyr::filter(SigmaTheta>=25.6) %>%
    summarize(Sigma25_6 = first(Depth))
  # remake CTD dataframe to include valsSigma26
  castSummary <- left_join(castSummary,valsSigma25_6,by = c("ID_Date"))
  datasetCTD_Vars <- left_join(datasetCTD_Vars,valsSigma25_6,by =  c("ID_Date"))
  
  ## Depth at Density (SigmaT) of 26
  valsSigma26 <- ctd_DF %>%
    arrange(ID_Date,Depth) %>%
    group_by(ID_Date) %>%
    # Must be deeper than SST
    dplyr::filter(SigmaTheta>=26) %>%
    summarize(Sigma26 = first(Depth))
  # remake CTD dataframe to include valsSigma26
  castSummary <- left_join(castSummary,valsSigma26,by = c("ID_Date"))
  datasetCTD_Vars <- left_join(datasetCTD_Vars,valsSigma26,by =  c("ID_Date"))
  # clean up
  rm(sst,vals10m,vals20m,vals30m,vals40m,vals50m,valsSigma25,valsSigma26,valsSigma25_5,valsSigma25_6)
  # Casts 20120618-6-E-6 and 20120618-6-EX-5 have the same timestamp but are 7nm apart. Drop both.
  castSummary <- castSummary %>% dplyr::filter(ID_Date %notin% c("20120618-6-E-6","20120618-6-EX-5"))
  
  ## Save castSummary to reload later
  save(castSummary,file=paste0("./dataProcessed/castSummary.RData"))
  ## Save the cast dataframe datasetCTD_Vars
  save(datasetCTD_Vars,file=paste0("./dataProcessed/datasetCTD_Vars.RData"))

#### Import and Process Cetacean Data ####
  #  Rorqual data from the ACCESS cruises 2012-2018. 
  #    Notes: Unidentified whales are included in this list; If data are missing for Distance from ship, Perpendicular distance from ship, Sighting Latitude and/or Sighting 
  #           Longitude, either bearing or reticle were missing from the original data, and these values could not be 
  #           calculated OR the observation was off effort and the bearing and reticle were not recorded.			
  #   
  ## Field names and descriptions:	
  # Date	Date	YYMMDD	
  # LocalTime	Time	HHMMSS	
  # Lat	Ship's Latitude	DD.DDDDDD	
  # Long	Ship's Longitude	DDD.DDDDDD	
  # ShipHeading	Ship's Heading	Degrees	
  # ShipSpeed	Ship's Speed	Knots	
  # EventCode	Event Code		4 - On effort; 7 - Off effort
  # SpCode	Species Code		BLWH - Blue Whale; HUWH - Humpback Whale; FIWH - Fin Whale; GRWH - Gray Whale; MIWH - Minke Whale
  # DisFrShipNm	Distance from ship	nautical miles	
  # PerpDistNm	Perpendicular distance from ship to animal (how far the animal is from the transect line)	nautical miles	
  # SightLat	Sighting Latitude - Animal's Latitude	DD.DDDDDD	
  # SightLong	Sighting Longitude - Animal's Longitude	DDD.DDDDDD	
  # CommentCode	Comment Code		1 - no; 2 - yes
  # CommentText	Comment Text	
  
  # Import Raw Cetacean data
  cetacean_DF <- read_csv("./dataRaw/ACCESSdata/Cetaceans/cetacean_data.csv",
                          col_types = list(col_character(), col_character(), 
                                           col_double(), col_double(), col_double(), 
                                           col_double(), col_double(), col_character(),
                                           col_double(), col_double(), col_double(),
                                           col_double(), col_double(), col_character()
                          )) %>% 
    rename("SightLong"="SiteLong") 
  cetacean_DF$LocalTime <- if_else(as.integer(cetacean_DF$LocalTime)>=100000,cetacean_DF$LocalTime,
                                   paste0("0",cetacean_DF$LocalTime))
  # create a datetime column
  cetacean_DF$dttz <- ymd_hms(paste0(cetacean_DF$Date," ", cetacean_DF$LocalTime), tz = 'Etc/GMT+7')
  # make a GMT Datetime column
  cetacean_DF$dt <- cetacean_DF$dttz
  attr(cetacean_DF$dt, "tzone") <- 'GMT' 

  ## Link cetacean data to the krill Cruise lines 
  cetacean_DF <- fuzzy_left_join(cetacean_DF, lineTimes, 
                                 by=c("dttz"="start", "dttz"="end"),
                                 match_fun=list(`>=`, `<=`)) %>% 
    dplyr::select(Type,Cruise,Line,dt,dttz,Lat,Long,everything(),-start,-end,-Date,-LocalTime)

  # Filter to dates overlapping krill data
  cetacean_DF <- cetacean_DF %>% 
    dplyr::filter(Cruise %in% unique(acoustics_DF$Cruise)) %>% 
    dplyr::filter(dttz <= max(acoustics_DF$dttz),
                  dttz >= min(acoustics_DF$dttz),
                  EventCode == 4) %>% 
    # Remove Sightings not associated with a cruise and without a location
    dplyr::filter(!is.na(Cruise),!is.na(SightLat)) %>% 
    dplyr::select(-duration) %>% 
    # add a sighting number
    arrange(Cruise, dttz) %>% 
    group_by(Cruise) %>% 
    mutate(SightingNum = seq(1,n())) %>% 
    ungroup()

  # Add 600meter Acoustics Column for analysis
  cetacean_Cols <- data.frame()
  cetacean_DF$Col600 <- NA
  for(i in 1:length(cetacean_DF$dttz)){
    # sighting must be within +-2.5 min of a column to be linked
    tempOut <- acoustics600cols %>% 
      dplyr::filter(dttz>= cetacean_DF$dttz[i]-(2.5*60), dttz < cetacean_DF$dttz[i]+(2.5*60) ) %>% 
      # calculate the time difference (in minutes) and arrange smallest first
      mutate(tDif = abs(as.numeric(difftime(cetacean_DF$dttz[i], .$dttz, units = "secs")/60))) %>% 
      arrange(tDif)
    if(dim(tempOut)[1] != 0 ){
      cetacean_Cols <- rbind(cetacean_Cols,tempOut %>% 
                               slice(1) %>% 
                               dplyr::select(Type,Cruise,Line,Col600,tDif) %>% 
                               mutate(SightingNum = cetacean_DF$SightingNum[i]))
      #Add the 600 meter column to the Cetacean Dataset
      cetacean_DF$Col600[i] <- tempOut$Col600[1]
    }
    rm(tempOut)
  }
  rm(i)
  
  # Save Processed Cetacean Data
  save(cetacean_DF,file=paste0("./dataProcessed/cetacean_DF.RData"))
  save(cetacean_Cols,file=paste0("./dataProcessed/cetacean_Cols.RData"))
  
  # Create a summary of whale sightings
  cetacean_Sum <- cetacean_DF %>% 
    group_by(Cruise,Line) %>% 
    summarize(NumSightings = n())

#### Download Surface Current Data ####
  ## Download HF Radar data
  # NOTE: This script downloads SC data with a numDays buffer on either end of the Cruise dates and
  #       adds a degreeBuffer around the study area bounding box to avoid edge-effects
  # The bounding box is the same for all ACCESS Surveys in the Gulf of the Farallones
  
  useOldNC <- FALSE # if true, script will search for a local NC file before downloading a new one
  numDays <- 14 # number of days of HF Radar data before and after deployment
  degreeBuffer <- 1 # number of degrees around the min/max lat/longs of Cruise Locations
  
  # Bounding Box around each of the survey Types for HF Radar download
  surveyBB <- acoustics_DF %>% 
    group_by(Type) %>% 
    summarize(minLat = min(Lat),
              maxLat = max(Lat),
              minLong = min(Long),
              maxLong = max(Long))
  
  ncPattern <- paste0("_HFRADAR_US_West_Coast_", gridSize, "km_Resolution_Hourly.nc")
  for(i in 1:length(unique(acoustics_DF$Cruise))){
    # Get the CruiseID and time zone for this deployment
    Cruise <- unique(acoustics_DF$Cruise)[i]
    cat(paste0("\nProcessing: ",  unique(acoustics_DF$Cruise)[i],"\n"))
    # filter by Cruise
    data_Sub <- acoustics_DF %>% dplyr::filter(Cruise == unique(acoustics_DF$Cruise)[i],
                                               !is.na(Lat),!is.na(Long)) 
    Type <- unique(data_Sub$Type) # determine the type for bounding box
    # Determine the start and end times of deployment, +- 1 hour
    start <- min(data_Sub$dt)-(3600) 
    end <- max(data_Sub$dt)+(3600)
    
    # Only download new NC file if it doesn't exist locally
    if(useOldNC & file.exists(paste0("./dataRaw/",Cruise,ncPattern))){ 
      ncname <- paste0("./dataRaw/",Cruise,ncPattern)
      cat(paste0("Using previously downloaded HF Radar data:\n",ncname))
    }else{ # Download the appropriate NC file using parameters from tag 
      ncname <- hfRadar_Download(paste0("data/RawHFRadarData/", unique(acoustics_DF$Cruise)[i]),
                                 dStart = lubridate::date(start)-numDays,tStart = paste0(lubridate::hour(start),':00:00'),
                                 dEnd = lubridate::date(end)+numDays, tEnd = paste0(lubridate::hour(end),':00:00'), 
                                 bbox = c(round(surveyBB$maxLat[surveyBB$Type==Type]+degreeBuffer,2),
                                          round(surveyBB$minLat[surveyBB$Type==Type]-degreeBuffer,2),
                                          round(surveyBB$minLong[surveyBB$Type==Type]-degreeBuffer,2),
                                          round(surveyBB$maxLong[surveyBB$Type==Type]+degreeBuffer,2)), 
                                 grid = gridSize)
    }
    rm(ncname, data_Sub, start, end, Cruise, Type)
  }
  
#### Process FTLE for Acoustics ####
  # A study in the same region found good agreement between SST and LCS during summer months (Gough et al 2016). Here 
  # we use the 2 integration durations used in that study (1 and 5 days), the integration used in a blue whale paper 
  # in the same region (2 days, Fahlbusch et al 2022), as well as a longer period for comparison (10 days) to examine the interaction between 
  # krill density and distribution and LCS during summer months. 
  # References:
  #   Matt K Gough, Ad Reniers, M. Josefina Olascoaga, Brian K Haus, Jamie MacMahan, Jeff Paduan, Chris Halle,
  #   Lagrangian Coherent Structures in a coastal upwelling environment,
  #   Continental Shelf Research, 2016, https://doi.org/10.1016/j.csr.2016.09.007
  #   Fahlbusch, J. A., Czapanskiy, M. F., Calambokidis, J., Cade, D. E., Abrahms, B., Hazen, E. L., & Goldbogen, J. A.,
  #   Blue whales increase feeding rates at fine-scale ocean features, 
  #   Proceedings of the Royal Society B, 2022, https://doi.org/10.1098/rspb.2022.1180

  ## Import the externally processed FTLE data (netCDF) and Extract Values at Survey Locations
  ## Restoration Settings:
  #   Data was restored using automatic land detection and restored up to shoreline. 
  #   Concave Hull around available points with 10km alpha
  ## FTLE Settings: 
  #   Include Land and Apply a Free-Slip boundary condition
  #   Initial Time: 14 days before end of data 
  #   Number of Initial Times: One for each hour of interest (# of days of surveys * 24) 
  #   Initial Times Increment: 1 Hour
  #   Integration Duration: 1, 2, 5 and 10 days
  #   Reverse Direction of Time 
  ## FTLE is calculated with a grid of tracers having 10x the resolution of the input data.
  #   The final resolution varies slightly by deployment, but is roughly 590m or 6km HF radar resolution/10.
  #   All FTLE Files for a deployment (i.e. 24, 48, 120, 240 Hrs) have the same temporal and spatial resolution.
  ## FTLE data is in hourly increments (due to the temporal resolution of the HF Radar surface current measurements), 
  #   and some survey transect lines span multiple FTLE measurement increments (mean = 3.02). To avoid introducing
  #   artificial jumps in measurement values at hourly changes, for this analysis we calculated a spatial mean of 
  #   the FTLE layers for each line. 
  #     We use the lineTimes dataframe to determine which raster stack layers to average for each line

  FTLEdatapath <- "./dataRaw/FTLE/"
  FTLE_Integrations <- c("24","48","120", "240")

  locations_DF <- acoustics_DF %>% 
    dplyr::select(Cruise, dttz,dt, Lat,Long) %>%
    unique()
  
  # Create empty dataframes to hold all the extracted data
  lineBased_DF_FTLE <- data.frame()
  for(i in 1:length(unique(lineTimes$Cruise))){
    cruise <- unique(lineTimes$Cruise)[i]
    cat(paste0("\nProcessing: ",  cruise))
    cat(paste0("\nCruise ",  i, " of ", length(unique(lineTimes$Cruise))))
    ptm <- proc.time()
    # filter by Cruise 
    data_Sub <- locations_DF %>% dplyr::filter(Cruise == unique(lineTimes$Cruise)[i],
                                        !is.na(Lat),!is.na(Long))
    # We use the lineTimes dataframe to determine which raster stack layers to average 
    lines_Sub <- lineTimes %>% dplyr::filter(Cruise == unique(lineTimes$Cruise)[i])
    # Create an empty dataframe for extracted FTLE values
    ftle_ExtractAll <- data.frame()
    for(ii in 1:length(FTLE_Integrations)){
      ptm1<-proc.time()
      gc()
      cat(paste0("\nFTLE: ",  FTLE_Integrations[ii],"hrs"))
      ncname <- paste0(FTLEdatapath,cruise,"_FTLE_",FTLE_Integrations[ii],"Hrs.nc")
      ncin_FTLE <- nc_open(ncname, verbose=FALSE)
      # FTLE Time
      time_FTLE <- ncdate(ncin_FTLE)
      # convert to local time 
      time_FTLE_Local <- time_FTLE
      attr(time_FTLE_Local, "tzone") <- tzOffset # change the timezone to tzOffset
      # Create a FTLE Raster Stack and save for later use
      FTLE_stack <- stack(ncname, varname = "FTLE")

      ### Raster Extraction ###
      ## Create a spatial points dataframe of cruise locations
      xy <- coordinates(cbind(data_Sub$Long, data_Sub$Lat))
      dataSP <- SpatialPointsDataFrame(xy, data_Sub, proj4string = CRS("+proj=longlat +ellps=WGS84 +no_defs"), match.ID = FALSE)
      cat("\nExtracting FTLE Data\n")
      # Create an empty dataframe for extracted FTLE values
      ftle_Extract <- data.frame()
      for(j in 1:length(lines_Sub$Line)){
        cat(".") # progress bar...
        # Create a Mean Layer for the line
        lineDF = lines_Sub[j,]
        # Find the Stack layer that correspond to the Start/End time of the Line
        tempDf <- as.data.frame(time_FTLE_Local) 
        tempDf$layer <- seq(1,length(time_FTLE_Local),1)
        tempDf <- tempDf %>% 
          dplyr::filter(time_FTLE_Local > lineDF$start-timeDifS,
                        time_FTLE_Local <= lineDF$end+timeDifS) %>% 
          summarize(minJ = min(layer),
                    maxJ = max(layer))
        # Create a mean FTLE layer for the line
        FTLE_Mean <- mean(FTLE_stack[[tempDf$minJ:tempDf$maxJ]], na.rm = TRUE)
        # save(FTLE_Mean,file= paste0("./dataProcessed/FTLEstacks/",unique(data_Sub$Cruise),"-",lineDF$Line, "_FTLEmean-",
        #                             FTLE_Integrations[ii], "hrs.RData"))
        # Subset by line start/end times
        tempSP <- dataSP[(dataSP$dttz <= lineDF$end & dataSP$dttz >= lineDF$start),]
        if(dim(tempSP)[1] != 0 ){ # make sure there is data for this period
          # extract the FTLE values from the focal layer
          ftle_Values <- raster::extract(FTLE_Mean,tempSP,layer = 1, nl=1, cellnumbers = TRUE)
          colnames(ftle_Values) <- c('rasterCellFTLE_L',paste0('ftle',FTLE_Integrations[ii],"_L"))          
          # Calculate the overall mean of FTLE for the layer
          ftle_meanLayer <- cellStats(FTLE_Mean,'mean',na.rm=TRUE)
          # combine into a single dataframe
          tempSP <- as.data.frame(tempSP)
          if(ii==1){
            tempSP <- tempSP %>% mutate(Line = lineDF$Line,
                                        lineMinI = tempDf$minJ,
                                        lineMaxI = tempDf$maxJ)
          }
          tempSP <- cbind(tempSP,ftle_Values)
          tempSP$ftle_meanLayer <- ftle_meanLayer
          names(tempSP)[names(tempSP)=='ftle_meanLayer'] <- paste0('ftle_meanLayer',FTLE_Integrations[ii],"_L")
          # append it to an output dataframe
          ftle_Extract <- rbind(ftle_Extract, tempSP)
          rm(ftle_Values, ftle_meanLayer) # cleanup
        }
        rm(tempSP, lineDF,tempDf, FTLE_Mean)
      }
      rm(j)
      
      if(ii==1){
        ftle_ExtractAll <- ftle_Extract 
      } else{
        ftle_ExtractAll <- cbind(ftle_ExtractAll,dplyr::select(ftle_Extract, starts_with("ftle")))
      } 
      nc_close(ncin_FTLE)
      rm(ncin_FTLE, ncname,time_FTLE,time_FTLE_Local,ftle_Extract,xy,dataSP)
      tPTM1 <-  proc.time() - ptm1
      cat(paste0("\nElapsed time: ", round(tPTM1[3],1), " seconds"))
      rm(ptm1,tPTM1)
    }
    # Save to an output dataframe 
    if(i==1){
      lineBased_DF_FTLE <- ftle_ExtractAll 
    } else{
      lineBased_DF_FTLE <- rbind(lineBased_DF_FTLE, ftle_ExtractAll)
    } 
    tPTM <-  proc.time() - ptm
    cat(paste0("\nOverall elapsed time: ", round(tPTM[3]/60,1), " minutes"))
    rm(ptm, tPTM)
    gc()
  }
  # clean up unused variables
  rm(i,ii, cruise, ftle_ExtractAll, FTLE_Integrations,FTLEdatapath,data_Sub)
  
  # Remove the remnants of the Spatial Points DF conversion
  lineBased_DF_FTLE <- lineBased_DF_FTLE %>% dplyr::select(-c(coords.x1, coords.x2))   
  # Check for INF
  sapply(lineBased_DF_FTLE, function(x) sum(is.infinite(x))>0)
  # convert INF values to na (only seems to exist in Max/Min)
  is.na(lineBased_DF_FTLE) <- do.call(cbind,lapply(lineBased_DF_FTLE, is.infinite))
  is.na(lineBased_DF_FTLE) <- do.call(cbind,lapply(lineBased_DF_FTLE, is.nan)) 
  ##Save extracted FTLE for later import
  save(lineBased_DF_FTLE,file=paste0("./dataProcessed/lineBased_DF_FTLE.RData"))
  
  # Join Locations to Acoustics DF and Add 600 Meter Columns for future linking
  acoustics_DF_FTLE_W600 <- left_join(acoustics_DF, 
                                      lineBased_DF_FTLE, 
                                      by = c("Cruise","Line","dttz","dt","Long","Lat")) %>%
    #Join 600m Columns to Acoustics DF 
    left_join(dplyr::select(acoustics600columns, Cruise, Line, ColumnID, ColNorm, Col600),
              by = c("Cruise", "Line", "ColumnID")) %>% 
    arrange(Cruise, Line, ColumnID)
  ##Save Full Acoustics dataset (200m x 5m) with FTLE and 600 cols
  save(acoustics_DF_FTLE_W600,file=paste0("./dataProcessed/acoustics_DF_FTLE_Line_w600.RData"))
  
  ## Join FTLE Data to 600 meter columns for analysis
  acousticsFTLE600cols <- acoustics600cols %>%
    left_join(lineBased_DF_FTLE, by = c("Cruise","Line","dttz","dt","Long","Lat")) 
  ##Save 600m Acoustics with FTLE
  save(acousticsFTLE600cols,file=paste0("./dataProcessed/acousticsFTLE600cols.RData"))

  
  #### Process FTLE for Cetacean Sightings ####
  # Here we use the sighting location for FTLE data extraction
  
  FTLEdatapath <- "./dataRaw/FTLE/"
  FTLE_Integrations <- c("24","48","120", "240")
  whaleLocs_DF <- cetacean_DF %>% 
    dplyr::filter(!is.na(Cruise),!is.na(SightLat)) %>% 
    dplyr::select(Cruise,dttz,dt, SightLat,SightLong) %>%
    unique() %>% 
    rename("Lat" = "SightLat", "Long" = "SightLong")
  
  for(i in 1:length(unique(whaleLocs_DF$Cruise))){
    cruise <- unique(whaleLocs_DF$Cruise)[i]
    cat(paste0("\nProcessing: ",  cruise))
    cat(paste0("\nCruise ",  i, " of ", length(unique(lineTimes$Cruise))))
    ptm <- proc.time()
    # filter by Cruise 
    data_Sub <- whaleLocs_DF %>% dplyr::filter(Cruise == unique(whaleLocs_DF$Cruise)[i],
                                               !is.na(Lat),!is.na(Long))
    ftle_ExtractAll <- data.frame()
    for(ii in 1:length(FTLE_Integrations)){
      ptm1<-proc.time()
      gc()
      cat(paste0("\nFTLE: ",  FTLE_Integrations[ii],"hrs"))
      ncname <- paste0(FTLEdatapath,cruise,"_FTLE_",FTLE_Integrations[ii],"Hrs.nc")
      ncin_FTLE <- nc_open(ncname, verbose=FALSE)
      # FTLE Time
      time_FTLE <- ncdate(ncin_FTLE)
      # convert to local time 
      time_FTLE_Local <- time_FTLE
      attr(time_FTLE_Local, "tzone") <- tzOffset # change the timezone to tzOffset
      # Create a FTLE Raster Stack and save for later use
      FTLE_stack <- stack(ncname, varname = "FTLE")
      ### Raster Extraction ###
      ## Create a spatial points dataframe of cruise locations
      xy <- coordinates(cbind(data_Sub$Long, data_Sub$Lat))
      dataSP <- SpatialPointsDataFrame(xy, data_Sub, proj4string = CRS("+proj=longlat +ellps=WGS84 +no_defs"), match.ID = FALSE)
      cat("\nExtracting FTLE Data\n")
      # Create an empty dataframe for extracted FTLE values
      ftle_Extract <- data.frame()
      for(j in 1:length(time_FTLE)){ 
        cat(".") # progress bar...
        # subset by the time of FTLE observation (+- timeDifs)
        tempSP <- dataSP[(dataSP$dt <= time_FTLE[j]+timeDifS & dataSP$dt > time_FTLE[j]-timeDifS),]
        if(dim(tempSP)[1] != 0 ){ # make sure there is data for this period
          # extract the FTLE values from the focal layer
          ftle_Values <- raster::extract(FTLE_stack,tempSP,layer = j, nl=1, cellnumbers = TRUE)
          colnames(ftle_Values) <- c('rasterCellFTLE',paste0('ftle',FTLE_Integrations[ii]))          
          ftle_meanLayer <- cellStats(FTLE_stack[[j]],'mean',na.rm=TRUE)
          # combine into a single dataframe
          tempSP <- as.data.frame(tempSP)
          tempSP$extractNumFTLE <- j          
          tempSP <- cbind(tempSP,ftle_Values)
          tempSP$ftle_meanLayer <- ftle_meanLayer
          names(tempSP)[names(tempSP)=='ftle_meanLayer'] <- paste0('ftle_meanLayer',FTLE_Integrations[ii])
          # append it to an output dataframe
          ftle_Extract <- rbind(ftle_Extract, tempSP)
          rm(ftle_Values, ftle_meanLayer) # cleanup
        }
        rm(tempSP)
      }
      rm(j)
      if(ii==1){
        ftle_ExtractAll <- ftle_Extract 
      } else{
        ftle_ExtractAll <- cbind(ftle_ExtractAll,dplyr::select(ftle_Extract, starts_with("ftle")))
      } 
      nc_close(ncin_FTLE)
      rm(ncin_FTLE, ncname,time_FTLE,time_FTLE_Local,ftle_Extract,xy,dataSP)
      tPTM1 <-  proc.time() - ptm1
      cat(paste0("\nElapsed time: ", round(tPTM1[3],1), " seconds"))
      rm(ptm1,tPTM1)
    }
    # Save to an output dataframe 
    if(i==1){
      whaleLocs_DF_FTLE <- ftle_ExtractAll 
    } else{
      whaleLocs_DF_FTLE <- rbind(whaleLocs_DF_FTLE, ftle_ExtractAll)
    } 
    tPTM <-  proc.time() - ptm
    cat(paste0("\nOverall elapsed time: ", round(tPTM[3]/60,1), " minutes"))
    rm(ptm, tPTM)
    gc()
  }
  # clean up unused variables
  rm(i,ii, cruise, ftle_ExtractAll, FTLE_Integrations,FTLEdatapath,data_Sub)
  # Remove the remnants of the Spatial Points DF conversion
  whaleLocs_DF_FTLE <- whaleLocs_DF_FTLE %>% dplyr::select(-c(coords.x1, coords.x2))   
  # Check for INF
  sapply(whaleLocs_DF_FTLE, function(x) sum(is.infinite(x))>0)
  # convert INF values to na (only seems to exist in Max/Min)
  is.na(whaleLocs_DF_FTLE) <- do.call(cbind,lapply(whaleLocs_DF_FTLE, is.infinite))
  is.na(whaleLocs_DF_FTLE) <- do.call(cbind,lapply(whaleLocs_DF_FTLE, is.nan)) 
  ##Save extracted FTLE for later import
  save(whaleLocs_DF_FTLE,file=paste0("./dataProcessed/whaleLocs_DF_FTLE.RData"))
  # Join to Cetacean Dataset  
  cetacean_DF_FTLE <- left_join(cetacean_DF,whaleLocs_DF_FTLE, 
                                by = c("Cruise" = "Cruise","dttz"="dttz","dt"="dt",
                                       "SightLong"="Long","SightLat"="Lat"))
  save(cetacean_DF_FTLE,file=paste0("./dataProcessed/cetacean_DF_FTLE.RData"))
  
  
  #### Process FTLE for CTD Cast locations ####  
  FTLEdatapath <- "./dataRaw/FTLE/"
  FTLE_Integrations <- c("24","48","120", "240")
  for(i in 1:length(unique(castSummary$Cruise))){
    cruise <- unique(castSummary$Cruise)[i]
    cat(paste0("\nProcessing: ",  cruise))
    cat(paste0("\nCruise ",  i, " of ", length(unique(castSummary$Cruise))))
    ptm <- proc.time()
    # filter by Cruise 
    data_Sub <- castSummary %>% 
      dplyr::filter(Cruise == unique(castSummary$Cruise)[i],
                                               !is.na(Lat),!is.na(Long)) %>% 
      arrange(dt)
    ftle_ExtractAll <- data.frame()
    for(ii in 1:length(FTLE_Integrations)){
      ptm1<-proc.time()
      gc()
      cat(paste0("\nFTLE: ",  FTLE_Integrations[ii],"hrs"))
      ncname <- paste0(FTLEdatapath,cruise,"_FTLE_",FTLE_Integrations[ii],"Hrs.nc")
      ncin_FTLE <- nc_open(ncname, verbose=FALSE)
      # FTLE Time
      time_FTLE <- ncdate(ncin_FTLE)
      # convert to local time 
      time_FTLE_Local <- time_FTLE
      attr(time_FTLE_Local, "tzone") <- tzOffset # change the timezone to tzOffset
      # Create a FTLE Raster Stack and save for later use
      FTLE_stack <- stack(ncname, varname = "FTLE")
      ### Raster Extraction ###
      ## Create a spatial points dataframe of cruise locations
      xy <- coordinates(cbind(data_Sub$Long, data_Sub$Lat))
      dataSP <- SpatialPointsDataFrame(xy, data_Sub, proj4string = CRS("+proj=longlat +ellps=WGS84 +no_defs"), match.ID = FALSE)
      cat("\nExtracting FTLE Data\n")
      # Create an empty dataframe for extracted FTLE values
      ftle_Extract <- data.frame()
      for(j in 1:length(time_FTLE)){ 
        #cat(".") # progress bar...
        # subset by the time of FTLE observation (+- timeDifs)
        tempSP <- dataSP[(dataSP$dt <= time_FTLE[j]+timeDifS & dataSP$dt > time_FTLE[j]-timeDifS),]
        if(dim(tempSP)[1] != 0 ){ # make sure there is data for this period
          # extract the FTLE values from the focal layer
          ftle_Values <- raster::extract(FTLE_stack,tempSP,layer = j, nl=1, cellnumbers = TRUE)
          colnames(ftle_Values) <- c('rasterCellFTLE',paste0('ftle',FTLE_Integrations[ii]))          
          ftle_meanLayer <- cellStats(FTLE_stack[[j]],'mean',na.rm=TRUE)
          # combine into a single dataframe
          tempSP <- as.data.frame(tempSP)
          tempSP$extractNumFTLE <- j          
          tempSP <- cbind(tempSP,ftle_Values)
          tempSP$ftle_meanLayer <- ftle_meanLayer
          names(tempSP)[names(tempSP)=='ftle_meanLayer'] <- paste0('ftle_meanLayer',FTLE_Integrations[ii])
          # append it to an output dataframe
          ftle_Extract <- rbind(ftle_Extract, tempSP)
          rm(ftle_Values, ftle_meanLayer) # cleanup
        }
        rm(tempSP)
      }
      rm(j)
      
      if(ii==1){
        ftle_ExtractAll <- ftle_Extract 
      } else{
        ftle_ExtractAll <- left_join(ftle_ExtractAll,dplyr::select(ftle_Extract,dttz, starts_with("ftle")),by="dttz")
      } 
      nc_close(ncin_FTLE)
      rm(ncin_FTLE, ncname,time_FTLE,time_FTLE_Local,ftle_Extract,xy,dataSP)
      tPTM1 <-  proc.time() - ptm1
      cat(paste0("\nElapsed time: ", round(tPTM1[3],1), " seconds"))
      rm(ptm1,tPTM1)
    }
    # Save to an output dataframe 
    if(i==1){
      castSummary_FTLE <- ftle_ExtractAll 
    } else{
      castSummary_FTLE <- rbind(castSummary_FTLE, ftle_ExtractAll)
    } 
    tPTM <-  proc.time() - ptm
    cat(paste0("\nOverall elapsed time: ", round(tPTM[3]/60,1), " minutes"))
    rm(ptm, tPTM)
    gc()
  }
  # clean up unused variables
  rm(i,ii, cruise, ftle_ExtractAll, FTLE_Integrations,FTLEdatapath,data_Sub)
  # Remove the remnants of the Spatial Points DF conversion
  castSummary_FTLE <- castSummary_FTLE %>% dplyr::select(-c(coords.x1, coords.x2))   
  # Check for INF
  sapply(castSummary_FTLE, function(x) sum(is.infinite(x))>0)
  # convert INF values to na (only seems to exist in Max/Min)
  is.na(castSummary_FTLE) <- do.call(cbind,lapply(castSummary_FTLE, is.infinite))
  is.na(castSummary_FTLE) <- do.call(cbind,lapply(castSummary_FTLE, is.nan)) 
  ##Save extracted FTLE for later import
  save(castSummary_FTLE,file=paste0("./dataProcessed/castSummary_FTLE.RData"))
 
  
  #### Process GEBCO_2020 Bathy Data ####  
  # GEBCO_2020 Grid provides global coverage on a 15 arc-second grid (NOAA ETOPO1 was only 1 arc-minute)
  #   At 38 deg N, 15 arc-seconds is ~365m or .365km. 
  # This section imports a bathymmetry raster and extracts botDepth at Location
  bathyNCname <- "./dataRaw/BathymetryData/gebco_2020_n51.875711380012575_s21.52089312530333_w-132.52603158070204_e-108.35611507012007.nc"#file.choose()
  # Create a Raster Layer of HiRes Bathymetry
  bathy_raster <- raster(bathyNCname, varname = "elevation")
  
  ## Extract High-Resolution Bathymetry For Krill Data ##
  # Create an ID column to join later
  acousticsFTLE600cols$ID <- seq(1:length(acousticsFTLE600cols$Cruise))
  
  xy <- coordinates(cbind(acousticsFTLE600cols$Long, acousticsFTLE600cols$Lat))
  dataSP <- SpatialPointsDataFrame(xy, acousticsFTLE600cols, proj4string = CRS("+proj=longlat +ellps=WGS84 +no_defs"), 
                                   match.ID = FALSE)
  ### Extract Values ### 
  depth_Values <- raster::extract(bathy_raster,dataSP, cellnumbers = TRUE)
  colnames(depth_Values) <- c('rasterCellD','GebcoDepth')
  # combine into a single dataframe
  tempSP <- cbind(as.data.frame(dataSP),depth_Values)
  # Join Depth Data to cetacean_DF_FTLE_Bathy
  acousticsFTLE600cols_Bathy <- left_join(acousticsFTLE600cols,
                              dplyr::select(tempSP,ID,GebcoDepth),
                              by = c("ID")) %>% ungroup()
  # remove unused columns
  acousticsFTLE600cols_Bathy <- subset(acousticsFTLE600cols_Bathy, select = -c(ID)) 
  # Clean up
  rm(dataSP,xy,tempSP)
  rm(depth_Values,bathyNCname,bathy_raster)
  
  save(acousticsFTLE600cols_Bathy,file=paste0("./dataProcessed/acousticsFTLE600cols_Bathy.RData"))
  
 
#### Data Summary ####
  Cruise_Sum <- lineTimes %>% 
    full_join(CTD_Sum, by=c("Cruise","Line"))%>% 
    full_join(cetacean_Sum, by=c("Cruise","Line"))%>% 
    mutate(Year = year(start), Month = month(start))
  write_csv(Cruise_Sum,"./Output/Cruise_Sum.csv"   )
  
  ## Plot Acoustic Transects
  acoustics_DF %>%
    dplyr::select(Line, Lat,Long) %>%
    unique() %>%
    ggplot() +
    geom_sf(data = world) +
    coord_sf(xlim = c(min(acoustics_DF$Long,na.rm = TRUE)-.5, max(acoustics_DF$Long,na.rm = TRUE)+.5),
             ylim = c(min(acoustics_DF$Lat,na.rm = TRUE)-.5, max(acoustics_DF$Lat,na.rm = TRUE)+.5), expand = FALSE) +
    # add 200m contour
    geom_contour(data = bf,
                 aes(x=x, y=y, z=z),
                 breaks=c(-100,-200,-400,-800,-1600),
                 # breaks=c(-200),
                 size=c(0.4),
                 colour="darkgrey", show.legend = FALSE) +
    geom_text_contour(data = bf, aes(x=x, y=y,z = z),breaks=c(-100,-200,-400,-800,-1600),
                      # geom_text_contour(data = bf, aes(x=x, y=y,z = z),breaks=c(-200),
                      show.legend = FALSE, size = 2.2, alpha = .6, nudge_y = -.002) +
    geom_point(aes(Long,Lat, color = Line),size=1,show.legend = FALSE) +
    annotation_scale(location = "bl", width_hint = 0.5) +
    xlab("Longitude") +
    ylab("Latitude") +
    ggtitle("ACCESS Transect Lines (2012-2018)") +
    theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.5), panel.background = element_rect(fill = "aliceblue"))

  ## Plot Whale Sightings
  ggplot()+
    geom_sf(data = world) +
    coord_sf(xlim = c(min(acoustics_DF$Long,na.rm = TRUE)-.15, 
                      max(acoustics_DF$Long,na.rm = TRUE)+.15), 
             ylim = c(min(acoustics_DF$Lat,na.rm = TRUE)-.15, 
                      max(acoustics_DF$Lat,na.rm = TRUE)+.15), expand = FALSE) +
    # add 200m contour
    geom_contour(data = bf,
                 aes(x=x, y=y, z=z),
                 breaks=c(-100,-200,-400,-800,-1600),
                 # breaks=c(-200),
                 size=c(0.4),
                 colour="darkgrey", show.legend = FALSE) +
    geom_text_contour(data = bf, aes(x=x, y=y,z = z),breaks=c(-100,-200,-400,-800,-1600),
                      # geom_text_contour(data = bf, aes(x=x, y=y,z = z),breaks=c(-200),
                      show.legend = FALSE, size = 2.2, alpha = .6, nudge_y = -.002) +
    # geom_point(data=cetacean_DF[cetacean_DF$Cruise == "ACC1407",], aes(x=Long,y=Lat, group = SpCode,color = SpCode))+
    geom_point(data=cetacean_DF, aes(x=SightLong,y=SightLat, group = SpCode,color = SpCode))+
    facet_grid(~SpCode)
   
#### Save Workspace ####
save.image(file=paste0("Combined_alldata_",Sys.Date(), ".RData"))
  
  
