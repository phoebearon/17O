
########################### Supplement 7.  R script to calibrate raw triple oxygen IRMS data to the VSMOW-SLAP scale ###########################
# 
# Triple oxygen isotopes in the water cycle
# Author: Phoebe Aron
# Contact Phoebe Aron (paron@umich.edu) or Naomi Levin (nelevin@umich.edu) with questions and comments
# 
# This script shows how to calibrate raw triple oxygen IRMS data to the VSMOW-SLAP scale.
# This script generates three new dataframes, smow, slap, and usgs, with delta and delta prime values normalized to the VSMOW-SLAP scale.
#
# Before normalizing, raw data files should be 'cleaned up'. The 'cleaning' process is lab and instrument specific, but should 
#   consider whether individual analyses may be affected by memory and/or other analytical errors (e.g., errors during the fluorination or purification
#   process, problems with the mass spec, potentially compromised samples, etc.). We include Supplement 8 as an example of a clean file. 
#
# The example file (Supplement 8) is output from our Nu Perspective IRMS. The same basic principles of normalization work for data from any IRMS or laser absorption spectrometer. 

#   Column descriptions from example file (Supplement 8)
#   IPL num: this is a sequential, unique sample number of every triple oxygen isotope analysis in our lab
#   Type 1: sample identifier 
#   Type 2: sample identifier 
#   d33: m/z 33
#   d34: m/z 34
#   Date Time: timestamp of IRMS data reduction


## --------------------- Setup and Import Data --------------------- #

# define the path to the clean data file
#path.to.data <- "~User/paron/Desktop/"

# define files to import
files <- list.files(path=path.to.data, pattern="supplement_8") 

# import data
clean_data <- do.call("rbind", lapply(files, function(x) read.csv(paste(path.to.data, x, sep=''), stringsAsFactors = FALSE)))


## --------------------- VSMOW-SLAP Normalization --------------------- #

# define smow and slap data frames
smow <- subset(clean_data, Type.2=="SMOW")
slap <- subset(clean_data, Type.2=="SLAP")

#For SMOW and SLAP:  
# 1. define constants  
# 2. normalize to smow   
# 3. generate smow-slap transfer functions  
# 4. stretch to slap   
# 5. calculate dp (delta prime) and Dp17O values with smow-slap normalization  

# 1. DEFINE CONSTANTS
  d33.smow_ref <- mean(smow$d33) # observed
  d34.smow_ref <- mean(smow$d34) # observed

  lambda_ref <- 0.528 # Luz and Barkan, 2010

  d18O.slap <- -55.5 
  d17O.slap <- (exp(lambda_ref*log(d18O.slap/1000+1))-1)*1000 

# 2. NORMALIZE TO SMOW
  
  # observed SMOW
  smow$d33.smow_ref <- ((((smow$d33/1000)+1)/((d33.smow_ref/1000)+1))-1)*1000
  smow$d34.smow_ref <- ((((smow$d34/1000)+1)/((d34.smow_ref/1000)+1))-1)*1000

  #observed SLAP
  slap$d33.smow_ref <- ((((slap$d33/1000)+1)/((d33.smow_ref/1000)+1))-1)*1000
  slap$d34.smow_ref <- ((((slap$d34/1000)+1)/((d34.smow_ref/1000)+1))-1)*1000

# 3. GENERATE SMOW-SLAP TRANSFER FUNCTIONS
  
  #d17O linear model
  d17O.smow.accepted <- 0 #defined as 0
  d17O.smow.observed <- mean(smow$d33.smow_ref) #mean of observed smow
  d17O.slap.accepted <- d17O.slap #defined from slap d18O (-55.5) and lambda_ref
  d17O.slap.observed <- mean(slap$d33.smow_ref) #mean of observed slap
  lm.d17O.smow.slap <- lm(c(d17O.smow.accepted,d17O.slap.accepted) ~ c(d17O.smow.observed,d17O.slap.observed)) # linear model between smow and slap accepted and observed values
  d17O.smow.slap.slope <- lm.d17O.smow.slap$coefficients[2] #slope from linear model
  d17O.smow.slap.intercept <- lm.d17O.smow.slap$coefficients[1] #intercept from linear model

#d18O linear model 
  d18O.smow.accepted <- 0 
  d18O.smow.observed <- mean(smow$d34.smow_ref)
  d18O.slap.accepted <- d18O.slap 
  d18O.slap.observed <- mean(slap$d34.smow_ref) 
  lm.d18O.smow.slap <- lm(c(d18O.smow.accepted,d18O.slap.accepted) ~ c(d18O.smow.observed,d18O.slap.observed))
  d18O.smow.slap.slope <- lm.d18O.smow.slap$coefficients[2]
  d18O.smow.slap.intercept <- lm.d18O.smow.slap$coefficients[1]

# 4. STRETCH TO SLAP
  
  smow$d17O.slap <- smow$d33.smow_ref*d17O.smow.slap.slope
  smow$d18O.slap <- smow$d34.smow_ref*d18O.smow.slap.slope
  slap$d17O.slap <- slap$d33.smow_ref*d17O.smow.slap.slope
  slap$d18O.slap <- slap$d34.smow_ref*d18O.smow.slap.slope  

# 5. CALCULATE delta primes, Dp17O
  
  smow$dp17O.final <- log((smow$d17O.slap/1000)+1)*1000
  smow$dp18O.final <- log((smow$d18O.slap/1000)+1)*1000
  smow$D17O.final <- smow$dp17O.final - lambda_ref*smow$dp18O.final
  smow$D17O.per.meg <- smow$D17O.final * 1000

  slap$dp17O.final <- log((slap$d17O.slap/1000)+1)*1000
  slap$dp18O.final <- log((slap$d18O.slap/1000)+1)*1000
  slap$D17O.final <- slap$dp17O.final - lambda_ref*slap$dp18O.final
  slap$D17O.per.meg <- slap$D17O.final * 1000



## --------------------- Normalize Unknowns to the VSMOW-SLAP scale --------------------- #
  
#The example data file shows how to normalize USGS reference waters.  

  #Same steps as the standard data, but for unknowns. 
  # 1. normalize to smow  
  # 2. stretch to slap  
  # 3. calculate Dp17O
  # These steps  generate new columns with smow-slap normalized d18O, dp18O, d17O, dp17O, and Dp17O (in per mil and per meg)

# define usgs dataframe 
  usgs <- subset(clean_data, Type.1=="USGS")
  
# 1. Normalize to SMOW
  usgs$d33.smow_ref <- ((((usgs$d33/1000)+1)/((d33.smow_ref/1000)+1))-1)*1000
  usgs$d34.smow_ref <- ((((usgs$d34/1000)+1)/((d34.smow_ref/1000)+1))-1)*1000 

# 2. Stretch to SLAP
  usgs$d17O.slap <- usgs$d33.smow_ref*d17O.smow.slap.slope
  usgs$d18O.slap <- usgs$d34.smow_ref*d18O.smow.slap.slope

# 3. Calculate normalized values
  usgs$dp17O.final <- log((usgs$d17O.slap/1000)+1)*1000 # per mil
  usgs$dp18O.final <- log((usgs$d18O.slap/1000)+1)*1000 # per mil
  usgs$D17O.final <- usgs$dp17O.final - lambda_ref*usgs$dp18O.final # per mil
  usgs$D17O.per.meg <- usgs$D17O.final*1000 # per meg


