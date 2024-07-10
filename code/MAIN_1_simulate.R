# Description: 
# Runs the Gaussian puff atmospheric dispersion model to simulate methane
# concentrations at CMS sensor locations given wind data and potential source
# locations. 
# Author: William Daniels (wdaniels@mines.edu)
# Last Updated: December 2023

# Clear environment
rm(list = ls())

# Import necessary libraries
library(zoo)
library(lubridate)
library(foreach)
library(doParallel)
library(rstudioapi)

if (commandArgs()[1] == "RStudio"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

# START USER INPUT
#---------------------------------------------------------------------------

# Set path to simulation configuration file
config.file.dir <- '../input_data/simulation_config.txt'

# END OF USER INPUT - NO MODIFICATION NECESSARY BELOW THIS POINT
#---------------------------------------------------------------------------



# STEP 1: READ IN CONFIG FILE AND SET UP PARAMETERS AND DIRECTORY STRUCTURE
#---------------------------------------------------------------------------

# Read in config file
config <- read.table(config.file.dir)
config <- strsplit(config[,1], "=")

# Parse out config file
parameters <- sapply(config, function(X) X[1])
values <- sapply(config, function(X) X[2])

# Get parameter values
num.cores.to.use <- as.numeric(values  [parameters == "num.cores.to.use"])
dt               <- as.numeric(values  [parameters == "dt"])
cutoff.t         <- as.numeric(values  [parameters == "cutoff.t"])
ignore.dist      <- as.numeric(values  [parameters == "ignore.dist"])
chunk.size       <- as.numeric(values  [parameters == "chunk.size"])
emission.rate    <- as.numeric(values  [parameters == "emission.rate"])
run.mode         <- as.character(values[parameters == "run.mode"])
tz               <- as.character(values[parameters == "tz"])

# Get directories
raw.sensor.observations.path <-     as.character(values[parameters == "raw.sensor.observations.path"])
source.locations.path <-            as.character(values[parameters == "source.locations.path"])
sensor.locations.path <-            as.character(values[parameters == "sensor.locations.path"])
output.file.path <-                 as.character(values[parameters == "output.file.path"])
helper.distance.conversions.path <- as.character(values[parameters == "helper.distance.conversions.path"])
helper.gpuff.function.path <-       as.character(values[parameters == "helper.gpuff.function.path"])

# Source helper files which contain helper functions 
source(helper.distance.conversions.path)
source(helper.gpuff.function.path)

# Get start and stop times for simulation
start.time <- as_datetime(as.character(values[parameters == "start.time"]), tz = tz)
end.time   <- as_datetime(as.character(values[parameters == "end.time"]),   tz = tz)

# Start code timer
code.start.time <- Sys.time()



# STEP 2: READ IN CMS SENSOR DATA
#---------------------------------------------------------------------------

# Read in sensor data csv
raw.data <- read.csv(raw.sensor.observations.path)

# Set correct time zone
raw.data$time <- as_datetime(as_datetime(raw.data$time), tz = tz)

# Mask in data within start and end times
this.mask <- raw.data$time >= start.time & raw.data$time <= end.time
raw.data <- raw.data[this.mask, ]


# STEP 3: READ IN SITE GEOMETRY AND SETUP DATA STRUCTURE TO EXPORT
#---------------------------------------------------------------------------

# Read in sensor location csv
sensor.locs <- suppressWarnings(read.csv(sensor.locations.path))
sensor.locs <- sensor.locs[order(sensor.locs$name), ]
n.r <- nrow(sensor.locs)

# Read in source location csv
source.locs <- suppressWarnings(read.csv(source.locations.path))
n.s <- nrow(source.locs)

# Initialize data structure to store simulation output
data.to.save <- vector(mode = "list", length = 4 + n.s)
names(data.to.save) <- c("times", "obs", "WD", "WS", source.locs$name)


# STEP 4: COMBINE WIND DATA ACROSS SENSORS THAT MEASURE WIND
#---------------------------------------------------------------------------

# Figure out which units collect wind data
wind.unit <- vector(length = length(n.r))
for (r in 1:n.r){
  these.raw.data <- raw.data[raw.data$name == sensor.locs$name[r], ]
  wind.unit[r] <- any(!is.na(these.raw.data$wind.speed))
}


# create a sequence of minutes from start time to end time
min.seq <- seq(min(raw.data$time, na.rm = T), 
               max(raw.data$time, na.rm = T),
               by = "min")

# Fire up the parallel cluster
cl <- makeCluster(num.cores.to.use)
registerDoParallel(cl)

# Loop through time chunks
wind.list <- foreach(t = 1:length(min.seq)) %dopar% {
  
  these.wd <- these.ws <- vector(length = sum(wind.unit))
  ind <- 0
  for (r in which(wind.unit)){
    ind <- ind + 1
    this.mask <- raw.data$name == sensor.locs$name[r] & raw.data$time == min.seq[t]
    
    # if sum(this.mask) == 1, all is good, grab the observation from that sensor
    # if sum(this.mask) > 1, there are two obs at the same time stamp... weird... grab just the first
    # if sum(this.mask) == 0, no obs at that time stamp, use an NA
    if (sum(this.mask) > 0){
      
      these.wd[ind] <- raw.data$wind.direction[this.mask][1]
      these.ws[ind] <- raw.data$wind.speed[this.mask][1]
    } else {
      these.wd[ind] <- NA
      these.ws[ind] <- NA
    }
  }
  
  to.add <- list(WS = these.ws, WD = these.wd)
  to.add
}

stopCluster(cl)

# Parse output of foreach loop
out <- do.call(rbind, wind.list)

# Create matrix of wind speed and wind direction observations
# Rows are different minutes and columns are different sensors
WS <- do.call(rbind, out[,1])
WD <- do.call(rbind, out[,2])

# Gap fill the WS
WS <- na.approx(WS, maxgap = 5, na.rm = F)

# Take median of WS across sensors
WS <- apply(WS, 1, median, na.rm = T)

# Convert WD to radians and compute horizontal and vertical components
WD <- WD * pi/180
WD.x <- cos(WD)
WD.y <- sin(WD)

# Gap fill x and y components separately
WD.x <- na.approx(WD.x, maxgap = 5, na.rm = F)
WD.y <- na.approx(WD.y, maxgap = 5, na.rm = F)

# Take median across sensors of x and y components separately
WD.x <- apply(WD.x, 1, median, na.rm = T)
WD.y <- apply(WD.y, 1, median, na.rm = T)

# Recombine into an angle and clean up
WD <- atan2(WD.y, WD.x) * 180 / pi
WD <- ifelse(WD < 0, WD + 360, WD)

# Convert wind direction from clockwise with 0 deg at north to
# counterclockwise with 0 deg at east (traditional angle definition)
WD <- 90 - WD
WD <- ifelse(WD < 0, WD + 360, WD)

# Convert wind direction from direction wind is coming from to
# direction wind is going
WD <- 180 + WD
WD <- ifelse(WD >= 360, WD - 360, WD)

# Convert wind direction to radians
WD <- WD * pi/180


# STEP 5: MOVE DATA INTO TIME-ALIGNED TABLE WITH SENSORS AS COLUMNS
#---------------------------------------------------------------------------

# Set up data table to hold cleaned data
data <- matrix(NA, nrow = length(min.seq), ncol = n.r+3)
colnames(data) <- c(sensor.locs$name, "times", "WS", "WD")
data <- as.data.frame(data)
data$times <- min.seq
data$WS <- WS
data$WD <- WD

# Fire up the parallel cluster
cl <- makeCluster(num.cores.to.use)
registerDoParallel(cl)

# Loop through minute sequence and add data to the cleaned data table
row.list <- foreach(t = 1:length(min.seq)) %dopar% {
  
  this.row <- vector(length = n.r)
  
  for (r in 1:n.r){
    this.mask <- raw.data$name == sensor.locs$name[r] & raw.data$time == min.seq[t]
    
    # if sum(this.mask) == 1, all is good, grab the observation from that sensor
    # if sum(this.mask) > 1, there are two obs at the same time stamp... weird... grab just the first
    # if sum(this.mask) == 0, no obs at that time stamp, use an NA
    if (sum(this.mask) > 0){
      this.row[r] <- raw.data$methane[this.mask][1]
    } else {
      this.row[r] <- NA
    }
  }
  
  this.row
}

stopCluster(cl)

# Combine output of parallel foreach loop into one matrix
data[1:nrow(data), 1:n.r] <- do.call(rbind, row.list)

# Remove observations that still have no wind data
data <- data[!is.na(data$WS),]

# Pull out times
times <- data$times

# Pull out just the methane observations
obs <- data[, 1:n.r]

# Clean up
rm(raw.data, this.mask, wind.unit, WS, WD, row.list,
   wind.list, WD.x, WD.y, config, min.seq)


# STEP 6: LOOP THROUGH POTENTIAL SOURCES AND SIMULATE CONCENTRATION VALUES AT RECEPTORS
#---------------------------------------------------------------------------
for (s in 1:n.s){
  
  print(paste0("source: ", s, "/", n.s))
  
  # Get source location
  source.lon <- as.numeric(source.locs$lon[s])
  source.lat <- as.numeric(source.locs$lat[s])
  
  # Convert source lat/lon to utm
  source.x <- as.numeric(latlon.to.utm(source.lat, source.lon)[1])
  source.y <- as.numeric(latlon.to.utm(source.lat, source.lon)[2])
  source.z <- source.locs$height[s]
  
  # Get receptor locations
  recept.lon <- sensor.locs$lon
  recept.lat <- sensor.locs$lat
  z.r <- sensor.locs$height # receptor heights [m]
  
  # Convert receptor lat/lon to utm
  x.r <- y.r <- vector(length = n.r)
  for (j in 1:n.r){
    x.r[j] <- as.numeric(latlon.to.utm(recept.lat[j], recept.lon[j])[1])
    y.r[j] <- as.numeric(latlon.to.utm(recept.lat[j], recept.lon[j])[2])
  }
  
  # Get wind data
  WA <- data$WD
  WS <- data$WS
  
  # Set emission rate for all time steps
  Q.truth <- rep(emission.rate, length(times)) # [kg/s]
  
  # Compute number of time intervals based on simulation frequency set earlier
  n.low.res.ints <- length(WA) # number of one-minute data points
  n.ints <- (n.low.res.ints-1) * (60/dt) + 1 # number of intervals at simulation frequency
  
  # Get x and y components of wind angle
  WA.x <- cos(WA)
  WA.y <- sin(WA)
  
  # Interpolate one-minute-data to get it on simulation frequency
  est.on.fine.grid <- function(y, dt){
    out <- y[1]
    for (i in 2:length(y)){
      interp.y <- seq(y[i-1], y[i], length.out=60/dt)
      out <- c(out, interp.y)
    }
    return(out)
  }
  
  # Increase frequency of wind angle
  WA.x <- est.on.fine.grid(WA.x, dt)
  WA.y <- est.on.fine.grid(WA.y, dt)
  
  # Bring back to angle
  WA <- atan2(WA.y, WA.x) 
  WA <- ifelse(WA < 0,     WA + (2*pi), WA)
  WA <- ifelse(WA >= 2*pi, WA - (2*pi), WA)
  
  # Increase frequency of WS and Q
  WS      <- est.on.fine.grid(WS, dt)
  Q.truth <- est.on.fine.grid(Q.truth, dt)
  
  # Get x and y components of wind speed at simulation frequency
  WS.x <- cos(WA) * WS
  WS.y <- sin(WA) * WS
  
  # Compute number of time chunks based on "chunk.size"
  n.chunks <- ceiling(n.ints / (chunk.size*(60/dt)))
  
  # Fire up the parallel cluster
  cl <- makeCluster(num.cores.to.use)
  registerDoParallel(cl)
  
  # Loop through time chunks
  big.C.list <- foreach(h = 1:n.chunks) %dopar% {
    
    # Need to reinitialize libraries used within foreach loop..
    library(lubridate)
    
    # Standard chunk size at simulation frequency
    standard.chunk.size <- chunk.size*60/dt
    
    # First index of this chunk at simulation frequency
    this.chunk.start <- 1 + (h-1) * standard.chunk.size
    
    # Last index of this chunk at simulation frequency 
    this.chunk.end   <- min(h * standard.chunk.size,
                            n.ints)
    
    # Last index of the extended chunk to capture stragglers
    ext.chunk.end    <- min(h * standard.chunk.size + cutoff.t*60/dt - 1,
                            n.ints)
    
    # Mask for this chunk in simulation frequency
    this.chunk.mask <- this.chunk.start:this.chunk.end
    
    # Mask for extended chunk to capture stragglers in simulation frequency
    ext.chunk.mask  <- this.chunk.start:ext.chunk.end
    
    # Initialize matrices to hold puff locations and distance traveled
    puff.x.locs <- puff.y.locs <- total.dist <- matrix(NA, 
                                                       nrow = cutoff.t*(60/dt), 
                                                       ncol = length(ext.chunk.mask))
    
    # Compute first location and distance based on source location
    puff.x.locs[1,] <- source.x 
    puff.y.locs[1,] <- source.y 
    total.dist[1,]  <- 0
    
    # Fill in remainder of location and distance matrices using previous locations
    for (p in 1:ncol(puff.x.locs)){ # Loop through puffs
      for (t in 2:nrow(puff.x.locs)){ # Loop through time steps that this puff is alive
        
        # Compute movement
        if (run.mode == "dynamic"){
          
          x.movement <- WS.x[ext.chunk.mask][t+p-1] * dt
          y.movement <- WS.y[ext.chunk.mask][t+p-1] * dt
          
        } else if (run.mode == "constant"){
          
          x.movement <- WS.x[ext.chunk.mask][p] * dt
          y.movement <- WS.y[ext.chunk.mask][p] * dt
          
        } else {
          print("Invalid run mode")
        }
        
        # Compute new location and distance
        puff.x.locs[t,p] <- puff.x.locs[t-1, p] + x.movement
        puff.y.locs[t,p] <- puff.y.locs[t-1, p] + y.movement
        total.dist[t,p]  <- total.dist[t-1, p]  + sqrt(x.movement^2 + y.movement^2)
      }
    }
    
    if (run.mode == "constant"){
      
      # Initialize list to hold stability classes at each time step.
      # Puffs created at that time step will keep that stability class 
      # until they are eliminated
      stab.classes <- vector(mode = "list", length = length(ext.chunk.mask))
      
      # Loop through times, get time and WS, get stab class, save stab class
      for (j in 1:length(ext.chunk.mask)){
        
        # Get time
        time.to.use <- times[1] + seconds((h-1) * standard.chunk.size + dt*j)
        
        # Correction that only applies when there is a one hour gap in times due to day light savings.
        # This correction sets the time to be just after day light savings occurs and prevents NAs.
        # Note that the year is hard-coded, but this does not matter, as the time of day is the 
        # only aspect of the datetime object that is used to determine stability class.
        time.to.use <- as_datetime(ifelse(is.na(time.to.use),"2022-03-13T05:00:00",time.to.use), tz = tz)
        
        # Get wind speed
        WS.to.use <- WS[ext.chunk.mask][j]
        
        # Get stability class based on time and wind speed
        stab.classes[[j]] <- get.stab.class(WS.to.use, time.to.use)
      }
    }
    
    # Initialize matrix to hold concentration values for each receptor for each time step
    C <- matrix(NA, nrow = length(this.chunk.mask), ncol = n.r)
    
    # Initialize variable to hold previous time used to get stability class
    previous.time <- c()
    
    # Iterate through time steps at simulation frequency 
    for (j in 1:length(this.chunk.mask)){
      
      if (run.mode == "dynamic"){
        
        # Get time to use for stability class
        time.to.use <- times[1] + seconds((h-1) * standard.chunk.size + dt*j)
        
        # If time is NA (happens during daylight savings time transition) use previous time
        time.to.use <- as_datetime(ifelse(is.na(time.to.use), previous.time, time.to.use), tz = tz)
        
        # Get the stability class based on wind speed and time of day
        stab.class <- get.stab.class(U = WS[this.chunk.mask][j], 
                                     time = time.to.use)
        
        # Reset previous time
        previous.time <- time.to.use
        
      }
      
      # Vector to hold concentration predictions at the jth time step
      this.C <- rep(0, n.r)
      
      # Loop through puffs in existence at the jth time step
      for (k in 1:(cutoff.t*60/dt)){
        
        # Compute matrix indices for the kth puff in the pre-computed location matrices
        col.it <- j-k+1
        row.it <- k
        
        # End loop over puffs after all puffs are accounted for.
        # Will only be triggered early on when the number of puffs in existence is less than cutoff.t*60/dt
        if (col.it < 1) {break}
        
        # Get location of the kth puff and the total distance it has traveled
        this.puff.x <- puff.x.locs[row.it, col.it]
        this.puff.y <- puff.y.locs[row.it, col.it]
        this.total.dist <- total.dist[row.it, col.it]
        
        if (run.mode == "constant"){
          stab.class <- stab.classes[[col.it]]
        }
        
        # Compute distance between this puff and all receptors
        distances <- sqrt((x.r-this.puff.x)^2 + (y.r-this.puff.y)^2)
        
        # Skip this puff if smallest distance is greater than ignore.dist
        if (min(distances) > ignore.dist) {next}
        
        # Compute the concentration at each receptor location from the kth puff
        this.C <- this.C + dt * gpuff(Q = Q.truth[this.chunk.mask][j], 
                                      stab.class = stab.class,
                                      x.p = this.puff.x,
                                      y.p = this.puff.y,
                                      x.r.vec = x.r, 
                                      y.r.vec = y.r, 
                                      z.r.vec = z.r,
                                      total.dist = this.total.dist,
                                      H = source.z, 
                                      U = WS[this.chunk.mask][j])
        
      } # end loop over puffs
      
      C[j,] <- this.C
      
    } # end loop through time steps
    
    if (h != n.chunks){
      
      # Initialize matrix to hold concentration values for each receptor for each time step
      C.to.pass <- matrix(NA, nrow = ext.chunk.end - this.chunk.end, ncol = n.r)
      
      # Initialize variable to hold previous time used to get stability class
      previous.time <- c()
      
      this.it <- 1
      
      # Iterate through extra time steps at simulation frequency to get bottom right corner
      for (j in seq(length(this.chunk.mask)+1,length(ext.chunk.mask))){
        
        if (run.mode == "dynamic"){
          
          # Get time to use for stability class
          time.to.use <- times[1] + seconds((h-1) * standard.chunk.size + dt*j)
          
          # If time is NA (happens during daylight savings time transition) use previous time
          time.to.use <- as_datetime(ifelse(is.na(time.to.use), previous.time, time.to.use), tz = tz)
          
          # Get the stability class based on wind speed and time of day
          stab.class <- get.stab.class(U = WS[ext.chunk.mask][j], 
                                       time = time.to.use)
          
          # Reset previous time
          previous.time <- time.to.use
          
        }
        
        # Vector to hold concentration predictions at the jth time step
        this.C <- rep(0, n.r)
        
        # Loop through puffs in existence at the jth time step
        for (k in 1:(cutoff.t*60/dt)){
          
          # Compute matrix indices for the kth puff in the pre-computed location matrices
          col.it <- j-k+1
          row.it <- k
          
          if (col.it > length(this.chunk.mask)) {next}
          
          # End loop over puffs after all puffs are accounted for.
          # Will only be triggered early on when the number of puffs in existence is less than cutoff.t*60/dt
          if (col.it < 1) {break}
          
          # Get location of the kth puff and the total distance it has traveled
          this.puff.x <- puff.x.locs[row.it, col.it]
          this.puff.y <- puff.y.locs[row.it, col.it]
          this.total.dist <- total.dist[row.it, col.it]
          
          if (run.mode == "constant"){
            stab.class <- stab.classes[[col.it]]
          }
          
          # Compute distance between this puff and all receptors
          distances <- sqrt((x.r-this.puff.x)^2 + (y.r-this.puff.y)^2)
          
          # Skip this puff if smallest distance is greater than ignore.dist
          if (min(distances) > ignore.dist) {next}
          
          # Compute the concentration at each receptor location from the kth puff
          this.C <- this.C + dt * gpuff(Q = Q.truth[ext.chunk.mask][j], 
                                        stab.class = stab.class,
                                        x.p = this.puff.x,
                                        y.p = this.puff.y,
                                        x.r.vec = x.r, 
                                        y.r.vec = y.r, 
                                        z.r.vec = z.r,
                                        total.dist = this.total.dist,
                                        H = source.z, 
                                        U = WS[ext.chunk.mask][j])
          
        } # end loop over puffs
        
        C.to.pass[this.it,] <- this.C
        this.it <- this.it + 1
        
      } # end loop through time steps
      
    } else {
      
      C.to.pass <- NA
    }
    
    # Save pertinent data
    to.save <- list(C = C,
                    C.to.pass = C.to.pass)
    to.save
    
    
    
  } # end loop through time chunks
  
  stopCluster(cl)
  
  
  # Add the bottom right corner to the next time chunk
  if (n.chunks > 1){
    for (h in 2:n.chunks){
      correction <- big.C.list[[h-1]]$C.to.pass 
      original.C <- big.C.list[[h]]$C[1:nrow(correction), ]
      big.C.list[[h]]$C[1:nrow(correction), ] <- original.C + correction
    }
  }
  
  # Extract data from list object
  big.C.to.extract <- lapply(big.C.list, function(X) X$C)
  big.C <- do.call(rbind, big.C.to.extract)
  
  # Average the concentration predictions back up to a one-minute resolution
  C.avg <- matrix(NA, nrow = length(times), ncol = n.r)
  
  # First entry corresponds to last second of the first minute
  # This is because we define the minute time stamp to be the end of the minute
  # rather than the beginning of the minute
  C.avg[1, ] <- big.C[1, ]
  
  # Fill in the rest of the matrix
  for (j in 2:nrow(C.avg)){
    
    this.mask <- seq((j-2)*(60/dt)+2,
                     (j-1)*(60/dt)+1)
    
    if (length(this.mask) > 1){
      C.avg[j, ] <- apply(big.C[this.mask, ], 2, mean)
    } else {
      C.avg[j, ] <- big.C[j,]
    }
  }
  
  
  # Store simulation data
  data.to.save[[1]] <- times
  data.to.save[[2]] <- as.data.frame(obs)
  data.to.save[[3]] <- data$WD
  data.to.save[[4]] <- data$WS
  tmp <- as.data.frame(C.avg)
  colnames(tmp) <- sensor.locs$name
  data.to.save[[s+4]] <- tmp
  
} # end loop over sources

# Save simulation data
saveRDS(data.to.save, output.file.path)

# End code timer
code.stop.time <- Sys.time() 

# Print wall clock execution time
difftime(code.stop.time, code.start.time, units = "mins")
