# Description: 
# Performs emission event detection, localization, and quantification
# using output from the Gaussian puff dispersion model and CMS concentration
# observations.
# Author: William Daniels (wdaniels@mines.edu)
# Last Updated: December 2023

# Clear environment
if(!is.null(dev.list())){dev.off()}
rm(list = ls())

# Import necessary libraries
library(lubridate)
library(zoo)
library(rstudioapi)

if (commandArgs()[1] == "RStudio"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}


# START USER INPUT
#---------------------------------------------------------------------------

# Set path to simulation configuration file
config.file.dir <- '../input_data/DLQ_config.txt'

# END OF USER INPUT - NO MODIFICATION NECESSARY BELOW THIS POINT
#---------------------------------------------------------------------------



# STEP 1: READ IN CONFIG FILE AND SET UP PARAMETERS AND DIRECTORY STRUCTURE
#---------------------------------------------------------------------------

# Read in config file
config <- suppressWarnings(read.table(config.file.dir))
config <- strsplit(config[,1], "=")

# Parse out config file
parameters <- sapply(config, function(X) X[1])
values <- sapply(config, function(X) X[2])

# Get parameter values
gap.time         <- as.numeric(values  [parameters == "gap.time"])
length.threshold <- as.numeric(values  [parameters == "length.threshold"])

# Get directories
simulation.data.path <-            as.character(values[parameters == "simulation.data.path"])
output.file.path <-                as.character(values[parameters == "output.file.path"])
helper.spike.detection.alg.path <- as.character(values[parameters == "helper.spike.detection.alg.path"])

# Source helper files which contain helper functions 
source(helper.spike.detection.alg.path)

# Read in simulation data
data <- readRDS(simulation.data.path)

# Pull out sensor observations and replace NA's that are not on edge of 
# the time series with interpolated values
obs <- na.approx(data$obs, na.rm = F)

# Number of sensors
n.r <- ncol(obs)

# Pull out time stamps of observations and simulations
times <- data$times

# Pull out the simulation predictions
sims <- data[5:length(data)]

# Grab source info
n.s <- length(sims) 
source.names <- names(sims)

# Define function to create a logarithmic spaced sequence (used later)
lseq <- function(from, to, length.out) {
  exp(seq(log(from), log(to), length.out = length.out))
}





# STEP 2: REMOVE BACKGROUND FROM OBSERVATIONS
#---------------------------------------------------------------------------

# Skip sensors that have only NA values
to.use <- which(apply(obs, 2, function(X) !all(is.na(X))))

# Loop through sensor units
for (j in to.use){
  
  # Grab observations from this sensor
  this.raw.obs <- obs[,j]
  
  # Remove the NA's at the beginning and end of the time series
  # (some sensors start late or end early)
  to.keep <- !is.na(this.raw.obs)
  trimmed.obs <- this.raw.obs[to.keep]
  trimmed.times <- times[to.keep]
  
  # Flag spikes
  spikes <- find.spikes(obs = trimmed.obs,
                        times = trimmed.times,
                        amp.threshold = 0.75,
                        make.plot = F)
  
  # Add points immediately before and after the spike to the spike mask
  # This better captures the full event, as there is some delay between the 
  # start of an emission and when the sensors first see an enhancement
  for (i in na.omit(unique(spikes$events))){
    min.ind <- min(which(spikes$events == i))
    max.ind <- max(which(spikes$events == i))
    spikes$events[c(max(c(min.ind-1, 1)), min(c(max.ind+1,nrow(spikes))))] <- i
  }
  
  # Pull out integers that uniquely identify an event
  event.nums <- na.omit(unique(spikes$events))
  
  # If there are at least two events
  if (length(event.nums) > 1){
    
    # Loop through spikes and combine spikes that are separated by less than gap.time minutes
    for (i in 2:length(event.nums)){
      
      # Mask in this spike and last spike
      this.spike <- spikes$events == event.nums[i]
      previous.spike <- spikes$events == event.nums[i-1]
      
      # Clean up
      this.spike[is.na(this.spike)] <- F
      previous.spike[is.na(previous.spike)] <- F
      
      # Get start time of current spike and end time of previous spike
      this.spike.start.time <- spikes$time[this.spike][1]
      previous.spike.end.time <- spikes$time[previous.spike][length(spikes$time[previous.spike])]
      
      # Compute time difference
      time.diff <- difftime(this.spike.start.time, previous.spike.end.time, units = "mins")
      
      # Check gap
      if (minutes(time.diff) < minutes(gap.time)){
        spikes$events[this.spike] <- event.nums[i-1]
        event.nums[i] <- event.nums[i-1]
      }
    }
  }
  
  # Grab the new event numbers after combining events in the previous for loop
  event.nums <- na.omit(unique(spikes$events))
  
  # If there are any events
  if (length(event.nums) > 0){
    
    # Loop through events and (1) fill in any gaps between events with the correct
    # event number, and (2) estimate background as mean of first and last spike 
    # points, which occur before the sharp increase and after the spike has 
    # returned to return.threshold percent of the max value
    for (i in 1:length(event.nums)){
      
      # Fill in gaps
      first.ob <- min(which(spikes$events == event.nums[i]))
      last.ob <- max(which(spikes$events == event.nums[i]))
      this.mask <- first.ob:last.ob
      spikes$events[this.mask] <- event.nums[i]
      
      # Estimate background using observations before and after spike
      # This is just the first and last observation within the spike mask,
      # since we already added the points before and after the spike to the 
      # spike mask earlier on
      b.left <- trimmed.obs[first.ob]
      b.right <- trimmed.obs[last.ob]
      b <- mean(c(b.left, b.right))
      
      # Remove background from this spike
      trimmed.obs[this.mask] <- trimmed.obs[this.mask] - b
    }
    
  }
  
  # Remove background from all non-spike data. Since by definition these points
  # are not in an event, their concentration value is directly taken as the 
  # background estimate. Hence removing background is just setting the value to zero
  trimmed.obs[is.na(spikes$events)] <- 0
  
  # Set any negative values to zero. Negative value would arise if the 
  # background estimate was too large (greater than actual concentration value)
  trimmed.obs[trimmed.obs < 0] <- 0
  
  # Save background removed data
  obs[to.keep, j] <- trimmed.obs
}





# STEP 3: EVENT DETECTION
#---------------------------------------------------------------------------

# Compute maximum concentration value across sensors for each minute
# Note that since all non-spike observations have been set to zero, any 
# non-zero value in the max.obs time series should be considered an event 
# (and hence no need to run the spike detection algorithm again)
max.obs <- apply(obs, 1, max, na.rm = T)

# Create data frame with time steps and event mask
spikes <- data.frame(time = times, events = max.obs > 0)

# Find gaps between events that are shorter than gap.time and turn them into events
#   to.replace holds the indices that need to be switched from F to T
#   first.gap is an indicator for the first gap (which should not be replaced,
#   regardless of length)
#   false.seq holds the indices of each sequence of FALSEs (non-events)
to.replace <- c()
first.gap <- T
false.seq <- c()

# Loop through times
for (i in 1:length(times)){
  
  # If not a spike, add index to false.sequence
  if (!spikes$events[i]){
    false.seq <- c(false.seq, i)
    
    # Otherwise, check length of false sequence, if greater than gap time,
    # save those indices to replace later
  } else if (spikes$events[i]){
    
    if (length(false.seq) <= gap.time & !first.gap){
      to.replace <- c(to.replace, false.seq)
    }
    
    first.gap <- F
    false.seq <- c()
  }
}

# Replace gaps shorter than gap.time with T (meaning they are in an event)
spikes$events[to.replace] <- T

# Now we replace the T/F with an integer to distinguish between events
# Start by replacing F with NA and T with zero
spikes$events[!spikes$events] <- NA
spikes$events[spikes$events] <- 0

# Get indices of spikes
spike.points <- which(!is.na(spikes$events))
count <- 0

# Loop through spike points, if last point was not a spike, increase counter
for (i in spike.points){
  
  if (is.na(spikes$events[i-1])){
    count <- count + 1
    spikes$events[i] <- count
  } else {
    spikes$events[i] <- count
  }
}

# Get integers that uniquely define the different events
event.nums <- na.omit(unique(spikes$events))

# Filter events by the length threshold
for (i in 1:length(event.nums)){
  this.spike <- which(spikes$events == event.nums[i])
  
  if (length(this.spike) < length.threshold){
    spikes$events[this.spike] <- NA
  }
}

# Grab event number again after filtering by length
event.nums <- na.omit(unique(spikes$events))

# Number of events that we will perform localization / quantification on
n.ints <- length(event.nums)



# STEP 4: COMPUTE ALIGNMENT METRIC
#---------------------------------------------------------------------------

# Matrix to hold alignment metric for each event and potential source
metrics <- matrix(NA, nrow = n.ints, ncol = n.s)

# Loop through events
for (t in 1:n.ints){
  
  # Mask in this event
  this.mask <- seq(min(which(spikes$events == event.nums[t])),
                   max(which(spikes$events == event.nums[t])))
  
  # Loop through potential sources
  for (s in 1:n.s){
    
    # Grab simulation predictions from this source
    preds <- as.matrix(sims[[s]])
    
    # Initialize variables that will hold the predictions and observations
    # that we will compare to make the localization estimate
    all.obs.to.compare <- all.preds.to.compare <- c()
    
    # Loop through sensors
    for (r in 1:n.r){
      
      # Mask in only this event and this sensor
      obs.int <- obs[this.mask, r]
      preds.int <- preds[this.mask, r] 
      
      # Only compare non-NA observations
      to.compare <- !is.na(obs.int) 
      
      # Save predictions and observations to compare
      all.obs.to.compare <- c(all.obs.to.compare, obs.int[to.compare])
      all.preds.to.compare <- c(all.preds.to.compare, preds.int[to.compare])
      
    } # end loop through sensors
    
    # At least one non-NA observation is required to compute metric
    if (!all(is.na(all.obs.to.compare))){
      
      # Rename for brevity
      x <- all.preds.to.compare
      y <- all.obs.to.compare
      
      # Replace tiny values with zero
      x <- ifelse(x < 1e-30, 0, x)
      y <- ifelse(y < 1e-30, 0, y)
      
      # If there are nonzero values to compare, compute correlation
      if (all(x == 0) | all(y==0)){
        metrics[t,s] <- NA
      } else {
        metrics[t,s] <- cor(x, y, use = "complete.obs", method = "pearson")  
      }
      
    } # end if to check if there is at least one non-NA observation
  } # end loop through sources
} # end loop through events


# Convert NAs to zeros
metrics[is.na(metrics)] <- 0


# STEP 5: COMPUTE LOCALIZATION AND QUANTIFICATION ESTIMATES
#---------------------------------------------------------------------------

# Vectors to hold localization and quantification results
rate.est.all.events <- loc.est.all.events <-
  error.lower.all.events <- error.upper.all.events <- vector(length = n.ints)

all.preds.to.compare.all.events <- all.obs.to.compare.all.events <- c()

# Loop through events
for (t in 1:n.ints){
  
  print(paste0(t, "/", n.ints))
  
  # Get metric values for this event
  these.metrics <- metrics[t,]
  
  # Get source name corresponding to largest metric value
  loc.est.all.events[t] <- source.names[which.max(these.metrics)]
  
  # Mask in this event
  this.mask <- seq(min(which(spikes$events == event.nums[t])),
                   max(which(spikes$events == event.nums[t])))
  
  # Get predictions from most likely source and observations during this event
  event.preds <- sims [[which.max(these.metrics)]][this.mask, ]
  event.obs   <- obs  [this.mask, ]
  event.times <- times[this.mask]
  
  # Vectors to hold preds and obs when both are in a spike
  all.preds.to.compare <- all.obs.to.compare <- c()
  
  # Loop through sensors and save when both obs and preds are in a spike
  for (r in 1:n.r){
    
    # Get predictions and observations for this sensor
    these.preds <- event.preds[,r]
    these.obs <- event.obs[,r]
    
    if (all(is.na(these.obs))){ next }
    
    # Find spikes in both predictions and observations.
    # We will only use data in which both a observation and prediction is in a spike
    # to estimate the emission rate. This reduces impact of forward model inadequacies.
    preds.spikes <- find.spikes(event.times, these.preds, amp.threshold = 1, make.plot = F)
    obs.spikes   <- find.spikes(event.times, these.obs,   amp.threshold = 1, make.plot = F)
    
    # Mask for times in which both preds and obs are in a spike
    both.in.spike.mask <- !is.na(preds.spikes$events) & !is.na(obs.spikes$events)
    
    # If there are any such time steps..
    if (sum(both.in.spike.mask) > 0){
      
      # Save data to compare
      all.preds.to.compare <- c(all.preds.to.compare, these.preds[both.in.spike.mask])
      all.obs.to.compare   <- c(all.obs.to.compare,   these.obs  [both.in.spike.mask])
    } 
    
  } # End loop over sensors
  
  
  # If there is enough data to compare, compute rate estimate and 
  # percent difference between obs and preds (will be used for UQ)
  if (length(all.preds.to.compare) > 4){
    
    # Save obs and preds so that we can check that they are on the same order of magnitude later
    all.preds.to.compare.all.events <- c(all.preds.to.compare.all.events, all.preds.to.compare)
    all.obs.to.compare.all.events <- c(all.obs.to.compare.all.events, all.obs.to.compare)
    
    # Number of times to sample from predictions and observations
    n.samples <- 1000
    
    # Vector to hold emission rate estimate for each sample
    q.vals <- vector(length = n.samples)
    
    # Loop through samples
    for (o in 1:n.samples){
      
      # Sample from predictions and observations
      this.sample <- sample.int(n = length(all.preds.to.compare), 
                                size = length(all.preds.to.compare)/2,
                                replace = T)
      
      # Define a grid of emission rate values to optimize over
      q.grid <- lseq(0.0001, 3000, length.out = 2000)
      
      # Vector to hold sum of squared errors
      sse <- vector(length = length(q.grid))
      
      # Loop through possible emission rates
      for (z in 1:length(q.grid)){
        
        # Scale predictions by this emission rate
        qxp <- q.grid[z] * all.preds.to.compare
        
        # Compute root mean square error
        sse[z] <- sqrt( mean( (all.obs.to.compare[this.sample] - qxp[this.sample])^2, na.rm = T) )
        
      } # End loop through grid of emission rates
      
      # If RMSE is monotonically increasing, there is no optimal q 
      if (all(diff(sse) > 0)){
        q.vals[o] <- NA
        
        # Otherwise, save best rate and percent differences
      } else {
        q.vals[o] <- q.grid[which.min(sse)]
      }
      
    } # End loop through samples
    
    # Save mean of sampled emission rates and the 5th and 95th percentiles
    # Note: 3600 multiplier converts from g/s to g/hr
    rate.est.all.events[t] <- mean(q.vals, na.rm = T) * 3600
    error.lower.all.events[t] <- quantile(q.vals, probs = 0.05, na.rm = T) * 3600
    error.upper.all.events[t] <- quantile(q.vals, probs = 0.95, na.rm = T) * 3600
    
    # If there are not enough time steps in which both observations and
    # predictions are in a spike, then do not estimate a rate
  } else {
    rate.est.all.events[t]    <- NA
    error.lower.all.events[t] <- NA
    error.upper.all.events[t] <- NA
  }
  
} # End loop through events

# Check that predictions and observations are on the same order of magnitude
med.p <- median(all.preds.to.compare.all.events, na.rm = T)
med.o <- median(all.obs.to.compare.all.events,   na.rm = T)

if ( med.p > med.o * 10 & med.p < med.o / 10 ){
  if (med.p > med.o){
    print("WARNING: Simulation output and observations are on different order, recommend re-simulating with smaller rate")
  } else {
    print("WARNING: Simulation output and observations are on different order, recommend re-simulating with larger rate")
  }
}

# Package up event detection, localization, and quantification results
to.save <- list(event.mask = spikes,
                max.obs = max.obs,
                localization.estimates = loc.est.all.events,
                rate.estimates = rate.est.all.events,
                error.lower = error.lower.all.events,
                error.upper = error.upper.all.events,
                source.names = source.names,
                WD = data$WD,
                WS = data$WS)

# Save results
saveRDS(to.save, output.file.path)

