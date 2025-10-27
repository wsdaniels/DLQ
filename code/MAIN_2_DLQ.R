# Description: Performs emission event detection, localization, and quantification
#              using output from the Gaussian puff dispersion model and CMS concentration
#              observations.
# Author: William Daniels (wdaniels@mines.edu)
# Last Updated: September 5, 2025

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
values     <- sapply(config, function(X) X[2])

# Get parameter values
gap.time           <- as.numeric(values[parameters == "gap.time"])
length.threshold   <- as.numeric(values[parameters == "length.threshold"])
do.event.detection <- as.logical(values[parameters == "do.event.detection"])
first.sim.ind      <- as.numeric(values[parameters == "first.sim.ind"])
n.samples          <- as.numeric(values[parameters == "n.samples"])

# Get directories
forward.model.path              <- as.character(values[parameters == "forward.model.path"])
output.file.path                <- as.character(values[parameters == "output.file.path"])
helper.spike.detection.alg.path <- as.character(values[parameters == "helper.spike.detection.alg.path"])

# Source helper files which contain helper functions 
source(helper.spike.detection.alg.path)

# Read in simulation data
data <- readRDS(forward.model.path)

# Trim data so that they start and end at either the hour or half hour
first.clean.time <- min(which(minute(data$times) %in% c(0,30)))
last.clean.time <- max(which(minute(data$times) %in% c(0,30)))
to.keep <- first.clean.time:last.clean.time

data$times <- data$times[to.keep]
data$WD <- data$WD[to.keep]
data$WS <- data$WS[to.keep]
data$obs <- data$obs[to.keep, ]
data[first.sim.ind:length(data)] <- lapply(data[first.sim.ind:length(data)], function(X) X[to.keep, ])

# Pull out sensor observations and replace NA's that are not on edge of 
# the time series with interpolated values
obs <- na.approx(data$obs, na.rm = F)

# Number of sensors
n.r <- ncol(obs)

# Pull out time stamps of observations and simulations
times <- data$times

# Pull out the simulation predictions
sims <- data[first.sim.ind:length(data)]

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
  # Print status
  print(paste0("REMOVE BACKGROUND FROM OBSERVATIONS: Sensor ", names(to.use[j]), 
               " - ", j, "/", length(to.use)))
  
  # Grab observations from this sensor
  this.raw.obs <- obs[,j]
  
  # Remove the NA's at the beginning and end of the time series
  # (some sensors start late or end early)
  to.keep       <- !is.na(this.raw.obs)
  trimmed.obs   <- this.raw.obs[to.keep]
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
    } # End loop through spikes
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
      first.ob  <- min(which(spikes$events == event.nums[i]))
      last.ob   <- max(which(spikes$events == event.nums[i]))
      this.mask <- first.ob:last.ob
      spikes$events[this.mask] <- event.nums[i]
      
      # Estimate background using observations before and after spike
      # This is just the first and last observation within the spike mask,
      # since we already added the points before and after the spike to the 
      # spike mask earlier on
      b.left  <- trimmed.obs[first.ob]
      b.right <- trimmed.obs[last.ob]
      b <- mean(c(b.left, b.right))
      
      # Remove background from this spike
      trimmed.obs[this.mask] <- trimmed.obs[this.mask] - b
    } # End loop through events
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
} # End loop through sensor units





# STEP 3: EVENT DETECTION
#---------------------------------------------------------------------------

# Compute maximum concentration value across sensors for each minute
# Note that since all non-spike observations have been set to zero, any 
# non-zero value in the max.obs time series should be considered an event 
# (and hence no need to run the spike detection algorithm again)
max.obs <- apply(obs, 1, max, na.rm = T)

# Estimate emission start and end times 
if (do.event.detection){
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
  } # End loop through times
  
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
  } # End loop through spikes
  
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
  
  # Perform localization and quantification on non-overlapping 30-minute intervals
} else {
  # Get range of times, rounded to 30-minute intervals
  date.range <- range(data$times)
  date.range[1] <- ceiling_date(date.range[1], unit = paste0(gap.time," minutes")) 
  date.range[2] <- floor_date(date.range[2], unit = paste0(gap.time," minutes")) 
  
  # Times separating the 30-minute intervals
  int.breaks <- seq(date.range[1], date.range[2], by = paste0(gap.time," min"))
  
  # Number of intervals
  n.ints <- length(int.breaks) - 1
}





# STEP 4: COMPUTE ALIGNMENT METRIC
#---------------------------------------------------------------------------

# Matrix to hold alignment metric for each event and potential source
metrics <- matrix(NA, nrow = n.ints, ncol = n.s)

# Loop through events
for (t in 1:n.ints){
  # Print status
  print(paste0("COMPUTE ALIGNMENT METRIC: ", t, "/", n.ints))
  
  # If doing event-detection mode, use the interval breaks defined earlier
  if (do.event.detection){
    sub.mask <- times %within% interval(times[min(which(spikes$events == event.nums[t]))], 
                                        times[max(which(spikes$events == event.nums[t]))])
    
    # If doing fixed-interval mode, use fixed interval breaks separated by gap.time
    # Ensure no minutes are double counted between intervals
  } else {
    sub.mask <- times %within% interval(int.breaks[t], int.breaks[t+1] - minutes(1)) 
  }
  
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
      obs.int   <- obs[sub.mask, r]
      preds.int <- preds[sub.mask, r] 
      
      # Only compare non-NA observations
      to.compare <- !is.na(obs.int) 
      
      # Save predictions and observations to compare
      all.obs.to.compare   <- c(all.obs.to.compare,   obs.int[to.compare])
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
        metrics[t, s] <- NA
      } else {
        metrics[t, s] <- cor(x, y, use = "complete.obs", method = "pearson")  
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

# Initialize list of lists to hold quantification results for all events/intervals
big.out <- vector("list", n.ints)

# Loop through events
for (t in 1:n.ints){
  # Print status
  print(paste0("COMPUTE LOCALIZATION AND QUANTIFICATION ESTIMATES: ", t, "/", n.ints))
  
  # Initialize vectors to hold quantification results
  q.hat <- q.hat.lower <- q.hat.upper <- vector(length = n.s)
  
  # Get metric values for this event
  these.metrics <- metrics[t, ]
  
  # If doing event-detection mode, use the interval breaks defined earlier
  if (do.event.detection){
    sub.mask <- times %within% interval(times[min(which(spikes$events == event.nums[t]))], 
                                        times[max(which(spikes$events == event.nums[t]))])
    
    # If doing fixed-interval mode, use fixed interval breaks separated by gap.time
    # Ensure no minutes are double counted between intervals
  } else {
    sub.mask <- times %within% interval(int.breaks[t], int.breaks[t+1] - minutes(1)) 
  }
  
  # Create matrix that holds dispersion model output
  X <- matrix(NA, nrow = sum(sub.mask) * n.r, ncol = n.s)
  colnames(X) <- source.names
  
  # Fill X matrix with simulation output
  # Loop through sources
  for (s in 1:n.s){
    this.source <- c()
    # Loop through sensors
    for (r in 1:n.r){
      this.source <- c(this.source, sims[[s]][sub.mask, r])
    } # End loop through sensors
    X[, s] <- this.source
  } # End loop through sources
  
  # Determine which sources have downwind sensors (i.e., which sources have information)
  # Here, we define a source as having 'information' if there are
  # more than 4 minutes of simulated observations greater than 0.25 ppm
  info.mask <- apply(X, 2, function(this.col) sum(this.col > 0.25) > 4)
  q.hat[!info.mask] <- q.hat.lower[!info.mask] <- q.hat.upper[!info.mask] <- "no info"
  
  # If all sources have no information, move onto the next event/interval
  if(all(!info.mask)) {
    loc.est.all.events[t] <- rate.est.all.events[t] <- 
      error.lower.all.events[t] <- error.upper.all.events[t] <- "no info"
    
    next
  }
  
  # Get source name corresponding to largest metric value
  loc.est.all.events[t] <- source.names[which.max(these.metrics[info.mask])]
  
  # Get predictions from most likely source and observations during this event
  event.preds <- sims [[which.max(these.metrics[info.mask])]][sub.mask, ]
  event.obs   <- obs  [sub.mask, ]
  event.times <- times[sub.mask]
  
  # Vectors to hold preds and obs when both are in a spike
  all.preds.to.compare <- all.obs.to.compare <- c()
  
  # Loop through sensors and save when both obs and preds are in a spike
  for (r in 1:n.r){
    # Get predictions and observations for this sensor
    these.preds <- event.preds[, r]
    these.obs   <- event.obs[, r]
    
    # If all observations from this sensor is NA, move onto the next sensor
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
  } # End loop through sensors
  
  # If there is enough data to compare, compute rate estimate and 
  # percent difference between obs and preds (will be used for UQ)
  if (length(all.preds.to.compare) > 4){
    
    # Save obs and preds so that we can check that they are on the same order of magnitude later
    all.preds.to.compare.all.events <- c(all.preds.to.compare.all.events, all.preds.to.compare)
    all.obs.to.compare.all.events   <- c(all.obs.to.compare.all.events,   all.obs.to.compare)
    
    # Vector to hold emission rate estimate for each sample
    q.vals <- vector(length = n.samples)
    
    # Loop through samples
    for (o in 1:n.samples){
      # Sample from predictions and observations
      this.sample <- sample.int(n = length(all.preds.to.compare), 
                                size = length(all.preds.to.compare)/2,
                                replace = T)
      
      # --------- OLD GRID SEARCH CODE -------------------------
      
      # # Define a grid of emission rate values to optimize over
      # q.grid <- lseq(0.0001, 3000, length.out = 2000)
      # 
      # # Vector to hold sum of squared errors
      # sse <- vector(length = length(q.grid))
      # 
      # # Loop through possible emission rates
      # for (z in 1:length(q.grid)){
      # 
      #   # Scale predictions by this emission rate
      #   qxp <- q.grid[z] * all.preds.to.compare
      # 
      #   # Compute root mean square error
      #   sse[z] <- sqrt( mean( (all.obs.to.compare[this.sample] - qxp[this.sample])^2, na.rm = T) )
      # 
      # } # End loop through grid of emission rates
      # 
      # # If RMSE is monotonically increasing, there is no optimal q
      # if (all(diff(sse) > 0)){
      #   q.vals[o] <- NA
      #   # Otherwise, save best rate and percent differences
      # } else {
      #   q.vals[o] <- q.grid[which.min(sse)]
      # }
      
      # --------- END OLD GRID SEARCH CODE ----------------------
      
      # Binary search code written by Troy Sorensen. 
      # Implemented here on June 11, 2025. 
      
      # ------------ NEW BINARY SEARCH CODE ---------------------
      preds_s <- all.preds.to.compare[this.sample]
      obs_s   <- all.obs.to.compare[this.sample]
      
      # Set bounds for search (grid search had 0 to 3000 on a logarithmic grid)
      q_low  <- 0
      q_high <- 3000
      tol    <- 1e-3
      
      # Check derivative of upper and lower bounds of the search range
      deriv_low  <- sum(preds_s * (obs_s - q_low    * preds_s))
      deriv_high <- sum(preds_s * (obs_s - q_high * preds_s))
      
      if (all(preds_s == 0)) {
        # If all predictions are zero, no slope: cannot estimate q, set to NA
        q.vals[o] <- NA
      } else if (deriv_low * deriv_high > 0) {
        # Derivatives have same sign: thus no root lies in [q_low,q_high].
        # Set q to q_low or q_high accordingly.
        q.vals[o] <- if (deriv_low > 0) q_high else q_low
      } else {
        # Root lies in [q_low,q_high]. Find via binary search.
        # Stop when abs(deriv) < tol or 50 iterations are performed, whichever occurs first.
        for (iter in 1:50) {
          q_mid <- (q_low + q_high) / 2
          
          # derivative of SSE(q) w.r.t. q is sum(preds * (obs - q*preds))
          deriv <- sum(preds_s * (obs_s - q_mid * preds_s), na.rm = TRUE)
          
          if (is.na(deriv)) {
            q_mid <- NA
            break
          }
          if (abs(deriv) < tol) {
            break
          }
          # if deriv > 0, SSE is still decreasing as q increases → move lower bound up
          if (deriv > 0) {
            q_low <- q_mid
          } else {
            q_high <- q_mid
          }
        }
        
        q.vals[o] <- q_mid
      }
      # ------------ END BINARY SEARCH CODE ---------------------
      
    } # End loop through samples
    
    # Note: 3.6 multiplier converts from g/s to kg/hr
    q.vals <- q.vals * 3.6
    
    # Save mean of sampled emission rates and the 5th and 95th percentiles
    q.hat[which(info.mask)[which.max(these.metrics[info.mask])]] <- 
      rate.est.all.events[t] <- mean(q.vals, na.rm = TRUE)
    q.hat.lower[which(info.mask)[which.max(these.metrics[info.mask])]] <- 
      error.lower.all.events[t] <- quantile(q.vals, probs = 0.05, na.rm = T)
    q.hat.upper[which(info.mask)[which.max(these.metrics[info.mask])]] <- 
      error.upper.all.events[t] <- quantile(q.vals, probs = 0.95, na.rm = T)
    
    # Leaning into the single source assumption, set the 
    # emission rate for all other potential sources to 0
    q.hat[q.hat == "FALSE"] <- 0
    q.hat.lower[q.hat.lower == "FALSE"] <- 0
    q.hat.upper[q.hat.upper == "FALSE"] <- 0
    
    # If there are not enough time steps in which both observations and
    # predictions are in a spike, then do not estimate a rate.
    # NOTE: differentiate between events/intervals where we do not think emissions are happening,
    # i.e., when observations are all zero, and times when we think emissions are 
    # happening, but when we don't have enough alignment between simulation and 
    # observations to estimate an emission rate.
  } else {
    # No concentration enhancements, likely that no emissions are happening, so set rate to 0
    if (all(max.obs[sub.mask] == 0)) {
      rate.est.all.events[t]    <- 0
      error.lower.all.events[t] <- 0
      error.upper.all.events[t] <- 0
      
      q.hat[q.hat == "FALSE"] <- 0
      q.hat.lower[q.hat.lower == "FALSE"] <- 0
      q.hat.upper[q.hat.upper == "FALSE"] <- 0
      q.vals <- rep(0, n.samples)
      # Emissions could be happening, but not enough overlap to quantify
    } else {
      rate.est.all.events[t]    <- "no overlap"
      error.lower.all.events[t] <- "no overlap"
      error.upper.all.events[t] <- "no overlap"
      
      q.hat[q.hat == "FALSE"] <- "no overlap"
      q.hat.lower[q.hat.lower == "FALSE"] <- "no overlap"
      q.hat.upper[q.hat.upper == "FALSE"] <- "no overlap"
      q.vals <- "no overlap"
    }
  }
  
  # Save output
  out <- list(q.hat = q.hat,
              q.hat.lower = q.hat.lower, 
              q.hat.upper = q.hat.upper,
              rates = q.vals)
  big.out[[t]] <- out
  
} # End loop through events

# Check that predictions and observations are on the same order of magnitude
med.p <- median(all.preds.to.compare.all.events, na.rm = T)
med.o <- median(all.obs.to.compare.all.events,   na.rm = T)

if (med.p > med.o * 10 | med.p < med.o / 10){
  if (med.p > med.o){
    print("WARNING: Simulation output and observations are on different order, recommend re-simulating with smaller rate")
  } else {
    print("WARNING: Simulation output and observations are on different order, recommend re-simulating with larger rate")
  }
}

# Package up DLQ results
to.save <- c(list(event.mask = spikes,
                  max.obs = max.obs,
                  localization.estimates = loc.est.all.events,
                  rate.estimates = rate.est.all.events,
                  error.lower = error.lower.all.events,
                  error.upper = error.upper.all.events,
                  source.names = source.names,
                  WD = data$WD,
                  WS = data$WS, 
                  out = big.out),
             sims)

# Save results
saveRDS(to.save, output.file.path)
