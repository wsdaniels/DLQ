# Description: 
# Summarizes event detection, localization, and quantification results and 
# reproduces paper figures.
# Author: William Daniels (wdaniels@mines.edu)
# Last Updated: December 2023

# Clear environment
if(!is.null(dev.list())){dev.off()}
rm(list = ls())

# Import necessary libraries
library(lubridate)
library(scales)
# library(circular)
library(rstudioapi)

if (commandArgs()[1] == "RStudio"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}


# START USER INPUT
#---------------------------------------------------------------------------

# Directory to save plots
save.dir <- '../'

# Path to event detection, localization, quantification results
data.path <- '../output_data/DLQ_output.RData'

# Path to leak (truth) data
leak.path <- '../input_data/leak_data.csv'

# Set colors for plots
tank.color <- "#3062CF" #blue
wellhead.east.color <- "#9147B8" #purple
separator.east.color <- "#7EAD52" #green
separator.west.color <- "#F1C30E" #gold
wellhead.west.color <- "#C7383C" #red

# END OF USER INPUT - NO MODIFICATION NECESSARY BELOW THIS POINT
#---------------------------------------------------------------------------






# STEP 1: PARSE OUT RESULTS
#---------------------------------------------------------------------------

# Read in event detection, localization, and quantification results
data <- readRDS(data.path)

# Pull out times associated with event detection results
times <- data$event.mask$time

# Pull out event mask
events <- data$event.mask

# Pull out integers that uniquely identify events
event.nums <- na.omit(unique(events$events))

# Pull out source names and number of sources
source.names <- data$source.names
n.s <- length(source.names)

# Pull out time series of maximum concentration across sensors for each minute
max.obs <- data$max.obs



# STEP 2: READ IN TRUE LEAK DATA AND ELIMINATE MULTISOURCE LEAKS
#---------------------------------------------------------------------------

# Read in leak data
leak.data <- suppressWarnings(read.csv(leak.path))

# Determine unique equipment group IDs
unique.eg <- unique(leak.data$tc_EquipmentGroupID)

# Names to replace equipment group IDs. Order is set manually.
eg.names <- c("Tank", "Wellhead.East", "Wellhead.West", "Separator.East", "Separator.West")

# Replace equipment group ID values with names
for (i in 1:length(unique.eg)){
  this.mask <- leak.data$tc_EquipmentGroupID == unique.eg[i]
  leak.data$tc_EquipmentGroupID[this.mask] <- eg.names[i]
}

# Convert time strings to datetime objects
leak.data$tc_ExpStartDatetime <- as_datetime(as_datetime(leak.data$tc_ExpStartDatetime), tz = "America/Denver")
leak.data$tc_ExpEndDatetime <- as_datetime(as_datetime(leak.data$tc_ExpEndDatetime), tz = "America/Denver")

# Mask in leak data during simulation time frame
to.keep <- leak.data$tc_ExpStartDatetime >= times[1] & leak.data$tc_ExpStartDatetime <= times[length(times)]
leak.data <- leak.data[to.keep,]

# Sort by leak start time
leak.data <- leak.data[order(leak.data$tc_ExpStartDatetime), ]

# Figure out when multisource leaks are happening
multisource <- rep(F, nrow(leak.data))
for (i in 1:(nrow(leak.data)-1)){
  for (j in (i+1):nrow(leak.data)){
    
    first.leak.start <- leak.data$tc_ExpStartDatetime[i]
    first.leak.stop  <- leak.data$tc_ExpEndDatetime[i]
    first.leak.times <- seq(first.leak.start, first.leak.stop, by = "1 min")
    first.leak.times <- round_date(first.leak.times, unit = "min")
    
    second.leak.start <- leak.data$tc_ExpStartDatetime[j]
    second.leak.stop  <- leak.data$tc_ExpEndDatetime[j]
    second.leak.times <- seq(second.leak.start, second.leak.stop, by = "1 min")
    second.leak.times <- round_date(second.leak.times, unit = "min")
    
    if (any(first.leak.times %in% second.leak.times)){
      multisource[c(i,j)] <- T
    }    
    
  }
}


# Separate out the single and multisource leaks
singlesource.leaks <- leak.data[!multisource, ]
multisource.leaks <- leak.data[multisource, ]

# Number of METEC emission events
n.leaks <- nrow(singlesource.leaks)

# Plot time series of true emissions
plot.divs <- c(seq(min(times), max(times)+days(3), by = "5 days"))

axis.points <- seq(min(round_date(times, "day")-days(1)),
                   max(round_date(times, "day")+days(1)),
                   by = "day")
axis.labs <- paste0(month.abb[month(axis.points)],
                    "-", day(axis.points))

minor.axis.points <- seq(min(round_date(times, "day")-days(1)),
                         max(round_date(times, "day")+days(1)),
                         by = "12 hours")

for (i in 1:(length(plot.divs)-1)){
  
  png(paste0(save.dir, "leaks_", i, ".png"),
      res = 100, pointsize = 24, width = 1920, height = 1080)
  
  plot(leak.data$tc_ExpStartDatetime, leak.data$tc_C1MassFlow, pch = ".",
       xaxt= "n",
       xlim = c(plot.divs[i], plot.divs[i+1]))
  
  segments(x0 = leak.data$tc_ExpStartDatetime, y0 = leak.data$tc_C1MassFlow,
           x1 = leak.data$tc_ExpEndDatetime,   y1 = leak.data$tc_C1MassFlow,
           lwd = 2)
  
  axis(side = 1, at = axis.points, labels = axis.labs, lwd = 4)
  axis(side = 1, at = minor.axis.points, labels = NA)
  
  for (i in 1:length(multisource)){
    if (multisource[i]){
      abline(v = leak.data$tc_ExpStartDatetime[i], col = "red")
    }
  }
  
  dev.off()
  
}

# Initialize vector to hold estimated events that should be removed because
# they overlap with a multisource leak
to.remove <- vector(length = length(event.nums))

# Remove events that overlap with a multisource leak
for (i in 1:length(event.nums)){
  
  # Mask in this event
  this.spike <- which(events$events == event.nums[i])
  these.times <- events$time[this.spike]
  
  # Loop through multisource leaks
  for (j in 1:nrow(multisource.leaks)){
    
    leak.start <- multisource.leaks$tc_ExpStartDatetime[j]
    leak.stop  <- multisource.leaks$tc_ExpEndDatetime[j]
    
    # If there is any overlap between estimated event and leak, remove
    # the event
    if (any(these.times >= leak.start & these.times <= leak.stop)){
      events$events[this.spike] <- NA
      to.remove[i] <- T
      break
    }
  }
}


# Get event numbers after filtering out events that overlap with multisource leaks
event.nums <- na.omit(unique(events$events))

# Number of events after filtering multisource leaks
n.ints <- length(event.nums)

# Remove location, rate, and percent differences for events that overlapped with
# multisource leak
loc.est.all.events <- data$localization.estimates[!to.remove]
rate.est.all.events <- data$rate.estimates[!to.remove]
error.lower.all.events <- data$error.lower[!to.remove]
error.upper.all.events <- data$error.upper[!to.remove]


# Plot events
col.vals <- rep(c("blue", "orange", "forestgreen"), length(event.nums))

plot.divs <- c(seq(min(times), max(times)+days(3), by = "5 days"))

axis.points <- seq(min(round_date(times, "day")-days(1)),
                   max(round_date(times, "day")+days(1)),
                   by = "day")
axis.labs <- paste0(month.abb[month(axis.points)],
                    "-", day(axis.points))

minor.axis.points <- seq(min(round_date(times, "day")-days(1)),
                         max(round_date(times, "day")+days(1)),
                         by = "12 hours")

for (i in 1:(length(plot.divs)-1)){
  
  png(paste0(save.dir, "events_", i, ".png"),
      res = 100, pointsize = 24, width = 1920, height = 1080)
  
  plot(times, max.obs, type = "l", xaxt= "n",
       ylab = "Methane Concentration [ppb]",
       xlim = c(plot.divs[i], plot.divs[i+1]))
  
  for (i in 1:length(event.nums)){
    this.spike <- which(events$events == event.nums[i])
    lines(events$time[this.spike], max.obs[this.spike], col = col.vals[i])
  }
  
  axis(side = 1, at = axis.points, labels = axis.labs, lwd = 4)
  axis(side = 1, at = minor.axis.points, labels = NA)
  
  dev.off()
  
}



# STEP 3: PREP DETECTION RESULTS
#---------------------------------------------------------------------------

# Initialize matrix to hold event detection results.
# Rows are the METEC emission events, columns are our predicted events.
# True indicates an estimated event overlaps with a true leak.
overlap.matrix <- matrix(NA, nrow = n.leaks, ncol = n.ints)
overlap.matrix[is.na(overlap.matrix)] <- F

# Fill in the overlap matrix. First loop over METEC emissions.
for (l in 1:n.leaks){
  
  # Grab start and stop times of the METEC emission
  leak.start <- singlesource.leaks$tc_ExpStartDatetime[l]
  leak.end   <- singlesource.leaks$tc_ExpEndDatetime[l]
  
  # Create minute sequence over duration of METEC emission
  leak.times <- seq(leak.start, leak.end, by = "1 min")
  leak.times <- round_date(leak.times, unit = "min")
  
  # Loop over predicted events.
  for (t in 1:n.ints){
    
    # Create minute sequence over duration of this predicted event.
    this.mask <- seq(min(which(events$events == event.nums[t])),
                     max(which(events$events == event.nums[t])))
    event.times <- events$time[this.mask]
    
    # Store overlap status
    if (any(leak.times %in% event.times)){
      overlap.matrix[l,t] <- T
    }
  }
}


# Create mask for METEC emission that were identified
identified.leaks <- apply(overlap.matrix, 1, function(X) any(X))
n.identified.leaks <- sum(identified.leaks)

# Create mask for false positive predictions (predicted event when there was no METEC emission)
false.positives <- apply(overlap.matrix, 2, function(X) all(!X))
n.false.positives <- sum(false.positives)

# Event-level false positive rate
round(100 * n.false.positives / n.ints, 1)

# RESULT 1: START AND STOP TIME ERROR
#---------------------------------------------------------------------------

# Initialize vectors to hold errors and coverage for each identified leak
start.time.errors <- end.time.errors <- vector(length = n.identified.leaks)

# Initialize vector that will be TRUE for each minute of ALL METEC emissions
# where we predict an emission and FALSE for each minute of ALL METEC emissions
# where we do not predict and emission
count <- 1

# Loop through identified leaks only, as there is no way to compute timing error
# when we do not predict an event
for (l in which(identified.leaks)){
  
  # Get start and end times of METEC emission
  leak.start <- floor_date(singlesource.leaks$tc_ExpStartDatetime[l], unit = "min")
  leak.end   <- singlesource.leaks$tc_ExpEndDatetime[l]
  
  # Get start and end times of predicted emissions.
  # Note that if two predicted emissions fall within one METEC emission, the
  # start time of the first predicted emission and the end time of the second
  # predicted emission is taken.
  event.ind <- which(overlap.matrix[l, ])
  event.start <- events$time[which(events$events == event.nums[min(event.ind)])][1]
  event.end <- events$time[which(events$events == event.nums[max(event.ind)])]
  event.end <- event.end[length(event.end)]
  
  # Save error between start and stop times
  start.time.errors[count] <- difftime(event.start, leak.start, units = "mins")
  end.time.errors[count] <- difftime(event.end, leak.end, units = "mins")
  
  count <- count + 1
  
}


under.est.col <- "#00524F"
over.est.col <- "#29A385"
accent.col <- "#FB514B"


png(paste0(save.dir, "timing_error.png"),
    width = 1920, height = 600, res = 100, pointsize = 24)

par(mfrow = c(1,2))
par(mgp = c(2.25, 1, 0))
par(mar = c(3.5, 2, 1, 1))
par(oma = c(0,1.5,0,0))

under.est.mask <- start.time.errors < 0

hist(start.time.errors, breaks = 40, col = over.est.col,
     xlim = c(-200, 200),
     ylim = c(0,55),
     xlab = "Start Time Error [min]",
     right = F, main = "")

hist(start.time.errors[under.est.mask], breaks = 4, col = under.est.col, add = T,
     right = F)

abline(v = median(start.time.errors), col = accent.col, lwd = 4)

mtext("Number of Emission Events", side = 2, outer = T,
      line = 0.5)

under.est.mask <- end.time.errors > 0

hist(-end.time.errors, breaks = 40, col = over.est.col,
     xlim = c(-200, 200),
     ylim = c(0,55),
     xlab = "End Time Error [min]",
     right = F, main = "")

hist(-end.time.errors[under.est.mask], breaks = 10, col = under.est.col, add = T,
     right = F)

abline(v = median(-end.time.errors), col = accent.col, lwd = 4)

dev.off()


round(quantile(start.time.errors, probs = c(0, 0.25, 0.5, 0.75, 1)), 0)
round(quantile(-1*end.time.errors,   probs = c(0, 0.25, 0.5, 0.75, 1)), 0)


# RESULT 2: CONFUSION MATRIX
#---------------------------------------------------------------------------

# Initialize vector to hold times when a true leak is occurring
leak.times <- c()

# Loop through leaks
for (i in 1:nrow(singlesource.leaks)){
  
  # Get start and stop times of this leak
  leak.start <- singlesource.leaks$tc_ExpStartDatetime[i]
  leak.end   <- singlesource.leaks$tc_ExpEndDatetime[i]
  
  # Create minute sequence spanning this leak
  these.leak.times <- seq(leak.start, leak.end, by = "1 min")
  these.leak.times <- round_date(these.leak.times, unit = "min")
  
  # Save leak times
  leak.times <- c(leak.times, these.leak.times)
}

# Create minute-by-minute leak mask (TRUE = leak is occurring)
true.emission.mask <- rep(F, length(times))
true.emission.mask[times %in% leak.times] <- T

# Create minute-by-minute mask for predicted emission events
# (TRUE = emission event is predicted at that time step)
predicted.emission.mask <- !is.na(events$events)

# There is a leak and we predict a leak
tp <- sum(true.emission.mask & predicted.emission.mask)

# There is no leak and we predict a leak
fp <- sum(!true.emission.mask & predicted.emission.mask)

# There is a leak and we do not predict a leak
fn <- sum(true.emission.mask & !predicted.emission.mask)

# There is no leak and we do not predict a leak
tn <- sum(!true.emission.mask & !predicted.emission.mask)

# Compute summary statistics
true.positive.rate <- round(100 * tp / (tp+fn), 1)
true.negative.rate <- round(100 * tn / (fp+tn), 1)
accuracy <- round(100 * (tp+tn) / (tp+fp+fn+tn), 1)
npv <-      round(100 * tn / (fn+tn), 1)
ppv <-      round(100 * tp / (tp+fp), 1)
total <- tp + fn + fp + tn

# Create confusion matrix
rbind(c(tp, round(100*tp/total,1), fp, round(100*fp/total,1), ppv),
      c(fn, round(100*fn/total,1), tn, round(100*tn/total,1), npv),
      c(true.positive.rate, NA, true.negative.rate, NA, accuracy))



# RESULT 3: EVENT DETECTION, LOCALIZATION, AND AVAILABILITY OF RATE ESTIMATES
#---------------------------------------------------------------------------

# Initialize vectors to hold true source location and availability of rate estimate
# regardless of true source location.
all.sources.true.locs <- all.sources.rate.est.available <- vector(length = n.leaks)

# Initialize list to hold location estimates regardless of true source location.
all.sources.loc.est <- vector(mode = "list", length = n.leaks)

# Loop through identified leaks only, since we cannot create a localization estimate
# or rate estimate for events that we do not detect
for (l in 1:n.leaks){
  
  # Get index of all predicted events that overlap with this METEC emission 
  event.ind <- which(overlap.matrix[l, ])
  
  # Save all localization estimates for this METEC emission
  all.sources.loc.est[[l]] <- loc.est.all.events[event.ind]
  
  # Determine if there is a rate estimate available during this leak
  all.sources.rate.est.available[l] <- any(!is.na(rate.est.all.events[event.ind]))
  
  # Save true emission location for this METEC emission
  all.sources.true.locs[l] <- singlesource.leaks$tc_EquipmentGroupID[l]
  
}

# Initialize lists to hold summary statistics separated by true source location
source.separated.completely.correct <- source.separated.partially.correct <-
  source.separated.identified <- source.separated.rate.est.available <- 
  vector(mode = "list", length = n.s)

# Loop through sources
for (s in 1:n.s){
  
  # Mask in METEC emissions that occurred at this source
  this.mask <- all.sources.true.locs == source.names[s]
  
  # Mask in location estimates and rate availability for METEC emissions
  # occurring at this source
  these.loc.ests <- all.sources.loc.est[this.mask]
  these.rate.availabilities <- all.sources.rate.est.available[this.mask]
  
  # Initialize each component of the summary lists to hold results for each emission
  source.separated.completely.correct[[s]] <- source.separated.partially.correct[[s]] <- 
    source.separated.identified[[s]] <- source.separated.rate.est.available[[s]] <- 
    vector(length = sum(this.mask))
  
  # Loop through emissions occurring at this source
  for (p in 1:sum(this.mask)){
    
    # If all predicted locations overlapping with this true emission are correct, record a TRUE
    if (all(these.loc.ests[[p]] == source.names[s]) & length(these.loc.ests[[p]]) > 0){
      source.separated.completely.correct[[s]][p] <- T
    } 
    
    # If any predicted locations overlapping with this true emission are correct, record a TRUE
    if (any(these.loc.ests[[p]] == source.names[s]) & length(these.loc.ests[[p]]) > 0){
      source.separated.partially.correct[[s]][p] <- T
    }
    
    # If there is a rate estimate available during this true emission, record a TRUE
    if (these.rate.availabilities[p]){
      source.separated.rate.est.available[[s]][p] <- T
    }
    
    # If there is a location estimate overlapping with this true event, it means we detected
    # the event, so record a TRUE
    if (length(these.loc.ests[[p]]) > 0){
      source.separated.identified[[s]][p] <- T
    }
  }
}

# Put everything together
counts <- cbind(sapply(source.separated.completely.correct, function(X) sum(X)),
                sapply(source.separated.partially.correct,  function(X) sum(X)),
                sapply(source.separated.rate.est.available, function(X) sum(X)),
                sapply(source.separated.identified,         function(X) sum(X))) / n.leaks


png(paste0(save.dir, "detection_localization_results.png"),
    width = 1920, height = 750, res = 100, pointsize = 24)

par(mar = c(4,10,1,5))

barplot(counts,
        xlab = "Percent of METEC Emission Events",
        horiz = T, xlim = c(0,1),
        col = alpha(c(tank.color, separator.east.color, wellhead.east.color,
                      wellhead.west.color, separator.west.color), 0.85),
        las = 2, 
        xaxt = "n",
        names.arg = c("Location Estimate\nCompletely Correct",
                      "Location Estimate\nPartially Correct",
                      "Rate Estimate\nAvailable",
                      "Correctly Identified\nOccurrence of Emission"))

axis(side = 1, at = seq(0,1, by = 0.1), labels = seq(0,100, by = 10))

dev.off()




# RESULT 4: QUANTIFICATION
#---------------------------------------------------------------------------

# Initialize vectors to hold rate estimates, true rates, ratios of estimated to true rates,
# percent differences between observations and predictions, and uncertainty on true rates
rate.est.identified.events <- rate.est.ratios <- true.rates.identified.events <- 
  true.rate.uncertainty.identified.events <-  rate.est.errors <- 
  error.lower.identified.events <- error.upper.identified.events <- 
  rate.est.pdiffs <- vector(length = n.identified.leaks)
count <- 1

# Loop through identified leaks only
for (l in which(identified.leaks)){
  
  # Get true METEC emission rate (plus uncertainty) and average of predicted rates during this METEC emission
  event.ind <- which(overlap.matrix[l, ])
  true.rate <- singlesource.leaks$tc_C1MassFlow[l]
  rate.est <- mean(rate.est.all.events[event.ind], na.rm = T)
  true.rate.uncertainty <- singlesource.leaks$tc_C1MassFlowUncertainty[l]
  
  # Save true rate, estimated rate, ratio of estimate to truth, and truth uncertainty
  true.rates.identified.events[count] <- true.rate
  rate.est.identified.events[count] <- rate.est
  rate.est.ratios[count] <- rate.est / true.rate
  rate.est.errors[count] <- rate.est - true.rate
  rate.est.pdiffs[count] <- (rate.est - true.rate)/true.rate
  true.rate.uncertainty.identified.events[count] <- true.rate.uncertainty
  
  # Save percent differences
  error.lower.identified.events[count] <- mean(error.lower.all.events[event.ind], na.rm = T)
  error.upper.identified.events[count] <- mean(error.upper.all.events[event.ind], na.rm = T)
  
  count <- count + 1
  
}



# RESULT 4.1: PLOT TIME SERIES OF RESULTS
#---------------------------------------------------------------------------

# Get start, mid, and end times for each event
event.start.times <- event.end.times <- event.mid.times <- vector(length = length(event.nums))
for (i in 1:length(event.nums)){
  event.start.times[i] <- events$time[min(which(events$events == event.nums[i]))]
  event.end.times[i]   <- events$time[max(which(events$events == event.nums[i]))]
  event.mid.times[i]   <- mean(c(event.start.times[i], event.end.times[i]))
}

# Set plotting colors for leak location truth data
true.rate.cols <- vector(length = n.leaks)
true.rate.cols[singlesource.leaks$tc_EquipmentGroupID == "Tank"] <- tank.color
true.rate.cols[singlesource.leaks$tc_EquipmentGroupID == "Wellhead.West"] <- wellhead.west.color
true.rate.cols[singlesource.leaks$tc_EquipmentGroupID == "Separator.West"] <- separator.west.color
true.rate.cols[singlesource.leaks$tc_EquipmentGroupID == "Wellhead.East"] <- wellhead.east.color
true.rate.cols[singlesource.leaks$tc_EquipmentGroupID == "Separator.East"] <- separator.east.color

# Set plotting colors for estimated leak location data
rate.est.all.events.cols <- vector(length = n.ints)
rate.est.all.events.cols[loc.est.all.events == "Tank"] <- tank.color
rate.est.all.events.cols[loc.est.all.events == "Wellhead.West"] <- wellhead.west.color
rate.est.all.events.cols[loc.est.all.events == "Separator.West"] <- separator.west.color
rate.est.all.events.cols[loc.est.all.events == "Wellhead.East"] <- wellhead.east.color
rate.est.all.events.cols[loc.est.all.events == "Separator.East"] <- separator.east.color


plot.divs <- c(seq(min(times), max(times)+days(3), by = "5 days"))

axis.points <- seq(min(round_date(times, "day")-days(1)),
                   max(round_date(times, "day")+days(1)),
                   by = "day")
axis.labs <- paste0(month.abb[month(axis.points)],
                    "-", day(axis.points))

minor.axis.points <- seq(min(round_date(times, "day")-days(1)),
                         max(round_date(times, "day")+days(1)),
                         by = "6 hours")


png(paste0(save.dir, "result_time_series.png"), 
    res = 100, width = 1920, height = 2400, pointsize = 32)

par(mgp = c(2.5, 1, 0))
par(mar = c(1.75, 3.75, 1, 1))
par(mfrow = c(6,1))

for (i in 1:(length(plot.divs)-1)){
  
  ylim.max <- 8000
  
  plot(times, max.obs, col = "white", 
       ylim = c(0,ylim.max/1000),
       xlim = c(plot.divs[i], plot.divs[i+1]),
       ylab = "",
       xlab = "",
       xaxt= "n")
  
  axis(side = 1, at = axis.points, labels = rep("", length(axis.points)), lwd = 4)
  axis(side = 1, at = axis.points, labels = axis.labs, lwd = 0, line = -0.4)
  axis(side = 1, at = minor.axis.points, labels = NA)
  
  rect(xleft = event.start.times, 
       xright = event.end.times,
       ybottom = 0, ytop = ylim.max/1000, 
       col = alpha(rate.est.all.events.cols, 0.25))
  
  rect(xleft  = singlesource.leaks$tc_ExpStartDatetime, 
       xright = singlesource.leaks$tc_ExpEndDatetime,
       ybottom = 0, ytop = singlesource.leaks$tc_C1MassFlow/1000,
       col = alpha(true.rate.cols, 1))
  
  points(event.mid.times, rate.est.all.events/1000, pch = 19)
  
  segments(x0 = event.mid.times, 
           y0 = error.lower.all.events/1000,
           y1 = error.upper.all.events/1000,
           lwd = 2)
  
  if (i == length(plot.divs)-1){
    
    # legend("topright",
    #        legend = c("Tanks", "West Wellhead", "West Separator", "East Wellhead", "East Separator",
    #                   "Rate Estimate with Error Bound", NA),
    #        fill = c(tank.color, wellhead.west.color, separator.west.color, wellhead.east.color, separator.east.color,
    #                 NA, NA),
    #        border = c("black", "black", "black", "black", "black", NA, NA),
    #        pch = c(NA, NA, NA, NA, NA, 19, NA),
    #        lty = c(NA, NA, NA, NA, NA, NA, 1),
    #        bty = "n")
    
  }
  
}

mtext("Methane Emission Rate [kg/hr]", side = 2, outer = T, line = -1.25)

dev.off()




# RESULT 4.2: PARITY PLOT
#---------------------------------------------------------------------------

# Set shape for different true source locations
tank.pch <- 21 # circle
wellhead.east.pch <- 22 # square
separator.east.pch <- 23 # diamond
separator.west.pch <- 24 # up triangle
wellhead.west.pch <- 25 # down triangle

# Set color for identified events based on true source location
rate.est.identified.events.cols <- vector(length = n.identified.leaks)
rate.est.identified.events.cols[singlesource.leaks$tc_EquipmentGroupID[identified.leaks] == "Tank"] <- tank.color
rate.est.identified.events.cols[singlesource.leaks$tc_EquipmentGroupID[identified.leaks] == "Wellhead.West"] <- wellhead.west.color
rate.est.identified.events.cols[singlesource.leaks$tc_EquipmentGroupID[identified.leaks] == "Separator.West"] <- separator.west.color
rate.est.identified.events.cols[singlesource.leaks$tc_EquipmentGroupID[identified.leaks] == "Wellhead.East"] <- wellhead.east.color
rate.est.identified.events.cols[singlesource.leaks$tc_EquipmentGroupID[identified.leaks] == "Separator.East"] <- separator.east.color

# Shape shape for identified events based on true source location
rate.est.identified.events.pch <- vector(length = n.identified.leaks)
rate.est.identified.events.pch[singlesource.leaks$tc_EquipmentGroupID[identified.leaks] == "Tank"] <- tank.pch
rate.est.identified.events.pch[singlesource.leaks$tc_EquipmentGroupID[identified.leaks] == "Wellhead.West"] <- wellhead.west.pch
rate.est.identified.events.pch[singlesource.leaks$tc_EquipmentGroupID[identified.leaks] == "Separator.West"] <- separator.west.pch
rate.est.identified.events.pch[singlesource.leaks$tc_EquipmentGroupID[identified.leaks] == "Wellhead.East"] <- wellhead.east.pch
rate.est.identified.events.pch[singlesource.leaks$tc_EquipmentGroupID[identified.leaks] == "Separator.East"] <- separator.east.pch

# Pick out missed leaks
non.identified.leaks <- singlesource.leaks$tc_C1MassFlow[!identified.leaks]
n.non.identified.leaks <- length(non.identified.leaks)

# Set color for non-identified events based on true source location
rate.est.non.identified.events.cols <- vector(length = n.non.identified.leaks)
rate.est.non.identified.events.cols[singlesource.leaks$tc_EquipmentGroupID[!identified.leaks] == "Tank"] <- tank.color
rate.est.non.identified.events.cols[singlesource.leaks$tc_EquipmentGroupID[!identified.leaks] == "Wellhead.West"] <- wellhead.west.color
rate.est.non.identified.events.cols[singlesource.leaks$tc_EquipmentGroupID[!identified.leaks] == "Separator.West"] <- separator.west.color
rate.est.non.identified.events.cols[singlesource.leaks$tc_EquipmentGroupID[!identified.leaks] == "Wellhead.East"] <- wellhead.east.color
rate.est.non.identified.events.cols[singlesource.leaks$tc_EquipmentGroupID[!identified.leaks] == "Separator.East"] <- separator.east.color

# Set shape for non-identified events based on true source location
rate.est.non.identified.events.pch <- vector(length = n.non.identified.leaks)
rate.est.non.identified.events.pch[singlesource.leaks$tc_EquipmentGroupID[!identified.leaks] == "Tank"] <- tank.pch
rate.est.non.identified.events.pch[singlesource.leaks$tc_EquipmentGroupID[!identified.leaks] == "Wellhead.West"] <- wellhead.west.pch
rate.est.non.identified.events.pch[singlesource.leaks$tc_EquipmentGroupID[!identified.leaks] == "Separator.West"] <- separator.west.pch
rate.est.non.identified.events.pch[singlesource.leaks$tc_EquipmentGroupID[!identified.leaks] == "Wellhead.East"] <- wellhead.east.pch
rate.est.non.identified.events.pch[singlesource.leaks$tc_EquipmentGroupID[!identified.leaks] == "Separator.East"] <- separator.east.pch

y <- rate.est.identified.events/1000
x <- true.rates.identified.events/1000

this.fit <- lm(y ~ x + 0)

small.fit <- lm(y[x < 1] ~ x[x < 1] + 0)
small.fit

intercept <- this.fit$coefficients[1]
slope <- this.fit$coefficients[2]
r2 <- summary(this.fit)$r.squared

small.intercept <- small.fit$coefficients[1]
small.slope <- small.fit$coefficients[2]
small.r2 <- summary(small.fit)$r.squared


zoom.in <- 1

no.rate.factor <- 0.04

no.rate.est.mask <- is.na(rate.est.identified.events)





png(paste0(save.dir, 'parity.png'),
    res = 100, pointsize = 24, width = 1920, height = 1000)

par(mgp = c(2.25, 0.75, 0))
par(mar = c(3.5,2,1,1))
par(oma = c(0,1.5,0,0))
par(mfrow = c(1,2))

axis.max <- 10.5
no.rate.est.y <- -no.rate.factor*axis.max

################################################ START ZOOM OUT PLOT

plot(true.rates.identified.events/1000, rate.est.identified.events/1000, 
     asp = 1, 
     pch = rate.est.identified.events.pch, 
     col = rate.est.identified.events.cols,
     bg = rate.est.identified.events.cols,
     xlim = c(0,axis.max), ylim = c(no.rate.est.y,axis.max),
     ylab = "", las = 1,
     xlab = "METEC Emission Rate [kg/hr]")

segments(x0 = 0, y0 = 0, x1 = 20, y1 = 0, lwd = 3, col = "grey55") # false negative
segments(x0 = 0, y0 = 0, x1 = 0, y1 = 20, lwd = 3, col = "grey55") # false positive

points(true.rates.identified.events/1000, rate.est.identified.events/1000, 
       asp = 1, 
       pch = rate.est.identified.events.pch, 
       col = rate.est.identified.events.cols,
       bg = rate.est.identified.events.cols)

points(rep(0, n.false.positives), rate.est.all.events[false.positives]/1000,
       pch = '*', col= "black", cex = 2)

points(non.identified.leaks/1000, rep(0, n.non.identified.leaks),
       col = rate.est.non.identified.events.cols,
       bg = rate.est.non.identified.events.cols,
       pch = rate.est.non.identified.events.pch)

segments(x0 = true.rates.identified.events/1000,
         y0 = pmax(0,error.lower.identified.events/1000),
         y1 = error.upper.identified.events/1000,
         col = alpha(rate.est.identified.events.cols, 0.5),
         lwd = 2)

segments(x0 = true.rates.identified.events/1000 - true.rate.uncertainty.identified.events/1000,
         x1 = true.rates.identified.events/1000 + true.rate.uncertainty.identified.events/1000,
         y0 = rate.est.identified.events/1000,
         col = alpha(rate.est.identified.events.cols, 0.5),
         lwd = 2)

segments(x0 = 0, y0 = no.rate.est.y, x1 = 20, y1 = no.rate.est.y, lwd = 3,
         lty = 4, col = "gray55")

points(true.rates.identified.events[no.rate.est.mask]/1000, rep(no.rate.est.y, sum(no.rate.est.mask)),
       col = rate.est.identified.events.cols[no.rate.est.mask],
       pch = rate.est.identified.events.pch[no.rate.est.mask],
       bg = rate.est.identified.events.cols[no.rate.est.mask])

segments(x0 = 0, y0 = 0, x1 = 20, y1 = 20, lwd = 3) # 1:1

segments(x0 = 0, y0 = 0, x1 = 40, y1 = 20, lwd = 3, lty = 2) # factor of 2
segments(x0 = 0, y0 = 0, x1 = 20, y1 = 40, lwd = 3, lty = 2) # factor of 2

segments(x0 = 0, y0 = 0, x1 = 60, y1 = 20, lwd = 3, lty = 3) # factor of 3
segments(x0 = 0, y0 = 0, x1 = 20, y1 = 60, lwd = 3, lty = 3) # factor of 3

# segments(x0 = -intercept/slope, y0 = 0, x1 = (20-intercept)/slope, y1 = 20, lwd = 3, col = "magenta") # best fit
segments(x0 = 0, y0 = 0, x1 = (20)/intercept, y1 = 20, lwd = 3, col = "magenta") # best fit


rect(xleft = 0, ybottom = 0, xright = zoom.in, ytop = zoom.in, border = "black", lwd = 4)

# legend("topleft", c("1:1", "Factor of 2", "Factor of 3", "False positive & false negative", "No rate estimate", "Best fit"),
#        lwd = 5, cex = 1,
#        col = c("black", "black", "black", "gray55", "grey55", "magenta"),
#        lty = c(1,2,3,1,4,1),
# )

# legend("bottomright",
#        legend = c("Tanks", "West Wellhead", "West Separator", "East Wellhead", "East Separator", "False Positive"),
#        col = c(tank.color, wellhead.west.color, separator.west.color, wellhead.east.color, separator.east.color, "black"),
#        pt.bg = c(tank.color, wellhead.west.color, separator.west.color, wellhead.east.color, separator.east.color, NA),
#        pch = c(tank.pch, wellhead.west.pch, separator.west.pch, wellhead.east.pch, separator.east.pch, 42),
#        border = c("black", "black", "black", "black", "black", NA),
#        pt.cex = c(1,1,1,1,1,2))

mtext("Estimated Emission Rate [kg/hr]", side = 2, outer = T, line = 0.5)

################################################ START ZOOM IN PLOT

no.rate.est.y <- -no.rate.factor*zoom.in

plot(true.rates.identified.events/1000, rate.est.identified.events/1000, 
     asp = 1, las = 1,
     col = rate.est.identified.events.cols,
     bg = rate.est.identified.events.cols,
     pch = rate.est.identified.events.pch,
     xlim = c(0,zoom.in), ylim = c(no.rate.est.y,zoom.in),
     ylab = "", 
     xlab = "METEC Emission Rate [kg/hr]")

box(col = "black", lwd = 4)

segments(x0 = 0, y0 = 0, x1 = 20, y1 = 0, lwd = 3, col = "grey55") # false negative
segments(x0 = 0, y0 = 0, x1 = 0, y1 = 20, lwd = 3, col = "grey55") # false positive

points(true.rates.identified.events/1000, rate.est.identified.events/1000, 
       col = rate.est.identified.events.cols,
       bg = rate.est.identified.events.cols,
       pch = rate.est.identified.events.pch)

points(rep(0, n.false.positives), rate.est.all.events[false.positives]/1000,
       pch = '*', col= "black", cex = 2)

points(non.identified.leaks/1000, rep(0, n.non.identified.leaks),
       col = rate.est.non.identified.events.cols,
       pch = 19)

segments(x0 = true.rates.identified.events/1000,
         y0 = pmax(0,error.lower.identified.events)/1000,
         y1 = error.upper.identified.events/1000,
         col = alpha(rate.est.identified.events.cols, 0.5),
         lwd = 2)

segments(x0 = pmax(0,true.rates.identified.events/1000 - true.rate.uncertainty.identified.events/1000),
         x1 = true.rates.identified.events/1000 + true.rate.uncertainty.identified.events/1000,
         y0 = rate.est.identified.events/1000,
         col = alpha(rate.est.identified.events.cols, 0.5),
         lwd = 2)

segments(x0 = 0, y0 = no.rate.est.y, x1 = 20, y1 = no.rate.est.y, lwd = 3,
         lty = 4, col = "gray55")

points(true.rates.identified.events[no.rate.est.mask]/1000, rep(no.rate.est.y, sum(no.rate.est.mask)),
       col = rate.est.identified.events.cols[no.rate.est.mask],
       pch = rate.est.identified.events.pch[no.rate.est.mask],
       bg = rate.est.identified.events.cols[no.rate.est.mask])

segments(x0 = 0, y0 = 0, x1 = 20, y1 = 20, lwd = 3) # 1:1

segments(x0 = 0, y0 = 0, x1 = 40, y1 = 20, lwd = 3, lty = 2) # factor of 2
segments(x0 = 0, y0 = 0, x1 = 20, y1 = 40, lwd = 3, lty = 2) # factor of 2

segments(x0 = 0, y0 = 0, x1 = 60, y1 = 20, lwd = 3, lty = 3) # factor of 3
segments(x0 = 0, y0 = 0, x1 = 20, y1 = 60, lwd = 3, lty = 3) # factor of 3

# segments(x0 = -intercept/slope, y0 = 0, x1 = (20-intercept)/slope, y1 = 20, lwd = 3, col = "magenta")
segments(x0 = 0, y0 = 0, x1 = (20)/intercept, y1 = 20, lwd = 3, col = "magenta") # best fit
# segments(x0 = 0, y0 = 0, x1 = (20)/small.intercept, y1 = 20, lwd = 3, col = "magenta", lty = 2) # best fit

dev.off()



# RESULT 4.3: RATE ESTIMATE RATIOS
#---------------------------------------------------------------------------


# Create a sequence of factors
this.seq <- seq(1, max(1/rate.est.ratios, rate.est.ratios, na.rm = T) + 0.01, by = 0.01)

# Initialize vectors to hold the number of events whose ratio is below the corresponding factor in "this.seq"
to.save.all <- to.save.big <- to.save.small <- vector(length = length(this.seq))

# Loop through sequence of factors
for (i in 1:length(this.seq)){
  
  # Save number of emissions that have a ratio less than this factor
  to.save.all[i] <- sum(rate.est.ratios > 1/this.seq[i] & rate.est.ratios < this.seq[i], na.rm = T) / sum(!is.na(rate.est.ratios))
  
  # Save number of emissions > 1kg/hr that have a ratio less than this factor
  this.mask <- true.rates.identified.events > 1000
  to.save.big[i] <- sum(rate.est.ratios[this.mask] > 1/this.seq[i] & rate.est.ratios[this.mask] < this.seq[i], na.rm = T) / sum(!is.na(rate.est.ratios[this.mask]))
  
  # Save number of emissions < 1kg/hr that have a ratio less than this factor
  this.mask <- true.rates.identified.events <= 1000
  to.save.small[i] <- sum(rate.est.ratios[this.mask] > 1/this.seq[i] & rate.est.ratios[this.mask] < this.seq[i], na.rm = T) / sum(!is.na(rate.est.ratios[this.mask]))
}

# Set probabilities at which to print the empirical CDF value
probs <- c(1, 0.90, 0.75, 0.5)

# Initialize variables to hold the factors that occur at probs percent on the CDF
factors.all <- factors.big <- factors.small <- vector(length = length(probs))

# Loop through probabilities of interest
for (i in 1:length(probs)){
  factors.all[i] <- this.seq[which(to.save.all >= probs[i])[1]]
  factors.big[i] <- this.seq[which(to.save.big >= probs[i])[1]]
  factors.small[i] <- this.seq[which(to.save.small >= probs[i])[1]]
}

# Print number of emission events less than factor of 2 and 3 from truth
sum(rate.est.ratios > 1/2 & rate.est.ratios < 2, na.rm = T) / sum(!is.na(rate.est.ratios))
sum(rate.est.ratios > 1/3 & rate.est.ratios < 3, na.rm = T) / sum(!is.na(rate.est.ratios))


# Set plotting colors
under.est.col <- "#00524F"
over.est.col <- "#29A385"
line.col <- alpha("black", 0.475)
lwd.val <- 3
lty.val <- 3
axis.cex <- 0.82
label.line <- 1.8




png(paste0(save.dir, "over_under_est.png"),
    res = 100, pointsize = 24, width = 1920, height = 800)

par(mfrow= c(1,2))
par(mgp = c(2, 0.55, 0))
par(mar = c(1,0.5,1,0.5))
par(oma = c(2,2.5,2,0))

################################################ START < 1KG/HR HISTOGRAM

this.mask <- true.rates.identified.events <= 1000

rates.to.use <- rate.est.ratios[this.mask]

sum(rates.to.use > 1/2 & rates.to.use < 2, na.rm = T) / sum(!is.na(rates.to.use))
sum(rates.to.use > 1/3 & rates.to.use < 3, na.rm = T) / sum(!is.na(rates.to.use))

under.est.mask <- rates.to.use <= 1
under.est.mask[is.na(under.est.mask)] <- F

rates.to.use[under.est.mask] <- 1/rates.to.use[under.est.mask]

rates.to.plot <- rates.to.use
rates.to.plot[under.est.mask] <- -rates.to.plot[under.est.mask]

rates.to.plot[under.est.mask] <- rates.to.plot[under.est.mask] + 1
rates.to.plot[!under.est.mask] <- rates.to.plot[!under.est.mask] - 1


percent.of.events <- seq(1,100,1)/100

lower <- upper <- vector(length = length(percent.of.events))

for (i in 1:length(percent.of.events)){
  
  this.quant <- quantile(rates.to.plot, probs = c(0.5-percent.of.events[i]/2, 0.5+percent.of.events[i]/2), na.rm = T)
  
  lower[i] <- this.quant[1]
  upper[i] <- this.quant[2]
}



hist(rates.to.plot, breaks=seq(-6,6,1), col=over.est.col, xaxt = "n", xlab = "", 
     ylab = "", xpd = NA, ylim = c(-0.25,16.5), yaxt = "n", main = "")

mtext("Number of Emission Events", side = 2, line = label.line)

box()

abline(v = c(lower[50], upper[50]),
       col = line.col, lwd = lwd.val, lty = 1)

abline(v = c(lower[75], upper[75]),
       col = line.col, lwd = lwd.val, lty = 2)

abline(v = c(lower[90], upper[90]),
       col = line.col, lwd = lwd.val, lty = 3)

round(c(100*((-1/(lower[50]-1))-1), 100*(upper[50]+1-1)), 1)
round(c(100*((-1/(lower[75]-1))-1), 100*(upper[75]+1-1)), 1)
round(c(100*((-1/(lower[90]-1))-1), 100*(upper[90]+1-1)), 1)

round(c(lower[50]-1, upper[50]+1), 1)
round(c(lower[75]-1, upper[75]+1), 1)
round(c(lower[90]-1, upper[90]+1), 1)

tmp <- median(rates.to.plot,na.rm=T)
small.median.fd <- ifelse(tmp < 0, tmp-1, tmp+1)
print(round(small.median.fd, 2))

small.median.pd <- 100*ifelse(small.median.fd < 0, (-1/small.median.fd) - 1, small.median.fd)
print(round(small.median.pd, 1))

hist(rates.to.plot[under.est.mask], breaks=seq(-8,8,1), col=under.est.col, add=TRUE, main = "")
hist(rates.to.plot[!under.est.mask], breaks=seq(-8,8,1), col=over.est.col, add=TRUE, main = "")

segments(x0=median(rates.to.plot,na.rm=T), y0 = 0, y1 = 13, col = "#FB514B", lwd = lwd.val+2)

axis.factors <- c(-9:-2, 1:9)
converted.axis.factors <- ifelse(axis.factors < 0, -1/axis.factors, axis.factors)
axis.pdiff <- round(100 * (converted.axis.factors-1), 0)

axis(side = 1, at = seq(-8,8,1), labels = axis.factors, cex.axis = axis.cex)
axis(side = 2, at = seq(0,16,2), cex.axis = axis.cex)
axis(side = 3, at = seq(-8,8,2), labels = paste0(axis.pdiff[seq(1,length(axis.pdiff),2)],"%"), cex.axis = axis.cex)



################################################ START > 1KG/HR HISTOGRAM

this.mask <- true.rates.identified.events > 1000

rates.to.use <- rate.est.ratios[this.mask]

sum(rates.to.use > 1/2 & rates.to.use < 2, na.rm = T) / sum(!is.na(rates.to.use))
sum(rates.to.use > 1/3 & rates.to.use < 3, na.rm = T) / sum(!is.na(rates.to.use))

under.est.mask <- rates.to.use <= 1
under.est.mask[is.na(under.est.mask)] <- F

rates.to.use[under.est.mask] <- 1/rates.to.use[under.est.mask]

rates.to.plot <- rates.to.use
rates.to.plot[under.est.mask] <- -rates.to.plot[under.est.mask]

rates.to.plot[under.est.mask] <- rates.to.plot[under.est.mask] + 1
rates.to.plot[!under.est.mask] <- rates.to.plot[!under.est.mask] - 1

percent.of.events <- seq(1,100,1)/100

lower <- upper <- vector(length = length(percent.of.events))

for (i in 1:length(percent.of.events)){
  
  this.quant <- quantile(rates.to.plot, probs = c(0.5-percent.of.events[i]/2, 0.5+percent.of.events[i]/2), na.rm = T)
  
  lower[i] <- this.quant[1]
  upper[i] <- this.quant[2]
}



hist(rates.to.plot, breaks=seq(-8,8,0.25), col=over.est.col, xaxt = "n", xlab = "", 
     ylab = "", xlim = c(-2,2),
     ylim = c(-0.25,16.5), yaxt = "n", main = "")

box()

abline(v = c(lower[50], upper[50]),
       col = line.col, lwd = lwd.val, lty = 1)

abline(v = c(lower[75], upper[75]),
       col = line.col, lwd = lwd.val, lty = 2)

abline(v = c(lower[90], upper[90]),
       col = line.col, lwd = lwd.val, lty = 3)


round(c(100*((-1/(lower[50]-1))-1), 100*(upper[50]+1-1)), 1)
round(c(100*((-1/(lower[75]-1))-1), 100*(upper[75]+1-1)), 1)
round(c(100*((-1/(lower[90]-1))-1), 100*(upper[90]+1-1)), 1)

round(c(lower[50]-1, upper[50]+1), 1)
round(c(lower[75]-1, upper[75]+1), 1)
round(c(lower[90]-1, upper[90]+1), 1)


tmp <- median(rates.to.plot,na.rm=T)
big.median.fd <- ifelse(tmp < 0, tmp-1, tmp+1)
print(round(big.median.fd, 2))

big.median.pd <- 100*ifelse(big.median.fd < 0, (-1/big.median.fd) - 1, big.median.fd - 1)
print(round(big.median.pd, 1))


hist(rates.to.plot[under.est.mask], breaks=seq(-8,8,0.25), col=under.est.col, add=TRUE, main = "")
hist(rates.to.plot[!under.est.mask], breaks=seq(-8,8,0.25), col=over.est.col, add=TRUE, main = "")

segments(x0=median(rates.to.plot,na.rm=T), y0 = 0, y1 = 15, col = "#FB514B", lwd = lwd.val+2)


axis.factors <- c(seq(-3,-1.5,0.5), seq(1,3,0.5))
converted.axis.factors <- ifelse(axis.factors < 0, -1/axis.factors, axis.factors)
axis.pdiff <- round(100 * (converted.axis.factors-1), 0)

axis(side = 1, at = seq(-2,2,0.5), labels = axis.factors, cex.axis = axis.cex)
axis(side = 3, at = seq(-2,2,0.5), labels = paste0(axis.pdiff, "%"), cex.axis = axis.cex)


# legend("right", inset = c(0.1,0),
#        legend = c("Contains 50% of estimates",
#                   "Contains 75% of estimates",
#                   "Contains 90% of estimates",
#                   "Median factor difference"),
#        col = c(line.col, line.col, line.col, "#FB514B"),
#        lwd = 4, lty = c(1,2,3))

dev.off()


sum(rate.est.ratios >= 1/2 & rate.est.ratios <= 2, na.rm = T) / sum(!is.na(rate.est.ratios))
sum(rate.est.ratios >= 1/3 & rate.est.ratios <= 3, na.rm = T) / sum(!is.na(rate.est.ratios))





# RESULT 4.3: RATE ESTIMATE RATIOS VERSION 2!!!!!!
#---------------------------------------------------------------------------

# Set plotting colors
under.est.col <- "#00524F"
over.est.col <- "#29A385"
line.col <- alpha("black", 0.475)
lwd.val <- 3
lty.val <- 3
axis.cex <- 0.85
label.line <- 1.8




png(paste0(save.dir, "over_under_est.png"),
    res = 100, pointsize = 24, width = 1920, height = 800)

par(mfrow= c(1,2))
par(mgp = c(2, 0.55, 0))
par(mar = c(1,0.5,1,0.5))
par(oma = c(2,2.5,2,0))

################################################ START < 1KG/HR HISTOGRAM

this.mask <- true.rates.identified.events <= 1000

rates.to.plot <- rate.est.pdiffs[this.mask]

under.est.mask <- rates.to.plot <= 0
under.est.mask[is.na(under.est.mask)] <- F


round(100*sum(rates.to.plot<0, na.rm=T) / sum(!is.na(rates.to.plot)), 1)

percent.of.events <- seq(1,100,1)/100

lower <- upper <- vector(length = length(percent.of.events))

for (i in 1:length(percent.of.events)){
  
  this.quant <- quantile(rates.to.plot, probs = c(0.5-percent.of.events[i]/2, 0.5+percent.of.events[i]/2), na.rm = T)
  
  lower[i] <- this.quant[1]
  upper[i] <- this.quant[2]
}



hist(rates.to.plot, col=over.est.col,
     breaks = seq(-1,2.5, by = 0.25), 
     ylim = c(0, 16), 
     ylab = "", xlab = "",
     yaxt = "n", xaxt = "n", main = "")

axis(side = 1, at = seq(-1,2.5, by = 0.5), labels = paste0(100*seq(-1,2.5, by = 0.5),"%"),
     cex.axis = axis.cex)

axis(side = 2, at = seq(0,16, by = 2), cex.axis = axis.cex)

mtext("Number of Emission Events", side = 2, line = label.line)

box()

abline(v = c(lower[50], upper[50]),
       col = line.col, lwd = lwd.val, lty = 1)

# abline(v = c(lower[75], upper[75]),
#        col = line.col, lwd = lwd.val, lty = 2)

abline(v = c(lower[90], upper[90]),
       col = line.col, lwd = lwd.val, lty = 3)

round(100*c(lower[50], upper[50]), 1)
# round(100*c(lower[75], upper[75]), 0)
round(100*c(lower[90], upper[90]), 1)

lower <- ifelse(lower + 1 < 1, -1/(lower+1), lower + 1)
upper <- ifelse(upper + 1 < 1, -1/(upper+1), upper + 1)

round(c(lower[50], upper[50]), 1)
# round(c(lower[75], upper[75]), 2)
round(c(lower[90], upper[90]), 1)

small.median.pd <- mean(rates.to.plot,na.rm=T)
print(round(100*small.median.pd, 1))

100*median(rates.to.plot,na.rm=T)

hist(rates.to.plot[under.est.mask], col=under.est.col, breaks = seq(-1,3, by = 0.25), add=TRUE, main = "")
hist(rates.to.plot[!under.est.mask], col=over.est.col, breaks = seq(-1,3, by = 0.25), add=TRUE, main = "")

segments(x0 = small.median.pd, y0 = 0, y1 = 4, col = "#FB514B", lwd = lwd.val+2)

# axis.factors <- c(-3,-1,0,0.5,1,1.5,2,2.5,3)
# axis.factors <- ifelse(axis.factors < 0, axis.factors - 1, axis.factors + 1)
# axis.pdiff <- ifelse(axis.factors < 0, -1/axis.factors, axis.factors)
# axis.pdiff <- axis.pdiff - 1

# axis(side = 3, at = axis.pdiff, labels = paste0(axis.factors), cex.axis = axis.cex)




################################################ START > 1KG/HR HISTOGRAM

this.mask <- true.rates.identified.events > 1000

rates.to.plot <- rate.est.pdiffs[this.mask]

under.est.mask <- rates.to.plot <= 0
under.est.mask[is.na(under.est.mask)] <- F




round(100*sum(rates.to.plot<0, na.rm=T) / sum(!is.na(rates.to.plot)), 1)


percent.of.events <- seq(1,100,1)/100

lower <- upper <- vector(length = length(percent.of.events))

for (i in 1:length(percent.of.events)){
  
  this.quant <- quantile(rates.to.plot, probs = c(0.5-percent.of.events[i]/2, 0.5+percent.of.events[i]/2), na.rm = T)
  
  lower[i] <- this.quant[1]
  upper[i] <- this.quant[2]
}



hist(rates.to.plot, col=over.est.col,
     breaks = seq(-1,2.5, by = 0.25), 
     ylim = c(0,16),
     ylab = "", xlab = "",
     yaxt = "n", xaxt = "n", main = "")

axis(side = 1, at = seq(-1,2.5, by = 0.5), labels = paste0(100*seq(-1,2.5, by = 0.5),"%"),
     cex.axis = axis.cex)

box()

abline(v = c(lower[50], upper[50]),
       col = line.col, lwd = lwd.val, lty = 1)

# abline(v = c(lower[75], upper[75]),
#        col = line.col, lwd = lwd.val, lty = 2)

abline(v = c(lower[90], upper[90]),
       col = line.col, lwd = lwd.val, lty = 3)



round(100*c(lower[50], upper[50]), 1)
# round(100*c(lower[75], upper[75]), 0)
round(100*c(lower[90], upper[90]), 1)

lower <- ifelse(lower + 1 < 1, -1/(lower+1), lower + 1)
upper <- ifelse(upper + 1 < 1, -1/(upper+1), upper + 1)

round(c(lower[50], upper[50]), 1)
# round(c(lower[75], upper[75]), 2)
round(c(lower[90], upper[90]), 1)



large.median.pd <- mean(rates.to.plot,na.rm=T)
print(round(100*large.median.pd, 1))

100*median(rates.to.plot,na.rm=T)

hist(rates.to.plot[under.est.mask], col=under.est.col, breaks = seq(-1,3, by = 0.25), add=TRUE, main = "")
hist(rates.to.plot[!under.est.mask], col=over.est.col, breaks = seq(-1,3, by = 0.25), add=TRUE, main = "")

segments(x0 = large.median.pd, y0 = 0, y1 = 14, col = "#FB514B", lwd = lwd.val+2)

# axis.factors <- c(-3,-1,0,0.5,1,1.5,2,2.5,3)
# axis.factors <- ifelse(axis.factors < 0, axis.factors - 1, axis.factors + 1)
# axis.pdiff <- ifelse(axis.factors < 0, -1/axis.factors, axis.factors)
# axis.pdiff <- axis.pdiff - 1
# 
# axis(side = 3, at = axis.pdiff, labels = paste0(axis.factors), cex.axis = axis.cex)

# legend("right", inset = c(0.1,0),
#        legend = c("Contains 50% of estimates",
#                   "Contains 75% of estimates",
#                   "Contains 90% of estimates",
#                   "Average percent difference"),
#        col = c(line.col, line.col, line.col, "#FB514B"),
#        lwd = 4, lty = c(1,2,3))

dev.off()




# RESULT 4.4: CUMULATIVE RESULTS
#---------------------------------------------------------------------------

# Vector to hold true emission rate on minute-by-minute frequency
true.emission.rate <- vector(length = length(times)) # kg/hr

for (i in 1:nrow(singlesource.leaks)){
  
  # Get start and stop times of this leak
  leak.start <- singlesource.leaks$tc_ExpStartDatetime[i]
  leak.end   <- singlesource.leaks$tc_ExpEndDatetime[i]
  
  # Create minute sequence spanning this leak
  these.leak.times <- seq(leak.start, leak.end, by = "1 min")
  these.leak.times <- round_date(these.leak.times, unit = "min")
  
  # Mask in times during this leak
  this.mask <- times %in% these.leak.times
  
  # Apply rate in kg/hr to vector
  true.emission.rate[this.mask] <- singlesource.leaks$tc_C1MassFlow[i]/1000
}

# Convert from rate to mass (60 minutes per hour, and then each element is one minute long)
true.emissions.kg <- true.emission.rate/60 # kg

# Vector to hold estimated emission rate on minute-by-minute frequency
pred.emission.rate <- vector(length = length(times)) # kg/hr

for (i in 1:n.ints){
  
  # Get start and stop times of this event
  event.start <- as_datetime(event.start.times[i], tz = "America/Denver")
  event.end <- as_datetime(event.end.times[i], tz = "America/Denver")
  
  # Create minute sequence spanning this event
  these.event.times <- seq(event.start, event.end, by = "1 min")
  
  # Mask in times during this event
  this.mask <- times %in% these.event.times
  
  # Apply rate in kg/hr to vector
  pred.emission.rate[this.mask] <- rate.est.all.events[i]/1000
}


# Convert from rate to mass (60 minutes per hour, and then each element is one minute long)
pred.emissions.kg <- pred.emission.rate/60 # kg

# Set NAs (when we identify an event but don't produce rate estimate) to zero
na.times <- is.na(pred.emission.rate)
true.emissions.kg[na.times] <- 0 
pred.emissions.kg[na.times] <- 0
true.emission.rate[na.times] <- 0



png(paste0(save.dir, "cumulative.png"),
    width = 1920, height = 700, res = 100, pointsize = 30)

layout.matrix <- t(matrix(c(1, 1, 1, 1, 2, 2, 2)))

layout(mat = layout.matrix,
       heights = c(1),
       widths  = c(1)) 

par(mar = c(4,3.5,2,1))
par(mgp = c(2.4, 0.8, 0))
# par(mgp = c(2.4,2, 0))

this.order <- order(true.emission.rate)

plot(true.emission.rate[this.order], cumsum(true.emissions.kg[this.order]), type = "l",
     ylim = c(-100,420), lwd = 4, col = "gray25",
     ylab = "Cumulative Emissions [kg]", xlab = "True Emission Rate [kg/hr]")

abline(h = 0, lwd = 2)
abline(v = c(0,1), lty = 2, col = "gray55", lwd = 2)

lines(true.emission.rate[this.order], cumsum(pred.emissions.kg[this.order]), 
      col = "gray75", lwd = 4)

diff <- cumsum(pred.emissions.kg[this.order]) - cumsum(true.emissions.kg[this.order])

pdiff <- diff / cumsum(true.emissions.kg[this.order])

lines(true.emission.rate[this.order], diff, col = "#FB514B", lwd = 4)


 legend("topright", c("Truth", "Estimate", "Difference"),
       lty = 1, lwd = 5, col = c("gray25", "gray75", "#FB514B"), bty = "n",
       inset = c(0,0.5))


fp.mask <- true.emission.rate == 0
small.mask <- true.emission.rate <= 1 & true.emission.rate > 0
big.mask <- true.emission.rate > 1 


sum.true.fp <- sum(true.emissions.kg[fp.mask])
sum.pred.fp <- sum(pred.emissions.kg[fp.mask])

sum.true.small <- sum(true.emissions.kg[small.mask])
sum.pred.small <- sum(pred.emissions.kg[small.mask])

sum.true.big <- sum(true.emissions.kg[big.mask])
sum.pred.big <- sum(pred.emissions.kg[big.mask])

sum.true <- sum(true.emissions.kg[small.mask | big.mask])
sum.pred <- sum(pred.emissions.kg[small.mask | big.mask])

sum.true.fp.included <- sum(true.emissions.kg[small.mask | big.mask | fp.mask])
sum.pred.fp.included <- sum(pred.emissions.kg[small.mask | big.mask | fp.mask])

to.plot <- rbind(c(sum.true.fp, sum.true.small, sum.true.big, sum.true, sum.true.fp.included),
                 c(sum.pred.fp, sum.pred.small, sum.pred.big, sum.pred, sum.pred.fp.included))

b <- barplot(to.plot, beside = T, ylab = "Cumulative Emissions [kg]",
             # names.arg = c("0\nkg/hr\n(FP)", "(0,1]\nkg/hr\n(Small)", "> 1\nkg/hr\n(Large)", "Total:\nFP\nExcluded", "Total:\nFP\nIncluded"),
             names.arg = rep("", 5),
             col = c("gray25", "gray75"), ylim = c(0,420))


legend("topleft", c("Truth", "Estimate"),
       fill = c("gray25", "gray75"), bty = "n")


dev.off()


round(100*(sum.pred.small - sum.true.small)/sum.true.small, 1)

round(100*(sum.pred.big - sum.true.big)/sum.true.big, 1)

round(100*(sum.pred - sum.true)/sum.true, 1)

round(100*(sum.pred.fp.included - sum.true.fp.included)/sum.true.fp.included, 1)




big.weight <- sum.true.big / sum.true

small.weight <- sum.true.small / sum.true

fp.weight <- sum.pred.fp

sum.true - sum.pred.fp

fp.ratio <- sum.pred.fp / sum.true

round(100*fp.ratio, 1)


round(weighted.mean(c(-38.5, -0.2), c(small.weight, big.weight)), 1)

round(weighted.mean(c(-38.5, -0.2), c(small.weight, big.weight)) + 100*fp.ratio, 1)



round(weighted.mean(c(-38.5, -0.2), c(1-0.94, 0.94)), 1)

round(weighted.mean(c(-38.5, -0.2), c(1-0.94, 0.94)) + 100*fp.ratio, 1)



# RESULT 5: DEPENDENCE
#---------------------------------------------------------------------------

# Initialize vectors to hold covariates of interest
leak.durations.identified.events <- ws.identified.events <-
  wd.identified.events <- wd.sd.identified.events <- vector(length = n.identified.leaks)

count <- 1

# Loop through identified events and save covariates
for (i in which(identified.leaks)){
  
  # Emission duration
  leak.start.time <- singlesource.leaks$tc_ExpStartDatetime[i]
  leak.end.time <- singlesource.leaks$tc_ExpEndDatetime[i]
  leak.durations.identified.events[count] <- difftime(leak.end.time, leak.start.time, unit = "hours")
  
  # Mask in emission times
  this.int <- interval(leak.start.time, leak.end.time)
  this.mask <- times %within% this.int
  
  # Wind speed
  ws.identified.events[count] <- mean(data$WS[this.mask], na.rm = T)
  
  # Circular variance of WD
  these.wd <- data$WD[this.mask]
  wd.x <- cos(these.wd)
  wd.y <- sin(these.wd)
  wd.identified.events[count] <- atan2(mean(wd.y, na.rm = T), mean(wd.x, na.rm = T))
  # wd.sd.identified.events[count] <- var.circular(as.circular(these.wd,
  #                                                            units = "radians",
  #                                                            type = "angles",
  #                                                            template = "none",
  #                                                            zero = 0,
  #                                                            rotation = "counter",
  #                                                            modulo = "asis"), na.rm = T)
  
  count <- count + 1
}

# Set colors and points for plots
est.within.2 <- rate.est.ratios > 1/2 & rate.est.ratios < 2
est.within.3 <- rate.est.ratios > 1/3 & rate.est.ratios < 3

these.col <- rep("firebrick3", n.identified.leaks)
these.col[is.na(rate.est.ratios)] <- "black"
these.col[est.within.3] <- "goldenrod2"
these.col[est.within.2] <- "forestgreen"

these.pch <- rep(15, n.identified.leaks)
these.pch[is.na(rate.est.ratios)] <- 18
these.pch[est.within.3] <- 17
these.pch[est.within.2] <- 16

cex.val <- 1.2
ylim.max <- 7.2


# Mask in events with relatively similar emission rates (all less than 1 kg/hr)
this.mask <- true.rates.identified.events < 1000

# Get relative error for these events
y <- rate.est.ratios[this.mask]
y <- ifelse(y<1, 1/y, y)


# 
# png(paste0(save.dir, "dependence.png"),
#     res = 100, pointsize = 30, width = 1920, height = 1080*1.25)
# 
# par(mfrow = c(3,2))
# par(mar = c(3.25,3.25,1,1))
# par(mgp = c(2.25, 0.75, 0))
# par(oma = c(0,0,0.25,0))
# 
# 
# xlim.vals <- range(leak.durations.identified.events[this.mask])
# 
# plot(leak.durations.identified.events, true.rates.identified.events/1000,
#      col = these.col, pch = these.pch, cex = cex.val,
#      xlim = c(0, max(leak.durations.identified.events)),
#      ylab = "METEC Emission Rate [kg/hr]",
#      xlab = "METEC Emission Duration [hr]",
#      xpd = NA)
# 
# rect(xleft = xlim.vals[1] - 0.25, ybottom = 0 - 0.05,
#      xright = xlim.vals[2] + 0.25, ytop = 1 + 0.25,
#      border = "black", lwd = 4,
#      xpd = NA)
# 
# # legend("topright",
# #        legend = c("Est. within Factor of 2",
# #                   "Est. within Factor of 3",
# #                   "Est. outside Factor of 3",
# #                   "Quantification Unavailable"),
# #        col = c("forestgreen", "goldenrod2", "firebrick3", "black"),
# #        pch = c(16, 17, 15, 18),
# #        pt.cex = cex.val)
# 
# plot(leak.durations.identified.events[this.mask], y,
#      xlim = xlim.vals,
#      ylim = c(1,ylim.max),
#      col = these.col[this.mask], pch = these.pch[this.mask], cex = cex.val,
#      ylab = "Abs(Factor Difference)",
#      xlab = "METEC Emission Duration [hr]")
# 
# x <- leak.durations.identified.events[this.mask]
# 
# cor.val <- cor(x,y, use = "complete.obs", method = "pearson")
# 
# abline(lm(y~x), col = "magenta", lwd = 3)
# 
# box(col = "black", lwd = 4)
# 
# mtext(paste0("correlation = ", round(cor.val,2)),
#       adj = 0.025, cex = 0.8, line = -1.25)
# 
# 
# 
# #####
# 
# xlim.vals <- range(ws.identified.events[this.mask])
# 
# plot(ws.identified.events, true.rates.identified.events/1000,
#      col = these.col, pch = these.pch, cex = cex.val,
#      xlim = c(0, max(ws.identified.events)),
#      ylab = "METEC Emission Rate [kg/hr]",
#      xlab = "Average Wind Speed During Event [m/s]",
#      xpd = NA)
# 
# rect(xleft = xlim.vals[1] - 0.25, ybottom = 0 - 0.05,
#      xright = xlim.vals[2] + 0.25, ytop = 1 + 0.25,
#      border = "black", lwd = 4,
#      xpd = NA)
# 
# plot(ws.identified.events[this.mask],  y,
#      xlim = xlim.vals,
#      ylim = c(1,ylim.max),
#      col = these.col[this.mask], pch = these.pch[this.mask], cex = cex.val,
#      ylab = "Abs(Factor Difference)",
#      xlab = "Average Wind Speed During Event [m/s]")
# 
# box(col = "black", lwd = 4)
# 
# x <- ws.identified.events[this.mask]
# 
# cor.val <- cor(x,y, use = "complete.obs", method = "pearson")
# 
# abline(lm(y~x), col = "magenta", lwd = 3)
# 
# box(col = "black", lwd = 4)
# 
# mtext(paste0("correlation = ", round(cor.val,2)),
#       adj = 0.025, cex = 0.8, line = -1.25)
# 
# 
# #####
# 
# 
# xlim.vals <- range(wd.sd.identified.events[this.mask])
# 
# plot(wd.sd.identified.events, true.rates.identified.events/1000,
#      col = these.col, pch = these.pch, cex = cex.val,
#      xlim = c(0, 1),
#      ylab = "METEC Emission Rate [kg/hr]",
#      xlab = "Circular Variance of Wind Direction During Event",
#      xpd = NA)
# 
# rect(xleft = xlim.vals[1] - 0.025, ybottom = 0 - 0.05,
#      xright = xlim.vals[2] + 0.025, ytop = 1 + 0.25,
#      border = "black", lwd = 4,
#      xpd = NA)
# 
# plot(wd.sd.identified.events[this.mask], y,
#      xlim = xlim.vals,
#      ylim = c(1,ylim.max),
#      col = these.col[this.mask], pch = these.pch[this.mask], cex = cex.val,
#      ylab = "Abs(Factor Difference)",
#      xlab = "Circular Variance of Wind Direction During Event")
# 
# box(col = "black", lwd = 4)
# 
# x <- wd.sd.identified.events[this.mask]
# 
# cor.val <- cor(x,y, use = "complete.obs", method = "pearson")
# 
# abline(lm(y~x), col = "magenta", lwd = 3)
# 
# box(col = "black", lwd = 4)
# 
# mtext(paste0("correlation = ", round(cor.val,2)),
#       adj = 0.025, cex = 0.8, line = -1.25)
# 
# 
# dev.off()
