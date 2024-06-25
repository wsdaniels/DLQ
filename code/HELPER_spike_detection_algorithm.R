find.spikes <- function(
    times, obs,
    going.up.threshold = 0.25, # [ppm]
    return.threshold = 5, # [%], must be on (0,100)
    amp.threshold = 1, # [ppm]
    cont.diff.threshold = 0.25, # [ppm]
    cont.diff.num = 10, # [number of observations]
    make.plot = F # include a plot with the returned data?
){
  
  
  # Prep vector to hold event indices
  events <- rep(NA, length(obs))
  
  # Count is increased every time there is a new event
  count <- 0
  
  # Set the "in event" flag to false before looping through data
  in.event <- F
  
  # Skip leading NAs
  start.ind <- min(which(!is.na(obs))) + 1
  
  # Indicator for if loop has reached the last observation
  last.ob <- F
  
  # Running background estimate
  background <- NA
  
  # Loop through the observations
  for (i in start.ind:length(obs)){
    
    # Is this the last observation?
    if (i == length(obs)){last.ob <- T}
    
    # update current index and last index to deal with NA gaps
    if (is.na(obs[i]) & !is.na(obs[i-1])){
      last.ind <- i-1
      next
    } else if (is.na(obs[i]) & is.na(obs[i-1])){
      next
    } else if (!is.na(obs[i]) & is.na(obs[i-1])){
      current.ind <- i
    } else {
      current.ind <- i
      last.ind <- i-1
    }
    
    # Branch for when you are not currently in an event
    if (in.event == F){
      
      # Compute difference between this observation and the last
      current.diff <- obs[current.ind] - obs[last.ind]
      
      # If on last ob, only enter event if difference is greater than
      # both going.up.threshold and amp.threshold, as there will be no
      # other check for amplitude. Otherwise, just use going.up.threshold.
      threshold.to.use <- ifelse(last.ob,
                                 max(going.up.threshold, amp.threshold),
                                 going.up.threshold)
      
      # If difference is greater than threshold value, enter an event
      if (current.diff > threshold.to.use){
        in.event <- T
        count <- count + 1
        event.obs <- obs[current.ind]
        events[current.ind] <- count
        background <- obs[last.ind]
      }
      
      # Branch for when you are currently in an event
    } else {
      
      # Compute the maximum value of this event up to time step i and with
      # the background removed. We assume that the last observation before
      # the large difference is the background concentration.
      current.max <- max(event.obs) - background
      
      # Compute the current background corrected observation.
      current.ob <- obs[current.ind] - background
      
      # End event if the current observation has returned to return.threshold 
      # percent of the background corrected maximum value or if it is the 
      # last observation
      if ((current.ob < 2*background & current.ob < return.threshold * current.max / 100) | last.ob){
        
        # End event
        in.event <- F
        
        # Fill in any gaps due to NAs in observations
        event.seq <- seq(min(which(events == count), na.rm=T),
                         max(which(events == count), na.rm=T))
        events[event.seq] <- count
        
        # Compute background corrected event amplitude.
        # Use mean of first and last observation as the background estimate,
        # unless the event includes the last observation, than use just first ob
        if (last.ob){
          event.size <- max(event.obs) - background
        } else {
          event.size <- max(event.obs) - mean(c(background, obs[current.ind]))
        }
        
        # If the background corrected amplitude is less than amplitude threshold,
        # then remove this event
        if (event.size < amp.threshold){
          events[events == count] <- NA
          count <- count - 1
        }
        
        # Branch for if the event is not ended
      } else {
        
        # Add current observation to vector of observations in this event
        event.obs <- c(event.obs, obs[current.ind])
        
        # Add this event number to the record of events
        events[i] <- count
        
        # Check for long stretch of small differences. This indicates that the
        # return threshold was not triggered but should have been.
        if (length(event.obs) > cont.diff.num){
          
          # Compute first index of observations within this event that should
          # be considered in this check. Only consider cont.diff.num observations.
          window.start <- length(event.obs) - cont.diff.num
          
          # Grab the observations that should be considered
          obs.in.window <- event.obs[window.start:length(event.obs)]
          
          # If the absolute value of the differences are all below 
          # cont.diff.threshold, then end the event here.
          if (all(abs(diff(obs.in.window)) < cont.diff.threshold)){
            
            # End event
            in.event <- F
            
            # Fill in any gaps due to NAs in observations
            event.seq <- seq(min(which(events == count), na.rm=T),
                             max(which(events == count), na.rm=T))
            events[event.seq] <- count
            
            # Compute background corrected event amplitude
            # Use mean of first and last observation as the background estimate
            event.size <- max(event.obs) - mean(c(background, obs[current.ind]))
            
            # If the background corrected amplitude is less than amplitude threshold,
            # then remove this event
            if (event.size < amp.threshold){
              events[events == count] <- NA
              count <- count - 1
              
              # Otherwise, remove the observations that are not changing
            } else {
              events[seq(current.ind-cont.diff.num+1, current.ind)] <- NA
            }
            
          } # end if (all(abs(diff(obs.in.window)) < cont.diff.threshold))
        } # end if (length(event.obs) > cont.diff.num)
      } # end if (current.ob < return.threshold * current.max / 100)
    } # end if (in.event == F)
  } # end loop through observations
  
  
  # Save everything
  filtered.events <- data.frame(time = times,
                                events = events)
  
  
  
  if (make.plot){
    
    par(mar = c(3,4,1,1))
    
    plot(times, obs, type = "l", xlab = "",
         lwd = 2,
         ylab = "Methane [ppm]",
         ylim = c(0, max(obs, na.rm = T)))
    
    if (!all(is.na(events))){
      
      events <- filtered.events$events
      
      event.nums <- unique(na.omit(events))
      
      col.vals <- rep(c("red", "blue", "purple", "green"), length(event.nums))
      
      for (i in 1:length(event.nums)){
        this.spike <- which(events == event.nums[i])
        points(times[this.spike], obs[this.spike], pch = 19, col = col.vals[i])
        lines(times[this.spike], obs[this.spike], col = col.vals[i], lwd = 2)
      }
      
    }
  }
  
  return(filtered.events)
  
}