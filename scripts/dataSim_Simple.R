dataSim_Simple <- function(t_end, t_start = 0, replacement = TRUE, nFrames = 100, samples = NULL){
  
  # set up ----
  if(is.null(samples)){
  # set the controlling parameters
    sampSize <- nFrames * 99
    
    beta <- 1.15
    eta <- 5.253
    
    # simulate the data generating process ----
    # generate samples from weibull
    samples <- rweibull(sampSize, shape = beta, scale = eta)
    
  } else {  # Also allow user to pass in samples so that different observations of 
            # the same data can be achieved
    
    errorCondition("sample size must be equal to 100 x nFrames", length(samples) != 100 * nFrames)
    
  }
  
  # arrange the samples into a 100 x 99 matrix
  frameLifetimes <- matrix(samples, 
                           nrow = nFrames, 
                           byrow = TRUE)
  
  # take the cumulative sum across each row
  # this simulates 100 frames with 99 failures where time between failures
  # follows a weibull distribution
  frameFailureTimes <- apply(frameLifetimes, 1, cumsum) %>% t()
  
  # if no repetitive replacement of frames then data set is simple type one 
  # censoring so only need first column in matrix
  if(!replacement){
    
    dfType1 <- data.frame(lifetime = frameFailureTimes[, 1], 
                          cens = as.numeric(frameFailureTimes[, 1] > t_end))
    
    dfType1$lifetime[dfType1$cens == 1] <- t_end
    
    dfObs <- dfType1
    
  } else { 
    # now simulate observation of failures within window t_start -> t_end ----
    # find the failures in each row which fall within observation window
    observedFailureInd <- lapply(1:nFrames, 
                                 function(i) {
                                   which((frameFailureTimes[i, ] >= t_start) & 
                                           (frameFailureTimes[i, ] <= t_end))
                                 })
    
    # function to go get observed lifetimes for a frame and arrange into data frame
    extractFrameWindow <- function(dfRow, rowNum){
      
      startCensTrue <- ifelse(t_start == 0, 0, 1)
      
      if(length(dfRow) == 0){
        df <- data.frame(frame = rowNum, 
                         lifetime = t_end - t_start,
                         startCens = startCensTrue, 
                         endCens = 1, 
                         cens = 1)
      } else {
        failureTimes <- c(t_start, frameFailureTimes[rowNum, dfRow], t_end)
        censFrameLifetimes <- diff(failureTimes)
        numberOfLifetimes <- length(censFrameLifetimes)
        
        df <- data.frame(frame = rep(rowNum, numberOfLifetimes),
                         lifetime = censFrameLifetimes,
                         startCens = c(startCensTrue, 
                                       rep(0, numberOfLifetimes - 1)),
                         endCens = c(rep(0, numberOfLifetimes - 1), 
                                     1),
                         cens = rep(0, numberOfLifetimes))
        df$cens[(df$startCens == 1) | (df$endCens ==1)] <- 1 
      }
      return(df)
    }
    
    # apply function to all frames and colaps list into large data frame
    dfObsProcess <- lapply(1:nFrames, 
                           function(i) {
                             extractFrameWindow(observedFailureInd[[i]], i)
                           }) %>% bind_rows()
    dfObs <- dfObsProcess
  }
  
  return(dfObs)
  
}
