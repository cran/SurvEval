stCPE <- function(Time, Event, Markers, starttime, tau){
  
  if(any(is.na(Time))| any(is.na(Event))) stop("Please remove missing values in your data")
  
  if(starttime > max(Time,na.rm=TRUE)) stop("Error: start time beyond max observation time")
  if(tau > max(Time,na.rm=TRUE)) warning("End time beyond observation times; set to observed uncensored max")
  
  if(tau > max(Time)) tau <- max(Time[Event==1])
  
  mod1 <- coxph(Surv(Time, Event)~Markers) 
  
  Sstart <-  summary(survfit(mod1,data.frame(Markers)),  t=starttime)$surv
  Send <-  summary(survfit(mod1,data.frame(Markers)),  t=tau)$surv
  Xbetahat <- mod1$linear.predictors
  
  stCPE.estimate <- sstartconcordance(Xbetahat , Sstart, Send)
  
  return(stCPE.estimate)
}

