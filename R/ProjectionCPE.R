ProjectionCPE <- function(Time, Event,StandardMarkers,NewMarkers,  tau , Block=TRUE){
  
  NF <- length(Time)
  
  if(tau > max(Time,na.rm=TRUE)) warning("tau beyond observation times; set to observed uncensored max")
  if(tau > max(Time)) endtime <- max(Time[Event==1])
  
  if(is.vector(StandardMarkers)){
    StandardMarkers <- matrix(StandardMarkers)
    colnames(StandardMarkers) <- "StandardMarkers"
  }else{
    colnames(StandardMarkers) <- paste("SM", 1:dim(StandardMarkers)[2], sep="")
    
  }
  if(is.vector(NewMarkers)){
    NewMarkers <- matrix(NewMarkers)
    colnames(NewMarkers) <- "NewMarkers"
  }else{
    colnames(NewMarkers) <- paste("NM", 1:dim(NewMarkers)[2], sep="")
    
  }
  
  if( any(colnames(StandardMarkers) %in% colnames(NewMarkers)  )) stop("duplicate column names: please fix")
  
  datacombined <- data.frame(Time,Event,StandardMarkers,NewMarkers)
  mod1 <- coxph(Surv(Time, Event)~., data=datacombined)
  
  cumbasehazard <-  basehaz(mod1, centered=FALSE)
  cumbasehazardTau <-  max(cumbasehazard$hazard[cumbasehazard$time < tau])
  
  BX1.2 <- StandardMarkers %*% mod1$coefficients[colnames(StandardMarkers)]
  BZ1.2 <- NewMarkers %*% mod1$coefficients[colnames(NewMarkers)]
  
  BX1 <- BX1.2[order(BX1.2)]
  BZ1 <- BZ1.2[order(BX1.2)]
  
  bMatrix <- matrix(rep(BX1, NF), ncol=NF, byrow=FALSE)-matrix(rep(BX1, NF), ncol=NF, byrow=TRUE)
  gMatrix <- matrix(rep(BZ1,NF), ncol=NF, byrow=FALSE)-matrix(rep(BZ1, NF), ncol=NF, byrow=TRUE)
  
  sdBX1 =stats::sd(BX1)
  h = sqrt(2)*sdBX1*length(BX1)^(-1/5)    
  bKernel <- (exp(-.5*((bMatrix/h)^2))/sqrt(2*pi))/h
  
  
  aa <- matrix(rep(BX1, each= length(BX1)), nrow=length(BX1), byrow=TRUE)
  bb <- matrix(rep(BZ1, each=length(BZ1)), nrow=length(BZ1), byrow=FALSE)
  
  SurvMatrix <- exp(-1*cumbasehazardTau*exp(aa+bb))
  
  if(Block==FALSE){
  reducedCPENum <-  reducedCPE( bMatrix,  gMatrix,   bKernel, SurvMatrix)
  
  SurvDiag <- diag(SurvMatrix)
  aa <- matrix(rep(SurvDiag, each= length(SurvDiag)), nrow=length(SurvDiag), byrow=TRUE)
  bb <- matrix(rep(SurvDiag, each=length(SurvDiag)), nrow=length(SurvDiag), byrow=FALSE)
  cc <- 1-aa*bb
  diag(cc) <- 0 
  reducedCPE <-  reducedCPENum/(0.5*sum(cc))
  }else{
  
  blocknum <- floor(dim(bMatrix)[1]/25)
  randomgroups <- sample(1:dim(bMatrix)[1], 25*blocknum, replace=FALSE)
  
  reducedCPENum <- 0 
  denomsum <- 0 
  for(bb in 1:blocknum){
    subgrp <- randomgroups[((bb-1)*25+1):(bb*25)]
    
    subgrp <- sort(subgrp)
    
    reducedCPENum <- reducedCPENum+ reducedCPE( bMatrix[subgrp,subgrp],  gMatrix[subgrp,subgrp],   bKernel[subgrp,subgrp], SurvMatrix[subgrp,subgrp])
    
    
    SurvDiag <- diag(SurvMatrix[subgrp,subgrp])
    aa <- matrix(rep(SurvDiag, each= length(SurvDiag)), nrow=length(SurvDiag), byrow=TRUE)
    bb <- matrix(rep(SurvDiag, each=length(SurvDiag)), nrow=length(SurvDiag), byrow=FALSE)
    cc <- 1-aa*bb
    diag(cc) <- 0 
    
    denomsum <- denomsum+sum(cc)
  }
  
  reducedCPE <-  reducedCPENum/(0.5*denomsum)
  }
  
  
  
  return(list(projCPE=reducedCPE))
}


