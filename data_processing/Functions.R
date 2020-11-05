print("The functions were succesfully loaded")

substrRight <- function(x, n){
  sapply(x, function(xx)
    substr(xx, (nchar(xx)-n+1), nchar(xx))
  )
}

######################
#
#   Calculate ppm and ppm difference between two masses
#
#####################
ppmCal<-function(run,ppm) {
  return((run*ppm)/1000000)
}

ppmDiff<-function(mz1,mz2) {
  return(((mz1-mz2)/mz1)*1000000)
}

######################
#
#   Grep name in x
#
#####################
specgrep <- function(x,name) {
  x=x[,grep(name, names(x))]
  return(x)
}


######################
#
#   Remove name in x
#
#####################
removegrep <- function(x,name) {
  x=x[,-grep(name, names(x))]
  return(x)
}


######################
#
#   BlankFilter - find features which are "highly" present in blank (based on a ratio cutoff)
#   to.remove=BlankFilter(Blanks,samples,0.01)
#
#####################
BlankFilter <- function(blanks, samples, cutoff) {
  blanks[is.na(blanks)] <- 0
  samples[is.na(samples)] <- 0

  blanks <- apply(blanks,1,median,na.rm=TRUE)
  samples <- apply(samples,1,max,na.rm=TRUE)

  to.remove <- which(blanks/samples >= cutoff)
  return(to.remove)
}

######################
#
#   Outersect
#
#####################
outersect <- function(x, y) {
  sort(c(x[!x%in%y],
         y[!y%in%x]))
}

#####################
#
#   Normalize for runorder effects e.g. sensetivity decay
#   data.norm <- order.norm(data, span)
#
#####################
order.norm <- function(x, s) {
  n<-c(1:ncol(x))
  for(j in 1:nrow(x)){
    fit<-loessFit(as.numeric(x[j,]),n, span=s)
    x[j,] <-fit$residuals+mean(as.numeric(x[j,]), na.rm=TRUE)
  }
  return(x)
}
