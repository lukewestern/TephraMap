#' Estimate the isopachs and total volume given tephra
#' samples.
#'
#' @param fn String of the csv file containing the sampled tephra deposits. This must have columns containing the X coordinate, Y coordinate and thickness of the tephra deposit. There must be headers for all columns, with the sampled thickness having the units as the header: <mm>, <cm> or <m> (without the < >).
#' @param coordsys String of the coordinate system used. This must be either \code{"UTM"}, which is the UTM coordinates (where the UTM zone must also be given); \code{"dfs"}, which is distance from some origin in kilometres; or \code{"latlon"} for WG84 latitude and longitude, where the UTM zone must also be supplied as units are converted to UTM if \code{super=FALSE}.
#' @param UTMtz (optional) The UTM zone if using UTM or latlon coordinates (and \code{super==FALSE})
#' @param super (optional, default=\code{FALSE}) Set to true if this is a super eruption. The calculation is then done on a globe rather than a flat plain. Coordinates must be in latlon. (This only works for super eruptions, i.e. those spanning several degrees longitude/latitude).
#' @param prior.range (optional, default=\code{NULL}) Prior probability for the spatial range of the tephra field in km using a penalised complexity prior. This is a vector of length 2, \code{c(range0, Prange)}, where Pr(range < range0) = Prange. Setting Prange as \code{NA} fixes range to range0. Leave as \code{NULL} for the default prior, approximated from the measurement field.
#' @param prior.sigma (optional, default=\code{NULL}) Prior probability for the standard deviation of the log-tickness of the tephra field defined using a penalised complexity prior. This is a vector of length 2, \code{c(sigma0, Psigma)}, where Pr(sigma > sigma0) = Psigma. Setting Psigma as \code{NA} fixes sigma to sigma0. Leave as \code{NULL} for default prior.
#' @return A list containing (1) \code{$out}: a list as S3 object. The list contains the data.frame \code{$data} with the x and y coordinates and value of input data in mm; \code{$field}, a list containing the estimated tephra thickness field on a 500x500 grid in mm and the associated coordinates; and \code{$volume} the estimated total volume in km\eqn{^3}. (2) \code{$plotparams}, a list containing some parameters for plotting
tephra.estimate <- function(fn,coordsys,UTMtz=NaN, super=FALSE, prior.range=NULL, prior.sigma=NULL){
  library(INLA)
  library(raster)

  UTMtz <- ifelse(is.nan(UTMtz), NaN, as.integer(UTMtz))

  Data <- read.csv(file=fn, header=TRUE, sep=",")
  Xunit <- names(Data)[1]            #Get lat lon column names
  Yunit <- names(Data)[2]
  obsunit <- names(Data)[3]
  baddata <- is.na(Data[3])
  if(length(Data[is.na(Data) == TRUE]) > 0){   #Remove any rows with missing data
    Data <- Data[-baddata,]
  }
  zerodata <- which(Data[3] == 0, arr.ind = TRUE)   #Remove zeros
  if(length(zerodata) > 0){
    message("Removing 0 thicknesses from data set")
    Data <- Data[-zerodata,]
    }

  obsunit <- trimws(tolower(obsunit))
  if(obsunit == "mm"){
    obscale <- 10
  }else if(obsunit == "cm"){
    obscale <- 100
  }else if(obsunit == "m"){
    obscale <- 10000
  }else{
    stop("Not a valid measurement unit: Must be mm, cm or m")
  }
  obs<- as.vector(t(log(Data[3]*obscale)))

  Xin <- Data[Xunit]
  Yin <- Data[Yunit]

  if(coordsys == "latlon" & super == FALSE){ #Convert to UTM if latlon
    spll <- SpatialPoints(cbind(Xin,Yin), proj4string=CRS("+proj=longlat +datum=WGS84"))
    spcoord <- spTransform(spll, CRS(paste0("+proj=utm +zone=",UTMtz," ellps=WGS84")))
    coord <- cbind(spcoord@coords[,1],spcoord@coords[,2])
  }else{
    coord <- as.matrix(cbind(Xin,Yin))
  }

  if(super==TRUE){ #Project into 3D cartesian coordinates
    coord = inla.mesh.map(coord, projection = "longlat")
  }
  ms<-0
  for(ci in seq(dim(coord)[1])){
    dif <- apply(coord, 1, function(x){x-coord[ci,]})
    ms <- max(ms, sqrt(max(colSums(dif)**2)))
  }
  mesh = inla.mesh.create.helper(
    points=coord, offset=c(0.3*ms,2.5*ms), max.edge=c(0.1*ms,0.4*ms), min.angle=c(21,21), cutoff=1e-3)

  if(is.null(prior.range)){
    prior.range <- c(5*ms, 0.2) ### "Default": P(practic.range<5*ms)=0.2
  }else if(super==TRUE){
    prior.range[1] <- prior.range[1]/40075 * 2*pi   #Convert to projection on 3D cartesian
  }else if((coordsys == "UTM") | (coordsys == "latlon")){
    prior.range[1] <- prior.range[1]*1000 #Convert to m
  }

  if(is.null(prior.sigma)){
    prior.sigma <- c(1, 0.01) ### "Default": P(sigma>1)=0.01
  }
  spde <- inla.spde2.pcmatern(
    mesh=mesh, alpha=2, ### mesh and smoothness parameter
    prior.range=prior.range,
    prior.sigma=prior.sigma)

  Aop <- inla.spde.make.A(mesh, loc=coord)
  stack <- inla.stack(data=list(tephra=obs), A=list(Aop),
                      effects=list(i=1:spde$n.spde), tag='est')

  resi <- inla(tephra ~ -1 + f(i, model=spde),
               data=inla.stack.data(stack),
               control.predictor=list(A=inla.stack.A(stack)))#, compute=TRUE))

  pgrid0 <- inla.mesh.projector(mesh,  dims=c(500,500))
  prd0.m <- inla.mesh.project(pgrid0, resi$summary.ran$i$mode)
  prd0.mu <- inla.mesh.project(pgrid0, resi$summary.ran$i$mean)

  xUTM <- pgrid0$x
  yUTM <- pgrid0$y

  inds=which(exp(prd0.m) > 10, arr.ind = TRUE)
  if(super == TRUE){
    xlab = "Longitude"
    ylab = "Latitude"
    spgeo <- SpatialPoints(cbind(xUTM, yUTM), proj4string=CRS("+proj=longlat +datum=WGS84"))
    rast <- raster(spgeo,nrows=500, ncols=500)
    area <- area(rast)
  }else if(coordsys == "dfs"){
    xlab = "km"
    ylab = "km"
    xsize <- abs(xUTM[1] - xUTM[2])     # projected grid size in km
    ysize <- abs(yUTM[1] - yUTM[2])
    area <- matrix(rep(xsize*ysize, dim(prd0.m)[1]*dim(prd0.m)[2]), dim(prd0.m)[1], dim(prd0.m)[2])                 # Area of grid square in sq kms
  }else{
    xlab = "UTM"
    ylab = "UTM"
    sputm <- SpatialPoints(cbind(xUTM, yUTM), proj4string=CRS(paste0("+proj=utm +zone=",UTMtz," +datumWG84")) )
    spgeo <- spTransform(sputm, CRS("+proj=longlat +datum=WGS84"))
    rast <- raster(spgeo,nrows=500, ncols=500)
    area <- area(rast)
  }

  #Get rid of any really small values
  badpix <- which(prd0.m <= log(1.1))
  prd0.m[badpix] = NA
  prd0.mu[badpix] = NA
  goodpix <- which(is.finite(prd0.mu))
  scale <- 1e-7

  volume <- sum(as.matrix(area)[goodpix] * exp(prd0.mu[goodpix])*scale)   #Volume in square kms

  field <- list(x=xUTM, y=yUTM, mm=exp(prd0.mu)/10, units="mm")
  if(super == TRUE){
    data <- data.frame(x=Xin, y=Yin, mm=exp(obs)/10)
  }else{
    data <- data.frame(x=coord[,1], y=coord[,2], mm=exp(obs)/10)
  }
  out<-list(data=data, field=field, volume=volume)
  class(out) <- c("volume units=km^3", "thickness units=mm", paste0("coord system=",coordsys))
  plotparams <- list(xlim=c(xUTM[min(inds[,1])],xUTM[max(inds[,1])]),
                     ylim=c(yUTM[min(inds[,2])],yUTM[max(inds[,2])]),
                     xlab=xlab, ylab=ylab)
  return(list(out=out, plotparams=plotparams))
}

#'Plot an isopach map
#'
#'Plot an isopach map after estimating the thickness
#'
#' @param outlist The output from \code{tephra.estimate()}
tephra.plot <- function(outlist){
  ticks <- (c(1, 2,5, 10, 20, 50, 100, 200, 500, 1000, 2000))
  contour(outlist$out$field$x , outlist$out$field$y, log(outlist$out$field$mm*10),
          labels=ticks, levels=log(ticks*10),
          asp=1,
          xlim = outlist$plotparams$xlim,
          ylim = outlist$plotparams$ylim,
          xlab=outlist$plotparams$xlab, ylab=outlist$plotparams$ylab)
  points(outlist$out$data$x,outlist$out$data$y, pch=20)
  title(main="Median thickness (mm)")

}

#' Estimate isopachs and total volume
#'
#' Estimate the isopachs and total volume given tephra
#' samples and optionally plot the result. The estimated
#' tephra field is returned on a 500x500 rasterised grid
#' with the corresponding coordinates and an estimate of
#' the total volume and the input data.
#'
#' @param datafile String of the csv file containing the sampled tephra deposits. This must have columns containing the X coordinate, Y coordinate and thickness of the tephra deposit. There must be headers for all columns, with the sampled thickness having the units as the header: <mm>, <cm> or <m> (without the < >).
#' @param coord String of the coordinate system used. This must be either \code{"UTM"}, which is the UTM coordinates (where the UTM zone must also be given); \code{"dfs"}, which is distance from some origin in kilometres; or \code{"latlon"} for WG84 latitude and longitude, where the UTM zone must also be supplied as units are converted to UTM if \code{super=FALSE}.
#' @param UTMzone (optional) The UTM zone if using UTM or latlon coordinates (and \code{super==FALSE})
#' @param plot (optional, default=\code{TRUE}) whether to plot the output
#' @param super (optional, default=\code{FALSE}) Set to \code{TRUE} if this is a super eruption. The calculation is then done on a globe rather than a flat plain. Coordinates must be in latlon. (This only works for super eruptions, i.e. those spanning several degrees longitude/latitude)
#' @param prior.range (optional, default=\code{NULL}) Prior probability for the spatial range of the tephra field in km using a penalised complexity prior. This is a vector of length 2, \code{c(range0, Prange)}, where Pr(range < range0) = Prange. Setting Prange as \code{NA} fixes range to range0. Leave as \code{NULL} for the default prior, approximated from the measurement field.
#' @param prior.sigma (optional, default=\code{NULL}) Prior probability for the standard deviation of the log-tickness of the tephra field defined using a penalised complexity prior. This is a vector of length 2, \code{c(sigma0, Psigma)}, where Pr(sigma > sigma0) = Psigma. Setting Psigma as \code{NA} fixes sigma to sigma0. Leave as \code{NULL} for default prior.
#' @return List as S3 object. List contains: \code{$data}: a data.frame  with the x and y coordinates and value of input data in mm; \code{$field}: a list containing the estimated tephra thickness on a 500x500 grid in mm and the associated coordinates; and \code{$volume}: the estimated total volume in km\eqn{^3}.
#' @examples
#'  tephrafield <- tephra.map("CerroNegro.csv", "UTM", UTMzone=16)
tephra.map <- function(datafile, coord, UTMzone=NaN, plot=TRUE, super=FALSE,prior.range=NULL, prior.sigma=NULL){
  coord <- ifelse(coord=="utm", "UTM", coord)
  coord <- ifelse(coord=="lonlat", "latlon", coord) #Let them off it they write lonlat
  if(!file.exists(datafile)) stop("Data file does not exist")
  if(!((coord =="UTM") | (coord == "dfs") | (coord == "latlon"))) stop("Not a valid coordinate type: Must be 'UTM', 'dfs' or 'latlon'")
  if((coord == "UTM") & !(is.finite(UTMzone))) stop("UTM must have valid UTM zone")
  if((coord == "latlon") & !(is.finite(UTMzone)) & (super==FALSE)) stop("latlon must have valid UTM zone for conversion")
  if((super == TRUE) & (coord != "latlon")) stop("Super eruptions must be calculated with coordinate type 'latlon'")
  if(super == TRUE) message("super=TRUE: If the measurement points are not distributed over a large area this may freeze.")
  if(!is.null(prior.range)){
    if(length(prior.range)>2) stop("prior.range must be a vector of length 2 or NULL")
    if((prior.range[2]>1) & (!is.na(prior.range[2]))) stop("Prange must be less than 1")
    if((prior.range[2]<0) & (!is.na(prior.range[2]))) stop("Prange must be greater than 0")
  }
  if(!is.null(prior.sigma)){
    if(length(prior.sigma)>2) stop("prior.sigma must be a vector of length 2 or NULL")
    if((prior.sigma[2]>1) & (!is.na(prior.sigma[2]))) stop("Psigma must be less than 1")
    if((prior.sigma[2]<0) & (!is.na(prior.sigma[2]))) stop("Psigma must be greater than 0")
  }
  outlist <- tephra.estimate(datafile, coord, UTMzone, super=super, prior.range=prior.range, prior.sigma=prior.sigma)
  if (plot == TRUE)
  {
    tephra.plot(outlist)
  }
  out <- outlist$out
  return(out)
}
