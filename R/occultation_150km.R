library(raster)
library(sp)
library(bioRad)
library(rhdf5)
library(dplyr)
library(magrittr)
library(assertthat)

#################################################
# User functions
#################################################

radarpos=function(x){
  stopifnot(inherits(x,"list"))
  lat=sapply(vpl, function(x) x$attributes$where$lat)
  lon=sapply(vpl, function(x) x$attributes$where$lon)
  height=sapply(vpl, function(x) x$attributes$where$height)
  radar=sapply(vpl, function(x) x$radar)
  data.frame(radar=radar,lat=lat,lon=lon,height=height)
}

write.pvol=function(pvol,file,nodata=-9999,undetect=-9998){
  h5createFile(file)
  fid=H5Fopen(file)
  for (i in 1:length(pvol$scans)){
    h5createGroup(file,paste("dataset",i,sep=""))
    for (j in 1:length(pvol$scans[[i]]$params)){
      h5createGroup(file,paste("dataset",i,"/data",j,sep=""))
      data=pvol$scans[[i]]$params[[j]]
      data=replace(data,is.nan(data),undetect) # replace NaN before NA!
      data=replace(data,is.na(data),nodata)
      h5write(data,file,paste("dataset",i,"/data",j,"/data",sep=""))
      group=paste("dataset",i,"/data",j,"/what",sep="")
      h5createGroup(file,group)
      gid=H5Gopen(fid,group)
      h5writeAttribute(attributes(data)$param,gid,"quantity")
      h5writeAttribute(1,gid,"gain")
      h5writeAttribute(nodata,gid,"nodata")
      h5writeAttribute(0,gid,"offset")
      h5writeAttribute(undetect,gid,"undetect")
      H5Gclose(gid)
    }

    # write scan attributes
    attrgroupnames=names(pvol$scans[[i]]$attributes)
    for(k in 1:length(attrgroupnames)){
      group=paste("dataset",i,"/",attrgroupnames[k],sep="")
      h5createGroup(file,group)
      gid=H5Gopen(fid,group)
      attrgroup=pvol$scans[[i]]$attributes[[k]]
      attribnames=names(attrgroup)
      for (l in 1:length(attribnames)){
        h5writeAttribute(attrgroup[[l]],gid,attribnames[l])
      }
      H5Gclose(gid)
    }
  }

  # write volume attributes
  attrgroupnames=names(pvol$attributes)
  for(k in 1:length(attrgroupnames)){
    group=paste(attrgroupnames[k],sep="")
    h5createGroup(file,group)
    gid=H5Gopen(fid,group)
    attrgroup=pvol$attributes[[k]]
    attribnames=names(attrgroup)
    for (l in 1:length(attribnames)){
      h5writeAttribute(attrgroup[[l]],gid,attribnames[l])
    }
    H5Gclose(gid)
  }

  H5Fclose(fid)
}

occultationVolume=function(elevs,radinfo,dem,rscale=250,nazim=720,nrang=50){
  scans=list()
  for(elev in sort(elevs)){
    param=occultation(dem,radinfo$lat,radinfo$lon,radinfo$antenna,elev,rscale,nazim,nrang)
    scans=append(scans,list(params2scan(param)))
    # no need to calculate higher elevation angles when occultation-free
    if(length(which(c(param)!=0))==0) break
  }
  scans2pvol(scans,radinfo$radar,Sys.time())
}

scans2pvol=function(scans,radar,datetime,wavelength=10){
  if (FALSE %in% sapply(scans,is.scan)){
    stop("only scans expected as input")
  }
  vol=list()
  vol$radar=radar
  vol$datetime=datetime
  vol$scans=scans
  vol$attributes=list()
  vol$attributes$how$wavelength=wavelength
  vol$attributes$what$date=format(datetime, "%Y%m%d")
  vol$attributes$what$object="PVOL"
  vol$attributes$what$source=paste("RAD:",radar,sep="")
  vol$attributes$what$time=format(datetime, "%H%M%S")
  vol$attributes$where$height=scans[[1]]$geo$height
  vol$attributes$where$lat=scans[[1]]$geo$lat
  vol$attributes$where$lon=scans[[1]]$geo$lon
  vol$geo$height=scans[[1]]$geo$height
  vol$geo$lat=scans[[1]]$geo$lat
  vol$geo$lon=scans[[1]]$geo$lon
  class(vol)="pvol"
  vol
}

params2scan=function(...,startdate=Sys.time(),enddate=Sys.time()){
  if (FALSE %in% sapply(list(...),is.param)){
    stop("only scan parameters expected as input")
  }
  params=list(...)
  scan=list()
  scan$params=params
  param_names=sapply(params,function(x) attributes(x)$param)
  names(scan$params)=param_names
  scan$attributes=list()
  elangles=sapply(params,function(x) attributes(x)$geo$elangle)
  if(length(unique(elangles))!=1) stop("scan parameters with unequal elevations")
  rscales=sapply(params,function(x) attributes(x)$geo$rscale)
  if(length(unique(rscales))!=1) stop("scan parameters with unequal range dimension")
  ascales=sapply(params,function(x) attributes(x)$geo$ascale)
  if(length(unique(ascales))!=1) stop("scan parameters with unequal azimuth dimension")
  lats=sapply(params,function(x) attributes(x)$geo$lat)
  lons=sapply(params,function(x) attributes(x)$geo$lon)
  if(length(unique(lats))!=1 || length(unique(lons))!=1) stop("scan parameters have different radar latitude and/or longitude")
  heights=sapply(params,function(x) attributes(x)$geo$height)
  if(length(unique(ascales))!=1) stop("scan parameters with unequal antenna height")
  nbins=sapply(params,function(x) dim(x)[1])
  nazims=sapply(params,function(x) dim(x)[2])
  if(length(unique(nbins))!=1 || length(unique(nazims))!=1) stop("scan parameters have different dimensions")
  NIs=sapply(params,function(x) attributes(x)$how$NI)
  if(length(unique(NIs))!=1) stop("scan parameters with different nyquist intervals")

  scan$attributes=list()
  scan$attributes$how$NI=NIs[1]

  scan$attributes$where$elangle=elangles[1]
  scan$attributes$where$nbins=nbins[1]
  scan$attributes$where$nrays=nazims[1]
  scan$attributes$where$rscale=rscales[1]

  scan$attributes$what$enddate=format(enddate, "%Y%m%d")
  scan$attributes$what$endtime=format(enddate, "%H%M%S")
  scan$attributes$what$product="SCAN"
  scan$attributes$what$startdate=format(startdate, "%Y%m%d")
  scan$attributes$what$starttime=format(startdate, "%H%M%S")

  scan$geo=attributes(params[[1]])$geo
  class(scan)="scan"
  scan
}

getElev=function(elevmap,xy,crs){
  coordinates(xy)<-c("x","y")
  proj4string(xy)<-crs
  output=suppressWarnings(raster::extract(elevmap,xy))
  dimnames(output)=NULL
  output
}

dem2occultation=function(dem,elev,height){
  ranges=as.numeric(rownames(dem))
  beam_agl=height+beam_height(ranges,elev)
  beam_width_value=beam_width(ranges)
  occ=apply(dem, 2, function(x) 2*(beam_agl-x)/beam_width_value)
  occ[occ>1]=1     # no blockage
  occ[occ< -1]=-1  # full blockage
  occ=1-(occ+1)/2  # transform [-1,1] domain to [0,1]: 0 = no blockage, 1 = full blockage
  for (i in (1:(dim(occ)[2]))){
    maxocc=0
    for (j in (1:(dim(occ)[1]))){
      if(occ[j,i]>maxocc){
        maxocc=occ[j,i]
      }
      occ[j,i]=maxocc
    }
  }
  occ
}

occultation=function(elevmap,lat,lon,height,elev,rscale=1000,nazim=360,nrang=30){
  proj4string=CRS(paste("+proj=aeqd +lat_0=",lat," +lon_0=",lon," +ellps=WGS84 +datum=WGS84 +units=m +no_defs",sep=""))
  grid=lapply((0:(nazim-1))*360/nazim,function(azim) getElev(elevmap,polar2cartesian(range=seq(1,nrang)*rscale,azim=azim,elev),proj4string))
  DEM=do.call(cbind,grid)
  colnames(DEM)=(0:(nazim-1))*360/nazim
  rownames(DEM)=seq(1,nrang)*rscale
  data=dem2occultation(DEM,elev,height)
  geo=list()
  geo$lat=lat
  geo$lon=lon
  geo$height=height
  geo$elangle=elev
  geo$rscale=rscale
  geo$ascale=360/nazim
  geo$rstart=0
  attributes(data)$geo=geo
  attributes(data)$param="OCCULT"
  attributes(data)$how$NI=0
  class(data)=c("param",class(data))
  data
}

# wgs2proj is a wrapper for spTransform
# proj4string should be an object of class 'CRS', as defined in package sp.
# returns an object of class SpatialPoints
wgs2proj<-function(lon,lat,proj4string){
  xy <- data.frame(x = lon, y = lat)
  coordinates(xy) <- c("x", "y")
  proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")
  res <- spTransform(xy, proj4string)
  return(res)
}

# proj2wgs is a wrapper for spTransform
# proj4string should be an object of class 'CRS', as defined in package sp.
# returns an object of class SpatialPoints
proj2wgs<-function(x,y,proj4string){
  xy <- data.frame(lon=x, lat=y)
  coordinates(xy) <- c("lon", "lat")
  proj4string(xy) <- proj4string
  res <- spTransform(xy, CRS("+proj=longlat +datum=WGS84"))
  return(res)
}

polar2cartesian=function(range,azim,elev=0){
  y=range*sin(0.5*pi-azim*pi/180)*cos(elev*pi/180)
  x=range*cos(0.5*pi-azim*pi/180)*cos(elev*pi/180)
  data.frame(x=x,y=y)
}

summaryElev=function(elevmap,lat,lon,range=25000){
  radposxy=wgs2proj(lon,lat,proj4string(elevmap))@coords
  ext=raster::extent(radposxy[1]-range,radposxy[1]+range,radposxy[2]-range,radposxy[2]+range)
  extraction=raster::extract(elevmap,ext)
  data.frame(asl.min=min(extraction),asl.max=max(extraction),asl.mean=mean(extraction))
}

#############################################
# elevations to consider (from VCPs)
#############################################
# elevation angles used in the various VCPs
VCP11=c(0.5, 1.5, 2.4, 3.4, 4.3, 5.3, 6.2, 7.5, 8.7, 10, 12, 14, 16.7, 19.5)
VCP12=c(0.5, 0.9, 1.3, 1.8, 2.4, 3.1, 4.0, 5.1, 6.4, 8.0, 10.0, 12.5, 15.6, 19.5)
VCP21=c(0.5, 1.5, 2.4, 3.4, 4.3, 6.0, 9.9, 14.6, 19.5)
VCP31=c(0.5, 1.5, 2.5, 3.5, 4.5)
VCP215=c(0.5, 0.9, 1.3, 1,8, 2.4, 3.1, 4.0, 5.1, 6.4, 8.0, 10.0, 12, 14, 16.7, 19.5)
VCP_elev=sort(unique(c(VCP11,VCP12,VCP21,VCP31,VCP215)))
seq_elev = seq(0,3,0.1)

#VCP_elev=c(0.5,1.0,1.5,1.8,2.4,3.1,3.5,4.0,4.5,5.1,6.0,6.4,7.5,8.0,8.7,10.0,12.0,12.5,14.0,14.6,15.6,16.7,19.5)


#############################################
# get radar positions (radpos data.frame)
#############################################
#fs=list.files("~/Dropbox/radar/NEXRAD/coverage/pvol_ex_lower48",full.names=T,pattern="*.h5")
# vpl=read_vpfiles(fs)
# lead elevation map and
load(file="~/Dropbox/radar/NEXRAD/all/test_all_radars.RData")
vpl <- lapply(vpl,convert_legacy)
radpos=radarpos(vpl)
antenna=read.csv("~/Dropbox/radar/NEXRAD/occultation/NEXRAD_antenna_height.txt",stringsAsFactors = F)
radpos %>% left_join(antenna, by=c("radar"="site")) -> radpos
radpos$antenna=radpos$antenna*0.3048

#############################################
# making occultation maps
#############################################

# # lead digital elevation map
# dem=brick("~/Documents/gis/elev/elev48i0100a.tif")
# dem_alaska=brick("~/Documents/gis/elev/alaska/elevaki0100a.tif")
# dem_puerto_rico=brick("~/Documents/gis/elev/puerto_rico/EarthEnv-DEM90_N15W070.bil")
# dem_south_korea=brick("~/Documents/gis/elev/south_korea/EarthEnv-DEM90_N35E125.bil")
# dem_guam=merge(brick("~/gis/elev/guam/EarthEnv-DEM90_N10E140.bil"),brick("~/gis/elev/guam/EarthEnv-DEM90_N10E145.bil"))

# load radar positions
load_dem_data <- function(region){
  assert_that(region %in% c("lower48","alaska","puerto_rico","south_korea","guam","japan","hawaii"))
  if(region == "lower48"){
    return(brick("gis/elev/elev48i0100a.tif"))
  }
  if(region == "alaska"){
    return(brick("~/Documents/gis/elev/alaska/elevaki0100a.tif"))
  }
  if(region == "puerto_rico"){
    return(merge(brick("~/Documents/gis/elev/puerto_rico/EarthEnv-DEM90_N15W070.bil"),brick("~/Documents/gis/elev/puerto_rico/EarthEnv-DEM90_N20W070.bil"),brick("~/Documents/gis/elev/puerto_rico/EarthEnv-DEM90_N15W065.bil"),brick("~/Documents/gis/elev/puerto_rico/EarthEnv-DEM90_N20W065.bil")))
  }
  if(region == "south_korea"){
    return(merge(brick("~/Documents/gis/elev/south_korea/EarthEnv-DEM90_N35E125.bil"),brick("~/Documents/gis/elev/south_korea/EarthEnv-DEM90_N40E125.bil"),brick("~/Documents/gis/elev/south_korea/EarthEnv-DEM90_N35E120.bil"),brick("~/Documents/gis/elev/south_korea/EarthEnv-DEM90_N40E120.bil"),brick("~/Documents/gis/elev/south_korea/EarthEnv-DEM90_N30E120.bil"),brick("~/Documents/gis/elev/south_korea/EarthEnv-DEM90_N30E125.bil")))
  }
  if(region == "guam"){
    return(merge(brick("~/Documents/gis/elev/guam/EarthEnv-DEM90_N10E140.bil"),brick("~/Documents/gis/elev/guam/EarthEnv-DEM90_N10E145.bil")))
  }
  if(region == "japan"){
    return(merge(brick("~/Documents/gis/elev/japan/EarthEnv-DEM90_N20E125.bil"),brick("~/Documents/gis/elev/japan/EarthEnv-DEM90_N25E120.bil"),brick("~/Documents/gis/elev/japan/EarthEnv-DEM90_N30E120.bil"),brick("~/Documents/gis/elev/japan/EarthEnv-DEM90_N25E125.bil"),brick("~/Documents/gis/elev/japan/EarthEnv-DEM90_N30E125.bil"),brick("~/Documents/gis/elev/japan/EarthEnv-DEM90_N25E130.bil"),brick("~/Documents/gis/elev/japan/EarthEnv-DEM90_N30E130.bil")))
  }
  if(region == "hawaii"){
    return(merge(brick("~/Documents/gis/elev/hawaii/EarthEnv-DEM90_N15W155.bil"),brick("~/Documents/gis/elev/hawaii/EarthEnv-DEM90_N15W160.bil"),brick("~/Documents/gis/elev/hawaii/EarthEnv-DEM90_N15W165.bil"),brick("~/Documents/gis/elev/hawaii/EarthEnv-DEM90_N20W155.bil"),brick("~/Documents/gis/elev/hawaii/EarthEnv-DEM90_N20W160.bil"),brick("~/Documents/gis/elev/hawaii/EarthEnv-DEM90_N20W165.bil"),brick("~/Documents/gis/elev/hawaii/EarthEnv-DEM90_N25W155.bil"),brick("~/Documents/gis/elev/hawaii/EarthEnv-DEM90_N25W160.bil"),brick("~/Documents/gis/elev/hawaii/EarthEnv-DEM90_N25W165.bil")))
  }
}

load_antenna_data <- function(region){
  if(missing(region)){
    load(file="~/Dropbox/radar/NEXRAD/occultation/NEXRAD_antenna_height.RData")
    return(radpos)
  }
  assert_that(region %in% c("lower48","alaska","puerto_rico","south_korea","guam","japan","hawaii"))
  if(region == "alaska"){
    load(file="~/Dropbox/radar/NEXRAD/occultation/NEXRAD_ground_antenna_height_alaska.RData")
    return(radarInfo)
  }
  if(region == "puerto_rico"){
    load(file="~/Dropbox/radar/NEXRAD/occultation/NEXRAD_ground_antenna_height_puerto_rico.RData")
    return(radarInfo)
  }
  if(region == "south_korea"){
    load(file="~/Dropbox/radar/NEXRAD/occultation/NEXRAD_ground_antenna_height_korea.RData")
    return(radarInfo)
  }
  if(region == "guam"){
    load(file="~/Dropbox/radar/NEXRAD/occultation/NEXRAD_ground_antenna_height_guam.RData")
    return(radarInfo)
  }
  if(region == "japan"){
    load(file="~/Dropbox/radar/NEXRAD/occultation/NEXRAD_antenna_height.RData")
    return(radpos %>% dplyr::filter(radar=="RODN"))
  }
  if(region == "hawaii"){
    load(file="~/Dropbox/radar/NEXRAD/occultation/NEXRAD_antenna_height.RData")
    return(radpos %>% dplyr::filter(radar %in% c("PHKI","PHMO","PHKM","PHWA")))
  }
}

dem=load_dem_data("lower48")
radpos = NEXRAD_antenna_height

target_radars <- c("KRIW")
rscale = 250
## make occultation maps for all NEXRAD radars
## only works for lower48 radars because of dem extent
for (idxrad in 1:nrow(radpos)){
  current_radar <- radpos[idxrad, "radar"]
    if (current_radar %in% target_radars) {
    print(current_radar)
    occvol=occultationVolume(seq_elev,radpos[idxrad,],dem, nrang=1201)
    fname=paste("occult/h5-150km/KCYS/",occvol$radar,".h5",sep="")
    write.pvol(occvol,fname)
  }
}

## make occultation maps for all ex-lower48 radars:
for(region in c("alaska","puerto_rico","south_korea","guam","japan","hawaii")){
  dem=load_dem_data(region)
  radpos=load_antenna_data(region)
  for (idxrad in 1:nrow(radpos)){
   print(radpos[idxrad,"radar"])
   occvol=occultationVolume(VCP_elev,radpos[idxrad,],dem, nrang=150)
   fname=paste("~/Dropbox/radar/NEXRAD/occultation/h5-150km/",occvol$radar,".h5",sep="")
   write.pvol(occvol,fname)
  }
}

# plot occultation maps
library(ggplot2)
fsocc=list.files("~/Dropbox/radar/NEXRAD/occultation/h5-150km",pattern="*.h5$",full.names = T)
for (focc in fsocc){
  occmap=read_pvolfile(focc,param="OCCULT")
  pocc=plot(project_as_ppi(occmap$scans[[1]], range_max=150000),zlim=c(0,1))+ggtitle(basename(focc))

  ggsave(paste(gsub("h5-150km","jpeg-150km",focc),".jpeg",sep=""),pocc)
}



