# global variables and functions to generate strings
library(raster) # needed for crs() function

# these next few lines are for file system specific variables
systemRootDir <- '/data/'

# should not need to modify below for file system changes, except for systemic changes
# that are affecting all systems.

MTBSDirName <- 'MTBS'
BoundariesDirName <- 'Boundaries'
sim30mDirName <- 'Simulation-30m'
cbseDirName <- 'CBSE'

dlRootDir <- paste0(systemRootDir, 'Downloads/')
analysisRootDir <- paste0(systemRootDir, 'Analysis/')
workspaceRootDir <- paste0(systemRootDir, 'Workspaces/')
backupRootDir <- paste0(systemRootDir, 'Backups/')

# PLEASE CHANGE when appropriate (and keep yyyy-mm-dd format for sorting)
sim30mExperimentDate <- '2021-02-25'

# can be used to verify projections. Extracted from MTBS data: MTBS.2020-04-20
crsMTBS <- crs("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0") 

# CRS for use with lan/long projectsions (like gridMET). verified against gridMET.2020-04-20 (metdata_elevation.nc)
crsGridMET <- crs("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

# CRS from CA State Boundary provided by Pyregence
crsTealeAlbers <- crs("+proj=aea +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 +y_0=-4000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs ")

# 30m was too tight for a single make grid so there isn't a Teale Albers version
gridSizeMetersList <- c(6000,4000,3000,2000,1500,1000,300,150,100, 30)

currentRootFmt <- "%sCurrent/"
archiveRootFmt <- "%sArchive/"

# directory construction functions to pull generalized directory names
# current is a link to the version of a data set that is currently being used generally
currentDir <- function(rootDir, dataName) {
  return(sprintf(paste0(currentRootFmt,"%s/"), rootDir, dataName))
}

# archive data sets by version, in case we want to pull a non-current version
archiveDir <- function(rootDir, dataName, archiveSuffix) {
  return(sprintf(paste0(archiveRootFmt,"%s.%s/"), rootDir, dataName, archiveSuffix))
}

MTBSYearRange <- 1984:2019

# used in MTBS Sampling
MTBSMaskingDirectory <- 'MaskedFires_crsTealeAlbers/'
MTBSMaskingCRS <- crsTealeAlbers
MTBSGridSize <- 30
MTBSMaskingGrid <- sprintf("%sCalifornia/Grids/TealeAlbers/grid-%04d-meters.tif ",currentDir(analysisRootDir,BoundariesDirName),MTBSGridSize)
MTBSSamplingVersion <- '2020_09_29_crsTealeAlbers'
MTBSSamplingSizePercent <- 10

# general use variables
maskingGridTA30mName <- sprintf("%sCalifornia/Grids/TealeAlbers/grid-%04d-meters.tif ",currentDir(analysisRootDir,BoundariesDirName), 30)
maskingTA30mDirectory <- 'MaskedFires_crsTealeAlbers30m/'
sampling30mSizePercent <- MTBSSamplingSizePercent
sampling30mVersion <- '2022_02_03_crsTealeAlbers'

# create a matrix (square) filter based on 30m pixel distances
# center is 0,0 offset and goes up/down from there
createInvDistFilter <- function(filterSize, pixelSize=30) {
  stopifnot((filterSize%%2)==1)
  offsetSize <- floor(filterSize/2)
  offsetDist <- offsetSize+1
  resfilter <- 
    sapply(1:filterSize,function(i) {
      sapply(1:filterSize, function(j) {
        1/(sqrt((i-offsetDist)**2 + (j-offsetDist)**2)*pixelSize) # assume distance in 30m increments for testing WJB
      })
    })
  resfilter <- ifelse(!is.infinite(resfilter),resfilter,0)
  return(resfilter)    
}

# inverse square distance filter
createInvSquareFilter <- function(filterSize, pixelSize=30) {
  stopifnot((filterSize%%2)==1)
  offsetSize <- floor(filterSize/2)
  offsetDist <- offsetSize+1
  resfilter <- 
    sapply(1:filterSize,function(i) {
      sapply(1:filterSize, function(j) {
        1/((sqrt((i-offsetDist)**2 + (j-offsetDist)**2)**2)*pixelSize) # assume distance in 30m increments for testing WJB
      })
    })
  resfilter <- ifelse(!is.infinite(resfilter),resfilter,0)
  return(resfilter)    
}

# inverse exponential distance(ish) filter (inverse exponential pixel index distance combined with inverse pixel size)
createInvExpDistFilter <- function(filterSize, pixelSize=30) {
  stopifnot((filterSize%%2)==1)
  offsetSize <- floor(filterSize/2)
  offsetDist <- offsetSize+1
  resfilter <- 
    sapply(1:filterSize,function(i) {
      sapply(1:filterSize, function(j) {
        1/(exp(sqrt((i-offsetDist)**2 + (j-offsetDist)**2))*pixelSize) # assume distance in 30m increments for testing WJB
      })
    })
  resfilter <- ifelse(!is.infinite(resfilter),resfilter,0)
  return(resfilter)    
}

#create wind filter so we have a single coherent function for all uses, 
# taking from Jonathan's model 2021-09-22-v2
# filterSize: size of matrix (square), assumed to be odd (e.g. 11 is most likely)
# WindAngle: the angle (in degrees)  of the wind (clockwise from North? met wind direction?)
# pixelSize: the size (in meters) of a pixel, default 30m
createWindFilter <- function(filterSize, WindAngle, pixelSize=30) {
  stopifnot((filterSize%%2)==1)
  offsetSize <- floor(filterSize/2)
  offsetDist <- offsetSize+1
  windCoords <- c(sin(-WindAngle*pi/180), cos(-WindAngle*pi/180)) # switch normal x/y intentionally
  resfilter <- 
    sapply(1:filterSize,function(i) {
      sapply(1:filterSize, function(j) {
        pixelCoords <- c((i-offsetDist), (j-offsetDist))
        pixelLength  <- sqrt(sum(pixelCoords**2))
        pixelCoords <- pixelCoords / pixelLength
        a <- acos(pixelCoords[1] * windCoords[1] + pixelCoords[2] * windCoords[2])
        return(a)
      })
    })
  resfilter <- ifelse(!is.na(resfilter),resfilter,pi/2.0)
  return(resfilter)    
}

# command line processing
getFromCommandLine <- function(varName, defaultValue) {
  if (!exists('commandLineArgs')) {
    commandLineArgs <- strsplit(commandArgs(trailingOnly = TRUE), split='=')
  }
  for(argIndex in seq_along(commandLineArgs)) {
    newArg <- commandLineArgs[[argIndex]]
    if (newArg[1] == varName) {
      newVal <- TRUE
      if (length(newArg)>1) {
        newVal <- newArg[2]
      }
      return(newVal)
    }
  }
  return(defaultValue)
}

# set raster dir temp files
rasterTempDir <- paste0(analysisRootDir, 'tmp/')
if (!dir.exists(rasterTempDir)) {
  dir.create(rasterTempDir, recursive=T)
}
rasterOptions(tmpdir=rasterTempDir)

