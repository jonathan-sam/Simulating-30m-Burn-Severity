# using different algorithms and parameters for simulation modeling

library(raster)
library(sf)
library(doParallel)

source("globals.R")
source("functions.R")

timeOut <- system.time({
  
# set up parallel parameters (number of cores)
numCores <- getFromCommandLine('numCores',max(detectCores()/2 -1 ,1))  # leave one core free for other things on the machine
registerDoParallel(cores=numCores)
cat('Setting up to use',numCores,'cores.\n')


# set up which algorithm/parameter set we are using
experimentNum <- 14

modelName <- 'JS-SN_Mods2022-01-22_vegTopo_v1'
#modelName <- 'JS-NC_Mods2022-01-22_vegTopo_v1'

singleRandom <- FALSE
setFocalMatrix <- FALSE
smoothLastIteration <- FALSE
useHistoricalLMH <- TRUE
runStartIndex <- 1
runEndIndex <- 1
useSubset <- FALSE
fireSubset <- c()

firelistName <- "Historical-Forested-cbi-v2"
maxIterationNum <- 3
setExpFocalMatrix <- TRUE
filterRadius <- 2
setWindFocalMatrix <- TRUE
useMTBSBase <- FALSE
useCBIBase <- TRUE
sizeTolerancePercent <- 0.02
experimentName <- sprintf("MultiPass%02d-Smooth%02dwind-ExpFocal-Tol%02d-CBI-Clean",maxIterationNum, filterRadius*2+1, round(sizeTolerancePercent*100))

# get the fires we are interested in
load(paste0(currentDir(analysisRootDir, MTBSDirName), "historicalMTBSForestedFireList_cbi_v2.rda"))

# set up the model environment
modelEnv <- createModelEnv(modelName)

useCBSEMaskedInput <- exists('useCBSEMaskedInput',envir=modelEnv) && modelEnv$useCBSEMaskedInput
# set up directories (input and output)
if (useCBSEMaskedInput) {
  baseModelInputDir <- paste0(currentDir(analysisRootDir,cbseDirName),maskingTA30mDirectory)
} else {
  baseModelInputDir <- paste0(currentDir(analysisRootDir,MTBSDirName),MTBSMaskingDirectory)
}
baseFireModelDir <- paste0(currentDir(analysisRootDir,sim30mDirName),firelistName,'/',modelName,'/')
baseExperimentOutputDir <- paste0(baseFireModelDir,experimentName,'/')
baseRandomInputDir <- paste0(currentDir(analysisRootDir,sim30mDirName),firelistName,'/Inputs/')
if (!dir.exists(baseExperimentOutputDir)) {
  dir.create(baseExperimentOutputDir, recursive=T)
}

dnbrReclassDF <- data.frame(id=0:6,v=c(0,0,1,2,3,0,0))

focalMatrix <- NULL
jitterFrac <- NULL

stopifnot(exists("filterRadius"))
pixelSize <- MTBSGridSize
filterBase <- createInvExpDistFilter(2*filterRadius+1, pixelSize)
circleFilter <- focalWeight(raster(ncols=30,nrows=30,crs=MTBSMaskingCRS,resolution=pixelSize), d=filterRadius*pixelSize+1, type='circle')
invDistMatrix <- filterBase*circleFilter#/circleFilter[6,6]
invDistMatrix <- invDistMatrix / sum(invDistMatrix) # normalize
focalMatrix <- invDistMatrix
save(focalMatrix, file=paste0(baseExperimentOutputDir,"focalMatrix.rda"))
gc()

# set up smoothing matrix as appropriate for experiment
# other algorithm parameters?
# save out full set of parameters

# for every fire in the fire list
fireIDList <- fireList$MTBSFireID

fireIDOut <- foreach(fireIndex = seq_along(fireIDList),.combine='c') %:% 
  when(!useSubset || (fireIDList[fireIndex] %in% fireSubset)) %:%
  foreach(runIndex = seq(from=runStartIndex,to=runEndIndex,by=1), .combine='c') %dopar%
  {
  currentFire <- fireList[fireIndex,]
  currentFireID <- tolower(fireIDList[fireIndex])
  print(sprintf("%s-%02d",currentFireID, runIndex))

  modelInputDir <- paste0(baseModelInputDir,currentFireID,'/')
  randomInputDir <- paste0(baseRandomInputDir, currentFireID,'/')
  fireOutputDir <- paste0(baseExperimentOutputDir, currentFireID,'/')
  if (runIndex > 1) {
    randomInputDir <- sprintf("%s%s-%02d/",baseRandomInputDir, currentFireID, runIndex)
    fireOutputDir <- sprintf("%s%s-%02d/",baseExperimentOutputDir, currentFireID,runIndex)
  }
  if (!dir.exists(randomInputDir)) {
    dir.create(randomInputDir, recursive = T)
  }
  if (!dir.exists(fireOutputDir)) {
    dir.create(fireOutputDir,recursive=T)
  }
  
  boundaryPoly <- read_sf(paste0(modelInputDir,'boundary.shp')) %>%
    st_transform(MTBSMaskingCRS)
  # the sample data is in Teale Albers in case it matters for x/y
  # and for some reason, the boundary poly is not
  #boundaryPoly <- st_transform(boundaryPoly, MTBSMaskingCRS) %>% st_geometry

  # random rasters for each iteration and model, save them to the 
  # get initial prediction matrices: createProbMapFromDir
  # this is the 0th iteration
  iterationNum <- 0
  
  # for further iterations (1-whatever): createProbMapFromPrev
  # iterationNum <- 1
  for(iterationNum in seq(from=0, to=maxIterationNum)) {
    print(iterationNum)

    fireIterOutputDir <- sprintf("%sIteration_%02d/",fireOutputDir, iterationNum)
    if (!dir.exists(fireIterOutputDir)) {
      dir.create(fireIterOutputDir, recursive=T)
    }
    randomIterInputDir <- sprintf("%sIteration_%02d/",randomInputDir, iterationNum)
    if (singleRandom) {
      randomIterInputDir <- sprintf("%sIteration_%02d/",randomInputDir, 0)
    }
    if (!dir.exists(randomIterInputDir)) {
      dir.create(randomIterInputDir, recursive=T)
    }
    if (iterationNum == 0) {
      probabilityMaps <- modelEnv$createProbMapFromDir(extent(st_buffer(boundaryPoly,1000)),
                                                       modelInputDir,
                                                       fireIterOutputDir)
      if (setWindFocalMatrix) {
        if (useCBSEMaskedInput) {
          wdir_mean <- subset(modelEnv$savedInputStack,"winddirearth_monthly.mean")
        } else {
          wdir_mean <- subset(modelEnv$savedInputStack,paste("wdir","monthly","mean",sep="_"))
        }
        # }
        windAngle <- mean(wdir_mean[], na.rm=T)
        windDistFilter <- modelEnv$createWindFilter(filterRadius*2+1, windAngle) * invDistMatrix
        focalMatrix <- windDistFilter/sum(windDistFilter)
        save(focalMatrix, file=paste0(fireIterOutputDir,"focalMatrix.rda"))
        
      }
    } else {
      probabilityMaps <- modelEnv$createProbMapFromPrevious(resStack,
                                                            modelInputDir,
                                                        fireIterOutputDir)
    }
    
    if (is.null(probabilityMaps)) {
      print("FAILURE (probabilityMaps)")
    } else {
      firePerimeter <- subset(probabilityMaps,1)
      firePerimeter[] <- 1
      firePerimeter <- mask(firePerimeter, boundaryPoly)
      
      fireFractions <- modelEnv$createFireFractions(currentFire, useHistoricalLMH)
      
      smoothFocal <- focalMatrix
      if (smoothLastIteration && (iterationNum != maxIterationNum)) {
        smoothFocal <- NULL
      }
      resStack <<- createFireSeverityMapsExp(firePerimeter, probabilityMaps, fireFractions, smoothFocal, dependentLayers = modelEnv$layersDependOnPrev, outputDir = fireIterOutputDir, randomInputDir=randomIterInputDir, sizeTolPer = sizeTolerancePercent)
      rm(probabilityMaps); gc()
      if (is.null(resStack)) {
        print("FAILURE (severity stack)")
      } else {
        severityRast <- modelEnv$createSeverityRast(resStack) %>% mask(firePerimeter)
        
        severityRast <- writeRaster(severityRast, paste0(fireIterOutputDir,'Severity'),
                                    overwrite=T, format='GTiff',
                                    options=c("TILED=YES","COPY_SRC_OVERVIEWS=YES","COMPRESS=LZW"))
        
        if (useMTBSBase) {
          if (iterationNum == 0) {
            observedRast <- dnbrRast %>%
              projectRaster(severityRast, method='ngb') %>%
              crop(severityRast)
            observedRast <- writeRaster(observedRast, filename = paste0(fireOutputDir,'ObservedImage.tif'), overwrite=T,
                                        format='GTiff', options=c("TILED=YES","COPY_SRC_OVERVIEWS=YES","COMPRESS=LZW"))
          }
          diffImage <- observedRast - severityRast
          writeRaster(diffImage, filename = paste0(fireIterOutputDir,'DiffImage.tif'), overwrite=T,
                      format='GTiff', options=c("TILED=YES","COPY_SRC_OVERVIEWS=YES","COMPRESS=LZW"))
          
          rm(diffImage)
          gc()
        } else if (useCBIBase) {
          if (iterationNum == 0) {
            cbiFilename <- paste0(modelInputDir,'cbi_severity.tif')
            observedRast <- NULL
            if (file.exists(cbiFilename)) {
              observedRast <- raster(cbiFilename) %>% reclassify(dnbrReclassDF) %>%
                projectRaster(severityRast, method='ngb') %>%
                crop(severityRast)
              observedRast <- writeRaster(observedRast, filename = paste0(fireOutputDir,'ObservedImage.tif'), overwrite=T,
                                          format='GTiff', options=c("TILED=YES","COPY_SRC_OVERVIEWS=YES","COMPRESS=LZW"))
              
            }
            if (exists('observedRast') && !is.null(observedRast)) {
              diffImage <- observedRast - severityRast
              writeRaster(diffImage, filename = paste0(fireIterOutputDir,'DiffImage.tif'), overwrite=T,
                          format='GTiff', options=c("TILED=YES","COPY_SRC_OVERVIEWS=YES","COMPRESS=LZW"))
              
              rm(diffImage)
              gc()
            }
            
          }
        }
        
      }
    }
  }
  
  }
}) #end of system.time
print(timeOut)
