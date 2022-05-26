# set up prediction model, must include support functions. Specifically:
#  createProbMapFromDir (boundaryPolygon, dataDir)
#  createProbMap (boundaryPolygon)

# uses an iterative model

# set up of prediction model from Jonathan Sam for Sierra Nevada (2021-09-22)

# pulled together by Jonathan Baldwin

library(raster)
library(mgcv)

#source('globals.R')

load(paste0(workspaceRootDir,'Baldwin/Models/JonathanSam/JS-Mods_2022-01-22/Downscaling_SNMods0122.vegTopo.v1.rda'))
LMH.mod.i0 <- LMHmod_SN
MH.mod.i0 <- MHmod_SN
H.mod.i0 <- Hmod_SN
rm(LMHmod_SN, MHmod_SN, Hmod_SN)

LMH.mod.i1 <- LMHmod_SN_wW_all
MH.mod.i1 <- MHmod_SN_wW_all
H.mod.i1 <- Hmod_SN_wW_all
rm(MHmod_SN_wW_all, Hmod_SN_wW_all, LMHmod_SN_wW_all)

requiredFields.i0 <- unique(c(attr(LMH.mod.i0$terms,'term.labels'),
                       attr(MH.mod.i0$terms,'term.labels'),
                       attr(H.mod.i0$terms,'term.labels')))

requiredFields.i1 <- unique(c(attr(LMH.mod.i1$terms,'term.labels'),
                              attr(MH.mod.i1$terms,'term.labels'),
                              attr(H.mod.i1$terms,'term.labels')))

requiredFields <- unique(c(requiredFields.i0, requiredFields.i1))

layersDependOnPrev <- T
layerNames <- c('LMH','MH','H')

# creates the wind filter based on the angle (clockwise) from north
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
  resfilter <- ifelse(!is.na(resfilter),resfilter,0)
  return(resfilter)    
}

# create the fractions for each severity class
createFireFractions <- function(currentFire, useOldLMH=TRUE) {
  if(useOldLMH) {
    highFrac <- currentFire$high
    modFrac <- currentFire$mod
    lowFrac <- currentFire$low
  } else {
    highFrac <- currentFire$pH
    modFrac <- currentFire$pM * (1.0 - highFrac)
    lowFrac <- currentFire$pL * (1.0 - highFrac - modFrac)
  }
  return(c(LMH=lowFrac + modFrac + highFrac,
                     MH=modFrac + highFrac,
                     H=highFrac))
}

# create the raster of severity: 1-low, 2-mod, 3-high
createSeverityRast <- function(severityStack) {
  return(sum(resStack,na.rm=T))
}

# creates a specific severity class prediction layer from the input variables and the model
# returns the stack of [0,1) probabilities
createPredLayer <- function(inputStack, modelObj, modelName, outputDir=NULL) {
  outputLayer <- predict(inputStack, modelObj, type='response')
  outputLayer[is.na(outputLayer[])] <- 0
  if (!is.null(outputDir)) {
    writeRaster(outputLayer, paste0(outputDir, 'Predict_', modelName),
                format='GTiff', overwrite=T,
                options=c("TILED=YES","COPY_SRC_OVERVIEWS=YES","COMPRESS=LZW"))
  }
  return(outputLayer)
}

savedInputStack <- stack()
savedInputDir <- NULL
savedTreatmentValue <- NA

# create the probability map from a given directory. Calls above for each layer
# mapExtents: Extent object specifying the extents of cropping
# dataDir: directory with input (previously clipped data), all accululated (biomass, vegetation, etc.)
# outputDir: (optional) directory to write probability layers or other temp files
createProbMapFromDir <- function(mapExtents, dataDir, outputDir=NULL) {

  # savedInputStack <- crop(savedInputStack, mapExtents)
  savedInputStack <<- getInputStack(mapExtents, dataDir)
  
  probStack <- stack(createPredLayer(savedInputStack, LMH.mod.i0, 'LMH', outputDir))
  probStack <- addLayer(probStack,createPredLayer(savedInputStack, MH.mod.i0, 'MH', outputDir))
  probStack <- addLayer(probStack,createPredLayer(savedInputStack, H.mod.i0, 'H', outputDir))
  names(probStack) <- c('LMH','MH','H')

  return(probStack)
}

# pull out an input stack
getInputStack <- function(mapExtents, dataDir) {

  inputStack <- stack()
  
  for(fieldName in requiredFields) {
    origFieldName <- fieldName
    
    # various calculatd or otherwise extracted values
    if (substr(fieldName,1,5)=="dist_") { 
      # calculatd distance from other fild
      cat("File",fieldName,"assuming calculated and grabbing input\n")
      fieldParts <- unlist(strsplit(fieldName,"_"))
      if (fieldParts[1] =="dist") {
        timeframe <- fieldParts[3]
        varname <- fieldParts[4]
        fieldName <- sprintf("%s_%s_mean",varname,timeframe)
      } else {
        cat("cannot decode field. Skipping\n")
        next
      }
    }
    newLayer <- raster(paste0(dataDir, fieldName,'.tif'))
    inputStack <- addLayer(inputStack, newLayer)
  }

  inputStack <- crop(inputStack, mapExtents)

  return(inputStack)
}

# creates a probability map (stack) based on a previous burnd stack, used for subsequent iterations
createProbMapFromPrevious <- function(prevBurnedStack,
                                      dataDir,
                                      outputDir=NULL) {
  
  # get an inverse distance filter to apply with wind distance
  invDistFilter11x11 <- createInvDistFilter(11, MTBSGridSize) * focalWeight(prevBurnedStack, 151, type='circle')
  invDistFilter11x11 <- invDistFilter11x11 / sum(invDistFilter11x11)
  invDistFilter21x21 <- createInvDistFilter(21, MTBSGridSize) * focalWeight(prevBurnedStack, 301, type='circle')
  invDistFilter21x21 <- invDistFilter21x21 / sum(invDistFilter21x21)

  # assume fields that only exists in subsequent iterations are calculated, just process those  
  calculatedFieldsList <- requiredFields.i1[!(requiredFields.i1 %in% requiredFields.i0)]
  for (calcField in calculatedFieldsList) {
    fieldParts <- unlist(strsplit(calcField,"_"))
    
    if (fieldParts[1] =="dist") {
      # the distance based fields
      inRastLayerName <- fieldParts[2]
      if (substr(inRastLayerName,1,3) == "cbi") {
        inRastLayerName <- substr(inRastLayerName,4,nchar(inRastLayerName))
      }
      inRast <- subset(prevBurnedStack,inRastLayerName)

      timeframe <- fieldParts[3]
      varname <- fieldParts[4]

      wdir_mean <- subset(savedInputStack,paste(varname,timeframe,"mean",sep="_"))
      windAngle <- mean(wdir_mean[], na.rm=T)
      
      filterSize <- fieldParts[5]
      if (filterSize == "f11") {
        windDistFilter11 <- createWindFilter(11, windAngle) * invDistFilter11x11
        windDistFilter11 <- windDistFilter11/sum(windDistFilter11)

        rastLayer <- focal(inRast, windDistFilter11, pad=T, padValue=0, na.rm=T)
        names(rastLayer) <- calcField
        savedInputStack <- addLayer(savedInputStack, rastLayer)

      } else {
        cat("filter",filterSize,"not implements.aborting\n")
        stop()
      }
    } else {
      cat("calcField",calcField,"not implemented.skipping\n")
    }
  }

  probStack <- stack(createPredLayer(savedInputStack, LMH.mod.i1, 'LMH', outputDir))
  probStack <- addLayer(probStack,createPredLayer(savedInputStack, MH.mod.i1, 'MH', outputDir))
  probStack <- addLayer(probStack,createPredLayer(savedInputStack, H.mod.i1, 'H', outputDir))
  names(probStack) <- c('LMH','MH','H')

  return(probStack)
  
}

# store a treatment value in the environment for future use
setTreatmentValue <- function(treatmentValue) {
  savedTreatmentValue <<- treatmentValue
}

# create a probability map stub
createProbMap <- function(mapExtents, outputDir=NULL) {
  cat("createProbMap not implemented\n")
  return(NULL)
}