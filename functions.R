# functions used in 30m simulation
# should be a collection of functions but no other code (other than library initializations)

library(raster)
library(sf)

#source("globals.R")

# this function is the heart of the algorithm - it generates a smoothed/grouped set of
# burned pixels
# it takes a raster (one layer) of random values (uniform is used above), a probability raster,
# a mask of where to calculate (the fire perimeter mask possibly combined with previous levels), 
# a filter matrix for the smoothing, and an estimate of the size of the burned area (in pixels)
# the smoothJitterRast is there to apply a jitter to the final smoothed burn raster, expected to be the same size as the others
# there is also the output from the previous iteration, a tolerance for the size of the generated
# severity class, a directory and suffix for output
# This function saves out most of the intermediate rasters into the output directory if it is included
calcIsBurnedExp <- function(probRast, 
                            randRast, 
                            maskRast, 
                            estSize, 
                            sizeTolerance=(estSize*0.05),
                            filterMat = NULL, 
                            prevBurned = NULL,
                            smoothJitterRast=NULL,
                            sizeTolPer=0.05, 
                            outputDir=NULL,
                            outputSuffix=NULL) {
  applySmoothing <- !is.null(filterMat)
  
  sizeTolerance <- max(estSize*sizeTolPer,1)
  
  smoothBurned <- (randRast <= probRast)
  
  # write out initial raster  
  if (!is.null(outputDir)) {
    writeRaster(smoothBurned, sprintf("%spreConv_%s", outputDir, outputSuffix),
                overwrite=T, format='GTiff',
                options=c("TILED=YES","COPY_SRC_OVERVIEWS=YES","COMPRESS=LZW"))
  }
  
  
  if (applySmoothing) {
    if (!is.null(prevBurned)) {
      focalInputRast <- (smoothBurned | prevBurned)
    } else {
      focalInputRast <- smoothBurned
    }
    smoothBurned <- focal(focalInputRast, filterMat)
    # write out smoothed raster
    if (!is.null(outputDir)) {
      writeRaster(smoothBurned, sprintf("%scolvolved_%s", outputDir, outputSuffix),
                  overwrite=T, format='GTiff',
                  options=c("TILED=YES","COPY_SRC_OVERVIEWS=YES","COMPRESS=LZW"))
    }
  }
  
  if (!is.null(smoothJitterRast)) {
    smoothBurned <- clamp(smoothBurn + smoothJitterRast , lower=0, upper=1, useValues=T)
  }
  
  smoothBurned <- mask(smoothBurned, maskRast)
  isBurned <- smoothBurned > 0
  
  if (applySmoothing) {
    sumBurned <- sum(isBurned[]>0, na.rm=T)
    cat('size', estSize,'(', estSize*(1-sizeTolPer), estSize*(1+sizeTolPer),') from ',sumBurned, '\n')
    if (estSize < 1) {
      cat("Size < 1. Returning unburned\n")
      isBurned[!is.na(isBurned[])] <- 0
    } else if (sumBurned < estSize) {
      cat('Not enough pixels for burning, returning what is there.\n')
    } else {
      rangeMax <- max(smoothBurned[], na.rm=T)
      rangeMin <- min(smoothBurned[], na.rm=T)
      cat('Iterating to generate burn map',rangeMax,rangeMin,'\n')
      # iterate while not within +/- sizeTolerance (default 5%)
      maxIterations <- 200
      iterCount <- 0 # need to cap the maximum number of iterations so we don't end up in an infinite loop
      while((abs(sumBurned-estSize)>sizeTolerance) && (iterCount < maxIterations)) {
        iterCount <- iterCount + 1
        if (iterCount >= maxIterations) {
          cat('max iteration reached, last step possible\n')
        }
        isBurned <- smoothBurned >= (rangeMax + rangeMin)/2.0
        sumBurned <- sum(isBurned[], na.rm=T)
        cat(rangeMax, rangeMin, sumBurned,'\n')
        if (sumBurned > (estSize+sizeTolerance)) {
          # too big more than sizeTolerance above estimate, set new lower bound
          rangeMin <- (rangeMax + rangeMin)/2
        } else if (sumBurned < (estSize-sizeTolerance)) {
          # too small, more than sizsTolerance below estimate, set new upper bound
          rangeMax <- (rangeMax + rangeMin)/2.0
        } else {
          #  within estimated size +/- sizeTolerance
          break
        }
      }
    }
    
  }
  
  return(isBurned)
}

# creates a set of fire severity maps from a given boundary and probability maps (as a brick or stack)
# also need the fractions array, names of layers and array should be the same and dependent
# on each other (layer 1 first, layer 2 from layer 1, etc.)
# can also have a single probability map for all fractions
# provides an external directory to save intermediate data files for repeatability
# fireMask to mask out the pixels of interest

createFireSeverityMapsExp <- function(fireMask, 
                                      probabilityMaps, 
                                      fractionArray, 
                                      focalMatrix = NULL, 
                                      dependentLayers = TRUE, 
                                      outputDir=NULL, 
                                      jitterSmooth=NULL, 
                                      randomInputDir=NULL,
                                      sizeTolPer=0.05) {
  writeTempFiles <- !is.null(outputDir)
  readRandomFiles <- !is.null(randomInputDir)
  if(readRandomFiles && !writeTempFiles) {
    cat("Cannot read random files, temp dir not specified. Ignoring\n")
    readRandomFiles <- F
  }

  numFractions <- length(fractionArray)
  
  # check valid number (matching) of layers for probabilities with fractions
  numProbLayers <- nlayers(probabilityMaps)
  if (!((numProbLayers == 1) || (numProbLayers == numFractions))) {
    cat('Number of probability layers is not equal to the number of fractions (or 1). Aborting\n')
    return(NULL)
  }
  
  # check valid bounds of probabilty maps and boundary. they have to match here
  if (extent(fireMask) != extent(probabilityMaps)) {
    cat("Extents of boundary and probability maps don't match. Aborting\n")
    return(NULL)
  }
  
  prevBurned <- fireMask
  prevBurned[] <- 0
  
  totalSize <- sum(!is.na(fireMask[]))
  
  returnStack <- stack()
  
  layerNames <- names(fractionArray)
  if (is.null(layerNames)) {
    layerNames <- sprintf("layer%02d",seq_along(fractionArray))
  }
  
  # loop over all fractions
  for(fractionIdx in seq_along(fractionArray)) {
    curFraction <- fractionArray[fractionIdx]
    curProb <- subset(probabilityMaps,min(numProbLayers, fractionIdx))
    
    randomFileName <- sprintf("%srand_%s",randomInputDir, layerNames[fractionIdx])
    generateRandom <- !readRandomFiles
    if (readRandomFiles) {
      if (file.exists(paste0(randomFileName,".tif"))) {
        curRand <- raster(paste0(randomFileName,".tif"))
        curRand <- crop(curRand, curProb)
      } else {
        generateRandom <- TRUE
      }
    }
    if (generateRandom) {
      curRand <- curProb
      curRand[] <- runif(ncell(curRand))
      if (readRandomFiles) {
        curRand <- writeRaster(curRand, randomFileName,
                               overwrite=T, format='GTiff',
                               options=c("TILED=YES","COPY_SRC_OVERVIEWS=YES","COMPRESS=LZW"))
      }
      if (writeTempFiles) {
        curRand <- writeRaster(curRand, sprintf("%srand_%s",outputDir, layerNames[fractionIdx]),
                               overwrite=T, format='GTiff',
                               options=c("TILED=YES","COPY_SRC_OVERVIEWS=YES","COMPRESS=LZW"))
      }
    }
    curJitter <- NULL
    if (!is.null(jitterSmooth) ) {
      curJitter <- curRand
      curJitter[] <- runif(ncell(curJitter),min=-1*jitterSmooth/2, max=jitterSmooth/2)
      if (writeTempFiles) {
        writeRaster(curJitter, sprintf("%ssmoothJitter_%s", outputDir, layerNames[fractionIdx]),
                    overwrite=T, format='GTiff',
                    options=c("TILED=YES","COPY_SRC_OVERVIEWS=YES","COMPRESS=LZW"))
      }
    }
    
    isBurned <- calcIsBurnedExp(curProb, curRand, fireMask, filterMat=focalMatrix, estSize=totalSize*curFraction, smoothJitterRast = curJitter, prevBurned=prevBurned, sizeTolPer = sizeTolPer, outputDir=outputDir, outputSuffix=layerNames[fractionIdx])
    if (writeTempFiles) {
      isBurned <- writeRaster(isBurned, sprintf("%sisBurned_%s", outputDir, layerNames[fractionIdx]),
                              overwrite=T, format='GTiff',
                              options=c("TILED=YES","COPY_SRC_OVERVIEWS=YES","COMPRESS=LZW"))
    }
    
    if (dependentLayers) {
      fireMask <- mask(fireMask, isBurned, maskvalue=0)
    } else {
      fireMask <- mask(fireMask, isBurned, maskvalue=1)
      prevBurned <- prevBurned | isBurned
    }
    
    returnStack <- addLayer(returnStack,isBurned)
  }
  
  if (!is.null(names(fractionArray))) {
    names(returnStack) <- names(fractionArray)
  }
  
  return(returnStack)
}

# create the model environment based on the name of the model. Model environment includes objects
# such as the required fields and saved input stack, the methods to generated the probability maps
# and final severity rastere, etc.
createModelEnv <- function(modelName) {
  modelEnv <- new.env()
  modelEnv$modelName <- modelName
  
  source(paste0(modelName,'.R'), local=modelEnv)
  requiredList <- c('createProbMap',
                    'createProbMapFromDir',
                    'createProbMapFromPrevious')
  for(reqName in requiredList) {
    if (!(exists(reqName, envir = modelEnv))) {
      cat("Model code doesn't include required object:",reqName,"Returning NULL\n")
      return(NULL)
    }
  }

  return(modelEnv)
}
