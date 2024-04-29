#Generate Weights using Chebyshev
generateWeight <- function(matrix, target_row, target_col, layer=0) {
  numrow <- nrow(matrix)
  numcol <- ncol(matrix)
  
  for (row in 1:numrow) {
    for (col in 1:numcol) {
      # Chebyshev distance
      distance <- max(abs(row - target_row), abs(col - target_col))
      
      if (distance > 0) {
        weight <- 1 / distance
        
        #Changing the threshold based on number of layers
        if (layer==0){ #All layers
          matrix[row, col] <- weight
        } else { #Custom layer
          if (weight < 1/layer){
            matrix[row, col] <- 0
          } else {
            matrix[row, col] <- weight
          }
        }
      } else {
        matrix[row, col] <- 0
      }
    }
  }
  
  return(c(t(matrix)))
}

#Weighted Mean Calculation
imputeWeightedMeanPerColumn <- function(column, weights){
  nonNaIdx <- !is.na(column)
  result <- sum(column[nonNaIdx] * weights[nonNaIdx] / sum(weights[nonNaIdx]))
  return(result)
}

#Weighted Impute Function
RSWeightedImpute <- function(data, rowDim, colDim, targetIndex, numWavelength=814, layer=0){
  #Convert to 2D matrix index
  targetrow <- (targetIndex-1) %/% colDim + 1
  targetcol <- (targetIndex-1) %% colDim + 1
  
  #Generate Weight
  weightPlaceholder <- matrix(0, rowDim, colDim)
  weightMatrix <- generateWeight(weightPlaceholder, targetrow, targetcol, layer)
  
  #Weighted Random Forest
  library(missRanger)
  imputed.wrf <- as.matrix(missRanger(as.data.frame(data), case.weights=weightMatrix, maxiter=5, sample.fraction=0.1, num.trees=80, seed=8566, verbose=0))
  
  #Weighted Mean:
  ImputedValues <- numeric(numWavelength)
  for(i in 1:numWavelength){
    ImputedValues[i] <- imputeWeightedMeanPerColumn(data[,i],weightMatrix)
  }
  imputed.wmean <- data
  imputed.wmean[targetIndex,] <- ImputedValues
  
  result <- list(imputed.wrf, imputed.wmean)
  return(result)
}