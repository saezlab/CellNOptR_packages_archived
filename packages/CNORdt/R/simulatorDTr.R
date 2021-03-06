#
#  This file is part of the CNO software
#
#  Copyright (c) 2011-2012 - EBI
#
#  File author(s): CNO developers (cno-dev@ebi.ac.uk)
#
#  Distributed under the GPLv2 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-2.0.html
#
#  CNO website: http://www.ebi.ac.uk/saezrodriguez/software.html
#
##############################################################################
# $Id$

simulatorDTr <- function(CNOlist, model, simList, indexList, boolUpdates) {

	nSp <- dim(model$interMat)[1]
	nReacs <- dim(model$interMat)[2]	
	nCond <- dim(CNOlist$valueStimuli)[1]
	yBool = array(dim=c(nCond, nSp, boolUpdates))
	yBool[,,1] = 0

	if(is.null(dim(model$interMat))) { 
		nSp <- length(model$interMat)
		nReacs <- 1
	}

	# this holds, for each sp, how many reactions have that sp as output
	# + other function definitions here
	endIx <- rep(NA,nSp)
		for(i in 1:nSp){
			endIx[i] <- length(which(simList$maxIx == i))
		}
	
	fillTempCube <- function(x) {
		cMatrix <- matrix(data=x, nrow=nReacs, ncol=nCond)
		cVector <- apply(cMatrix, 1, function(x){return(x)})
		return(cVector)
	}
		
	minNA <- function(x) {
		if(all(is.na(x))) {
			return(NA)
		} else {
			return(min(x, na.rm=TRUE))
		}
	}
				
	compOr <- function(x) {
		if(all(is.na(x[which(simList$maxIx == s)]))) {
			res <- NA
		} else {
			res <- max(x[which(simList$maxIx == s)], na.rm=TRUE)
		}
		return(res)
	}
	
	# this value is used to test the stop condition for difference between 2 iterations
	testVal <- 1E-3

	# create an initial values matrix	
	initValues <- matrix(data=NA, nrow=nCond, ncol=nSp)
	colnames(initValues) <- model$namesSpecies

	# set the initial values of the stimuli
	initValues[,indexList$stimulated] <- CNOlist$valueStimuli

	# flip the inhibitors so that 0 = inhibited / 1 = noninhibited
	valueInhibitors <- 1-CNOlist$valueInhibitors
	valueInhibitors[which(valueInhibitors == 1)] <- NA

	# set the initial values of the inhibited species: 0 if inhibited, untouched if not inhibited
	initValues[,indexList$inhibited] <- valueInhibitors

	# initialise main loop
	newInput <- initValues

	# main loop

	for(count in 2:boolUpdates) {
		
		outputPrev <- newInput
		# this is now a 2 column matrix that has a column for each input (column in finalCube)
		# and a set of rows for each reac (where a set contains as many rows as conditions)
		# all concatenated into one long column
				
		if(nReacs > 1) {
			tempStore <- apply(simList$finalCube, 2, function(x){return(outputPrev[,x])})
			tempIxNeg <- apply(simList$ixNeg, 2, fillTempCube)
			tempIgnore <- apply(simList$ignoreCube, 2, fillTempCube)
		} else {
			tempStore <- outputPrev[,simList$finalCube]
			tempIxNeg <- matrix(simList$ixNeg, nrow=nCond, ncol=length(simList$ixNeg), byrow=TRUE)
			tempIgnore <- matrix(simList$ignoreCube, nrow=nCond, ncol=length(simList$ignoreCube), byrow=TRUE)
		}

		# set to NA the values that are "dummies", so they won't influence the min
		tempStore[tempIgnore] <- NA

		# Flip the values that enter with a negative sign
		tempStore[tempIxNeg] <- 1-tempStore[tempIxNeg]
	
		# compute all the ands by taking, for each gate, the min value across the inputs of that gate
		if(nReacs > 1) {
			
			outputCube <- apply(tempStore, 1, minNA)

			# outputCube is now a vector of length (nCond*nReacs) that contains the input of each reaction in
			# each condition, concatenated as such allcond4reac1, allcond4reac2, etc...			# this is transformed into a matrix with a column for each reac and a row for each cond		
			outputCube <- matrix(outputCube, nrow=nCond, ncol=nReacs)
			# go through each species, and if it has inputs, then take the max across those input reactions
			# i.e. compute the ORs
		
			for(s in 1:nSp){
				if(endIx[s] != 0){
					newInput[,s] <- apply(outputCube, 1, compOr)
				}
			}
		
		} else {
			outputCube <- ifelse(all(is.na(tempStore)), NA, min(tempStore,na.rm=TRUE))
			newInput[,simList$maxIx] <- outputCube
		}
	
		# reset the inhibitors and stimuli
		for(stim in 1:length(indexList$stimulated)) {
			stimM <- cbind(CNOlist$valueStimuli[,stim], newInput[,indexList$stimulated[stim]])
			maxNA <- function(x) {
				return(max(x, na.rm=TRUE))
			}
			stimV <- apply(stimM, 1, maxNA)
			newInput[,indexList$stimulated[stim]] <- stimV
		}
		
		valueInhibitors <- 1-CNOlist$valueInhibitors
		newInput[,indexList$inhibited] <- valueInhibitors * newInput[,indexList$inhibited]
	
		# replace NAs with zeros to avoid having the NA penalty applying to unconnected species
		readout <- newInput
		readout[is.na(readout)] <- 0
		yBool[,,count] = readout		
	}

	return(yBool)
}

