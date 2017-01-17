# Original version


pknmodel<-readSIF("./edgeWeightDev/ToyPKNMMB.sif")
midas = readMIDAS("./edgeWeightDev/ToyDataMMB.csv")
cnolist = makeCNOlist(midas,subfield = FALSE)
plotModel(pknmodel, CNOlistToy)
# pknmodel = readSIF("ToyModel.sif")
# cnolist = CNOlist("ToyDataMMB.csv")


model = preprocessing(cnolist, pknmodel)
plotModel(model, CNOlistToy)

# ---------------------- perform the analysis
res = gaBinaryT1(cnolist, model, verbose=FALSE)

# ---------------------- plot the results
cutAndPlot(cnolist, model, list(res$bString))


# ------------   weighted version -------------------
# implementation of edge weighting for the edges

# INPUTing edge weights -------------------------------
# probably the simplest way to insert confidence scores is through the Sif files,
# especially if the prior knowledge network is generated outside of R (pypath)
# -> 4th column in the SIF file: confidence scores
# -> default confidence scores: 0.5  (will give 0 contribution to the cost function)


# Integrate confidence score to: readSif, plotModel, preprocessing -------------
# SIF file with confidence scores, but without AND gates
pknmodel<-readSIF("./edgeWeightDev/ToyPKNMMB_ew.sif")   

# SIF file with confidence scores and AND gates
pknmodel2<-readSIF("./edgeWeightDev/ToyPKNMMB_ew_and.sif")  

# TODO: test readSIF with sif file that contains complicated reactions:
# (more than 2 inputs to AND node and/or more than 1 output)


# Integrate with network plotting without preprocessing ------------------------
# plot model PKN with possible indication on confidence scores on the edges
pknmodel<-readSIF("./edgeWeightDev/ToyPKNMMB_ew.sif")   
midas = readMIDAS("./edgeWeightDev/ToyDataMMB_noP90RSK.csv")
cnolist = makeCNOlist(midas,subfield = FALSE)
plotModel(pknmodel, cnolist,edgeLabel = "confScore")


pknmodel2<-readSIF("./edgeWeightDev/ToyPKNMMB_ew_and.sif")  
midas = readMIDAS("./edgeWeightDev/ToyDataMMB_noP90RSK.csv")
cnolist = makeCNOlist(midas,subfield = FALSE)
plotModel(pknmodel2, cnolist,edgeLabel = "confScore")


# Integrate with preprocessing -------------------------------------------------

# preprocessing: take care of scores  of removed and added edges
# simpler toy model:
pknmodel<-readSIF("./edgeWeightDev/ToyPKNMMB_ew.sif")   
midas = readMIDAS("./edgeWeightDev/ToyDataMMB_noP90RSK.csv")
cnolist = makeCNOlist(midas,subfield = FALSE)
# cutNONC: confScores of non-controlable edges are simply removed. 
# compressModel: resulted edges are the  mean of orginal edges 
#                we could keep the lower weight or higher too
# expandGates: add weights for and gates.
plotModel(pknmodel, cnolist,edgeLabel = "confScore")
model = preprocessing(cnolist, pknmodel,verbose = TRUE)
plotModel(model, cnolist,edgeLabel = "confScore")


# tested until here.

# complex liver model:
pknmodel<-readSIF("PKN-LiverDREAM_confScore.sif")   
plotModel(pknmodel, cnolist,edgeLabel = "confScore")
midas = readMIDAS("MD-LiverDREAM.csv")
cnolist = makeCNOlist(midas,subfield = FALSE)

model = preprocessing(cnolist, pknmodel)
plotModel(model, cnolist,edgeLabel = "none")

plotModel(model, cnolist,edgeLabel = "confScore")

writeSIF(model,"ToyPKNMMB_ew_and.sif")

# ---------------------- perform the analysis
res = gaBinaryT1(cnolist, model, verbose=FALSE)

# ---------------------- plot the results
cutAndPlot(cnolist, model, list(res$bString))
