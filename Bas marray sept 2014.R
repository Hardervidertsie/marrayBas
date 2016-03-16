


rm(list=ls())
options(stringsAsFactors = FALSE )

library(limma) 
library(affy)
library( arrayQualityMetrics )
library("hgu133plus2.db")  # identical:hgu133plus2 array, except for a few extra control probesets, and the fact that they insisted on adding an extra _PM to all the probesets.
library("genefilter")

library("hthgu133pluspmcdf")

cellFilePath <- "H:/R_HOME/R_WORK/Bas/CEL files"
setwd(cellFilePath) # omdat ReadAffy automatisch alle cell files van working directory inlaadt.
dir()      
myData <- ReadAffy(   )
biocLite("hwriter")
setwd("H:/R_HOME/R_WORK/Bas/final analysis sept 2014")

myData@cdfName

myRawPheno <- read.delim("H:/R_HOME/R_WORK/Bas/official affy analysis/corrected sample annotation Bas.csv", sep = ',' )

# correct some typo's in file and split subjectID string
PD <- gsub( "inulin", "insulin", myRawPheno[ ,1 ] )
PD <- cbind( PD, myRawPheno[ ,2] ) 
splitStrings <- strsplit( PD[ ,1 ], "_"  )

myGroup <- gsub( "_[0-9]*$", "", PD[ ,1 ],  )
if (length( levels(factor( myGroup ) ) ) != 30)
{
  stop("grouping went wrong")
}
myPhenoData <- cbind( myGroup , PD[ , 1 ], myRawPheno[ , 2 ],  do.call("rbind", splitStrings) )
colnames( myPhenoData ) <- c( "group", "sampleName", "cellFile", "cell_line", "treatment", "time", "labelID")
myPhenoData <-as.data.frame( myPhenoData)
class( myPhenoData )
myPhenoData$group <- factor(myPhenoData$group)
myPhenoData$sampleName <- factor(myPhenoData$sampleName)
myPhenoData$cellFile <- factor(myPhenoData$cellFile)
myPhenoData$cell_line <- factor(myPhenoData$cell_line)
myPhenoData$treatment <- factor(myPhenoData$treatment)
myPhenoData$time <- factor(myPhenoData$time)
myPhenoData$replicate <- factor(myRawPheno$Bio_Rep)

print(myPhenoData)
sampleNames(phenoData(myData))

#check if level # is correct:
if ( length( levels( myPhenoData$cellFile ) ) != 96 | length( levels( myPhenoData$sampleName ) ) != 96 | 
       length( levels( myPhenoData$cell_line ) ) != 3 | length( levels( myPhenoData$treatment ) ) != 5 | 
       length( levels( myPhenoData$time ) ) != 2 ) 
{
  stop( "Incorrect phenoData loaded or phenodata processing error" )
}

mt <- match( row.names(pData(myData)), myPhenoData$cellFile)

for(i in 1: length(mt)) {
  print(mt[i+1] -mt[i])
}


if ( !all(rownames(pData(myData)) == myPhenoData[mt,'cellFile'])){
  stop(" ")
}

cbind(pData(myData),myPhenoData[mt,]) # check if correct

myPhenoData<-myPhenoData[mt,]
myPhenoData <- lapply(myPhenoData, factor)
myPhenoData$sampleName

cbind(pData(myData),myPhenoData) # check if is correct

pData(myData) <- cbind(pData(myData),myPhenoData)

sampleNames(phenoData(myData)) <- pData(myData)$sampleName
sampleNames(protocolData(myData)) <- pData(myData)$sampleName
sampleNames(assayData(myData)) <- pData(myData)$sampleName

# before pre processing:
intgroup = c(  "treatment", "time","cell_line")
arrayQualityMetrics(expressionset = myData,
                    outdir = "Report_for_rawData",
                    force = TRUE,
                    do.logtransform = TRUE,
                    intgroup=intgroup)


# after pre porcessing:
eSetrma <- rma( myData )

arrayQualityMetrics(expressionset = eSetrma,
                    outdir = "Report_for_rma",
                    force = TRUE,
                    intgroup=intgroup)

intgroup = c(  "cell_line","treatment", "time")
arrayQualityMetrics(expressionset = eSetrma,
                    outdir = "Report_for_rma_cell_lineColor",
                    force = TRUE,
                    intgroup=intgroup)

annotation(eSetrma) <- "hgu133plus2"

featureNames(eSetrma) <- gsub("_PM","", featureNames(eSetrma))

eSetrmaF <- nsFilter( eSetrma,  var.func = sd, filterByQuantile=TRUE, var.cutoff = 0.3, remove.dupEntrez = TRUE, require.entrez=TRUE )

# Annotation: hgu133plus2 
# 
# $filter.log
# $filter.log$numDupsRemoved
# [1] 21348
# 
# $filter.log$numLowVar
# [1] 6033
# 
# $filter.log$numRemoved.ENTREZID
# [1] 13247
# 
# $filter.log$feature.exclude
# [1] 11


# linear modeling
cellLine <- eSetrma$cell_line
treatment <- eSetrma$treatment
time <- eSetrma$time

allCombos <- paste( cellLine, treatment, time, sep = "."  )
allCombos <- factor( allCombos )
designAC <- model.matrix( ~0 + allCombos) # all combo's design matrix
cbind(colnames(designAC ), levels(allCombos))

colnames(designAC ) <- levels(allCombos)

eSetrmaF <- eSetrmaF$eset

# filter genes close to background level
row.means <- rowMeans(exprs(eSetrmaF))
filter.threshold <- quantile(row.means, 0.3)
eSetrmaFF <- eSetrmaF[row.means > filter.threshold, ]
mean(rowMeans(exprs(eSetrmaFF)))
dim(eSetrmaF)


fitAll_AC <- lmFit( eSetrmaFF, designAC )

contr.matrix_AC <- makeContrasts( IGF1R_glargVScontrol_1h = IGF1R.glargine.1h-IGF1R.control.1h,
                                  IGF1R_IGF1VScontrol_1h = IGF1R.IGF1.1h-IGF1R.control.1h,
                                  IGF1R_insulinVScontrol_1h = IGF1R.insulin.1h-IGF1R.control.1h,
                                  IGF1R_X10VScontrol_1h = IGF1R.X10.1h-IGF1R.control.1h,
                                  
                                  IGF1R_glargVScontrol_6h = IGF1R.glargine.6h-IGF1R.control.6h,
                                  IGF1R_IGF1VScontrol_6h = IGF1R.IGF1.6h-IGF1R.control.6h,
                                  IGF1R_insulinVScontrol_6h = IGF1R.insulin.6h-IGF1R.control.6h,  
                                  IGF1R_X10VScontrol_6h = IGF1R.X10.6h-IGF1R.control.6h,
                                  
                                  IRA_glargVScontrol_1h = IRA.glargine.1h-IRA.control.1h,
                                  IRA_IGF1VScontrol_1h = IRA.IGF1.1h-IRA.control.1h,
                                  IRA_insulinVScontrol_1h = IRA.insulin.1h-IRA.control.1h,
                                  IRA_X10VScontrol_1h = IRA.X10.1h-IRA.control.1h,
                                  
                                  IRA_glargVScontrol_6h = IRA.glargine.6h-IRA.control.6h,
                                  IRA_IGF1VScontrol_6h = IRA.IGF1.6h-IRA.control.6h,
                                  IRA_insulinVScontrol_6h = IRA.insulin.6h-IRA.control.6h,
                                  IRA_X10VScontrol_6h = IRA.X10.6h-IRA.control.6h,
                                  
                                  IRB_glargVScontrol_1h = IRB.glargine.1h-IRB.control.1h,
                                  IRB_IGF1VScontrol_1h = IRB.IGF1.1h-IRB.control.1h,
                                  IRB_insulinVScontrol_1h = IRB.insulin.1h-IRB.control.1h,
                                  IRB_X10VScontrol_1h = IRB.X10.1h-IRB.control.1h,
                                  
                                  IRB_glargVScontrol_6h = IRB.glargine.6h-IRB.control.6h,
                                  IRB_IGF1VScontrol_6h = IRB.IGF1.6h-IRB.control.6h,
                                  IRB_insulinVScontrol_6h = IRB.insulin.6h-IRB.control.6h,
                                  IRB_X10VScontrol_6h = IRB.X10.6h-IRB.control.6h,
                                  
                                  
                                  levels = designAC )

fitAll_AC1 <- contrasts.fit(fitAll_AC, contr.matrix_AC)
fitAll_AC_ <- eBayes(fitAll_AC1)

gatherList = list()
for (i in 1 : length(colnames(fitAll_AC_))) {
  
  print(i)
  
  output_AC <- topTable(fitAll_AC_[, colnames(fitAll_AC_)[ i ]], adjust.method = "BH", sort = "p", n = Inf)
  syms <- unlist( mget( rownames(output_AC), env = hgu133plus2SYMBOL, ifnotfound=NA ) )
  
  geneID <- unlist( mget( rownames(output_AC), env = hgu133plus2ENTREZID, ifnotfound=NA ) )
  length((geneID))
  
  if (nrow(output_AC) != length(syms) ){
    stop()
  }
  if (nrow(output_AC) != length(geneID) ){
    stop()
  }
  
  output_AC$syms <- syms
  output_AC$entrez_gi <- geneID
  head(output_AC)
  write.table(output_AC, file = paste(colnames(fitAll_AC_ )[i], "corrected_sampleAnot_sept2014.csv", sep = ""), sep = ",", col.names = NA ) 

  
  sortedOutput_AC <- output_AC[order(as.numeric(output_AC$entrez_gi)), ]

  gatherList[[i]] <- sortedOutput_AC[, "logFC", drop = FALSE]
}

gatherList2 <- do.call("cbind", gatherList)
colnames(gatherList2) <- colnames(fitAll_AC_)
write.table(gatherList2, file = "seperate_samples_log2FC.csv", sep =",", row.names = TRUE, col.names= NA)
write.table((exprs(eSetrmaFF)), file = "seperate_replicate_samples_log2 values.csv", sep =",", row.names = TRUE, col.names= NA)



