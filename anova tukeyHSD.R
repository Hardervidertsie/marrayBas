setwd("H:/R_HOME/R_WORK/Bas/official affy analysis")
rm(list=ls())


rawData <- read.table( file = "affy_defaul_rma_normalized_expression_data_Bas_ter_Braak.txt", sep = "\t", header = T)

rawData <- rawData[!is.na(rawData$entrezGeneID), ]
head(rawData)
rawData<- rawData[, -c(2,3)]

rawDataLong <- melt(rawData, id.vars = c("probeID") ) 
rawDataLong <- melt(rawData, id.vars = c("probeID", "entrezGeneID", "syms") ) 
head(rawDataLong)

rawDataLong$receptor <- str_extract( rawDataLong$variable, "^(I[RABGF1R]{2,4})") 
rawDataLong$time <- gsub("_", "", str_extract( rawDataLong$variable,  "_[16]{1}h_") )
rawDataLong$geneTime <- paste(rawDataLong$probeID,rawDataLong$time, sep = "_")

# remove control

indC <- grepl("control", rawDataLong$variable)
rawDataLongNoControl <- rawDataLong[!indC, ]
rawDataLongNoControl$variable <-  factor(rawDataLongNoControl$variable)
head(rawDataLongNoControl)
unique(rawDataLongNoControl$variable)

rawDataLong$replicate <-  str_match(rawDataLong$variable, "(_[1234]{1})$" )[,1]


rawDataLongNoControl$replicate <- str_match(rawDataLongNoControl$variable, "(_[1234]{1})$" )[,1]
head(rawDataLongNoControl)
lapply(rawDataLongNoControl, class)
unique(rawDataLongNoControl$replicate)

head(rawDataLongNoControl)

rawDataLongNoControlMeanTreats <- aggregate( value ~ geneTime + receptor + replicate, data = rawDataLongNoControl, FUN = mean)

head(rawDataLongNoControlMeanTreats)





myDataList <- split(rawDataLongNoControlMeanTreats, rawDataLongNoControlMeanTreats$geneTime   )
myDataList[[1]]
bas.out.anov <- lapply(myDataList, function(x) aov(value ~ receptor, data = x ))
bas.out.anov.s <- lapply(bas.out.anov, summary)

unlistPVals <- unlist(
  lapply(
  lapply(
    lapply(
      bas.out.anov.s, function(x) 
        {
        lapply(x, "[[", "Pr(>F)")
        }
      ), 
    function(xx) xx[[1]] 
    ), 
  function(xxx) xxx[[1]]
  )
  )

head(unlistPVals)

indSign <- unlistPVals <= 0.05
myDataSign <- bas.out.anov.s[indSign]
myDataSign <- lapply( myDataSign, unlist)
myDataSignDF <- do.call('rbind', myDataSign)
myDataSignDF <- as.data.frame(myDataSignDF)

length(bas.out.anov.s)
length(bas.out.anov)
signModelData <- bas.out.anov[indSign]

tukeyOutput <- lapply(signModelData, TukeyHSD)
tukeyOut.DF <- lapply(tukeyOutput, "[[", "receptor")
tukeyOut.DF <- do.call("rbind", tukeyOut.DF)
tukeyOut.DF<- as.data.frame(tukeyOut.DF)
ids <- unlist(names(tukeyOutput))
tukeyOut.DF$probeID <- rep(ids, each = 3)

tukeyOut.DF$tukeyTest <- rownames(tukeyOut.DF)
myDataSignDF$probeID <- rownames(myDataSignDF)

finalData <- merge(myDataSignDF, tukeyOut.DF, by = "probeID", sort = FALSE)
head(finalData)

head(rawDataLong)
finalData$probeIDTime <- finalData$probeID
finalData$probeID <-gsub("_[16]h", "", finalData$probeIDTime)
head(finalData)

rawDataLong$control <- as.numeric(grepl("control", rawDataLong$variable))
head(rawDataLong)
metaD <- rawDataLong[, 1:3]
summaryByMeanComp_Control <- aggregate( value ~ probeID + control + receptor + time, data = rawDataLong, FUN = mean)
head(metaD)
time1 <- summaryByMeanComp_Control[summaryByMeanComp_Control$time == "1h",]
time2 <- summaryByMeanComp_Control[ summaryByMeanComp_Control$time == "6h", ]

head(finalData)
head(metaD)
metaD <- unique(metaD)
almostFinalData<- merge(finalData, metaD, by = 'probeID', sort = FALSE)

almostFinalDataT1 <- almostFinalData[ grepl("_1h",almostFinalData$probeIDTime), ]

almostFinalDataT2 <- almostFinalData[ grepl("_6h",almostFinalData$probeIDTime), ]

head(almostFinalDataT1)
time1_treat <- time1[ time1$control == 0,]
time1_c <- time1[ time1$control == 1,]
time1<- merge(time1_treat, time1_c, by = c("probeID", "receptor"))

time2_treat <- time2[ time2$control == 0,]
time2_c <- time2[ time2$control == 1,]
time2<- merge(time2_treat, time2_c, by = c("probeID", "receptor"))
head(time2)

dim(almostFinalDataT1)
finalDataWithOrigData_time1[1:10,]

dim(time1)
time1.compounds <- dcast(time1, probeID~receptor, value.var = c("value.x"))
time1.control <- dcast(time1, probeID~receptor, value.var = c("value.y"))
colnames(time1.compounds)[2:4] <- paste0(colnames(time1.compounds)[2:4], "compounds")
colnames(time1.control)[2:4] <- paste0(colnames(time1.control)[2:4], "control")

time1Final <- merge(time1.compounds, time1.control, by = "probeID")

time2.compounds <- dcast(time2, probeID~receptor, value.var = c("value.x"))
time2.control <- dcast(time2, probeID~receptor, value.var = c("value.y"))
colnames(time2.compounds)[2:4] <- paste0(colnames(time2.compounds)[2:4], "compounds")
colnames(time2.control)[2:4] <- paste0(colnames(time2.control)[2:4], "control")

time2Final <- merge(time2.compounds, time2.control, by = "probeID")



dim(time1)
finalDataWithOrigData_time1 <- merge(almostFinalDataT1, time1Final, by = "probeID", all.x = TRUE, all.y = FALSE)
write.table(finalDataWithOrigData_time1, file = "anova_Tukey_time_1hrs.txt", sep = "\t", row.names=FALSE)
finalDataWithOrigData_time2 <- merge(almostFinalDataT2, time2Final, by = "probeID", all.x = TRUE, all.y = FALSE)
write.table(finalDataWithOrigData_time2, file = "anova_Tukey_time_6hrs.txt", sep = "\t", row.names=FALSE)

head(finalDataWithOrigData_time1)







aov.out.l <- lapply(INSList, function(x) aov( expression ~ receptor, data = x))
aov.out.l.s <- lapply( aov.out.l, summary )
unlist(lapply(lapply(lapply(aov.out.l.s, function(x) { lapply(x, "[[", "Pr(>F)")}), function(xx) xx[[1]] ), function(xxx) xxx[[1]]))



# p waarde van anova en d 3 tukey test p waarden
# plak ernaast gemmidelde expressie van treatments en vehicle

probe_ID  symbol  entrezID anovap tukey1 verschil1_t1, verschil2_t1 tukey2, verschil1_t2, verschil2_t2  tukey3, verschil1_t3, verschil2_t3 averExprTr averExprDMSO


?melt
melt

dataPlot24.w <- dcast(dataPlot24,  treatment~ symbols_t, value.var= "logFC")



aov.out.l.s[[1]][[1]][["Pr(>F)"]]


test <- aov.out.l.s[[1]]

test[[1]]$'Pr(>F)'[[1]]



lapply(aov.out.l.s , function(x) x[[1]]$'' )   


test <- str(aov.out.l.s)

lapply( aov.out.l.s, "[[", "Pr(>F)")
  
tests[[1]]$'Pr(>F)'

lapply( aov.out.l.s, "[[", "\'Pr(>F)\'")



test <- unlist(aov.out.l.s)
names(test)

summary(aov.out.l[[1]])[["Pr(>F)"]]
