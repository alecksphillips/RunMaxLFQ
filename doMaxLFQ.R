doMaxLFQ <- function(x, useDIANN = F){

  #Default return
  outputTable <- x[, c("groupID","sampleID", "Intensity")]
  setnames(outputTable, c("Intensity"), c("LFQ Intensity"))

  x[,Intensity := bit64::as.integer64(Intensity)]
  x <- x[!is.na(Intensity),]

  if (useDIANN){
    outputTable <- tryCatch({
      processedData <- diann::diann_maxlfq(
        x,
        sample.header = "sampleID",
        group.header = "groupID",
        id.header = "featureID",
        quantity.header = "Intensity"
      )
      outputTable <- data.table(
        groupID = x$groupID[1],
        sampleID = colnames(processedData),
        `LFQ Intensity` = bit64::as.integer64(processedData)
      )

      return(outputTable)

    }, error = function(cond){
      outputTable <- x[, c("groupID","sampleID", "Intensity")]
      setnames(outputTable, c("Intensity"), c("LFQ Intensity"))

      return(outputTable)
    })

  } else {
    #Use MaxLFQ from iq
    if (nrow(x) > 0){
      tmp <- unique(x)
      tmp[, Intensity := as.double(Intensity)]
      tmp <- dcast(tmp, featureID ~ sampleID, value.var = "Intensity", fun.aggregate = sum)
      sampleIDs <- colnames(tmp)[2:ncol(tmp)]
      tmp <- as.matrix(tmp[,2:ncol(tmp)])

      processedData <- iq::maxLFQ(
        tmp
      )
      outputTable <- data.table(
        groupID = x$groupID[1],
        sampleID = sampleIDs,
        `LFQ Intensity` = bit64::as.integer64(processedData$estimate)
      )
    }
  }

  outputTable[, `LFQ Intensity` := bit64::as.integer64(`LFQ Intensity`)]
  outputTable[, c("groupID", "sampleID", "LFQ Intensity")]
}
