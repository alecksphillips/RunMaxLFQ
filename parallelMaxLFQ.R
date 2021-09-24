#!Rscript --vanilla
library(data.table)
library(fst)
library(progress)
library(pbapply)
library(parallel)

source("processBatch.R")

if (!require("diann")){
  if(!require("devtools")){
    install.packages("devtools")
  }
  library(devtools)
  install_github("https://github.com/vdemichev/diann-rpackage")
}
library(diann)

args <- commandArgs(trailingOnly = T)

if (length(args) == 0){
  stop("Require path of input evidence.txt/peptides.txt at least")
}

if (length(args) < 2){
  args[2] = "proteinGroups.txt"
}

if (length(args) < 3){
  args[3] = "output.txt"
}

doSerial <- T
if (length(args) < 4){
  args[4] = "serial"
} else {
  if(!(args[4] %in% c("serial", "parallel"))){
    stop("Specify either `serial` or `parallel <numBatches> <numWorkers>` ")
  }
}

nBatches <- 1
nWorkers <- 1

if (args[4] == "parallel"){
  doSerial <- F
  if (length(args) < 6){
    stop("Specify <numBatches> <numWorkers>")
  } else {
    nBatches <- strtoi(args[5])
    nWorkers <- strtoi(args[6])
  }
}

#Read evidence.txt input
data <- fread(args[1], sep = "\t")

proteinGroups <- fread(args[2], sep = "\t")

#Throw away shared peptides
data <- data[!grepl(";", `Protein group IDs`),]


discardUnmodifiedCounterparts <- T
#Discard unmodified counterparts?
if (discardUnmodifiedCounterparts){
  tmp <- unique(data[, c("Sequence", "Modifications", "Modified sequence")])
  tmp <- tmp[, {
    keep <- c(T)
    if (.N > 1){

      if ("Unmodified" %in% .SD$Modifications){
        keep <- .SD$Modifications == "Unmodified"
      } else {
        keep <- rep(c(T), .N)
      }
    }
    list(
      "Modifications" = .SD[keep,]$`Modifications`,
      "Modified sequence" = .SD[keep,]$`Modified sequence`
    )
  }, by = c("Sequence")]
  data <- data[`Modified sequence` %in% tmp$`Modified sequence`,]
}

sampleTable <- unique(data[, c("Experiment")])
sampleTable[, sampleID :=seq(1,.N)]
data <- merge(data, sampleTable, by = "Experiment")

groupTable <- unique(data[ , c("Protein group IDs")])
groupTable[, groupID := seq(1,.N)]
data <- merge(data, groupTable, by = "Protein group IDs")


data[, feature := paste0(`Modified sequence`, `Charge`)]
featureTable <- unique(data[, c("Modified sequence", "Charge", "feature")])
featureTable[, featureID := seq(1,.N)]
data <- merge(data, featureTable, by = c("Modified sequence", "Charge", "feature"))

dataForProcessing <- data[, c("sampleID", "groupID", "featureID", "Intensity")]

output <- data.table()

if (doSerial){
  output <- rbindlist(lapply(seq(1,max(groupTable$groupID)), function(i){
    x <- dataForProcessing[groupID == i,]

    outputTable <- tryCatch({
      processedData <- diann::diann_maxlfq(
        x,
        sample.header = "sampleID",
        group.header = "groupID",
        id.header = "featureID",
        quantity.header = "Intensity"
      )
      outputTable <- as.data.table(
        reshape2::melt(
          processedData,
          value.name = "LFQ Intensity"
        )
      )
      setnames(outputTable, c("Var1","Var2"), c("groupID","sampleID"))
      return(outputTable)

    }, error = function(cond){
      outputTable <- x[, c("groupID","sampleID", "Intensity")]
      setnames(outputTable, c("Intensity"), c("LFQ Intensity"))
      return(outputTable)
    })

    outputTable[, c("groupID", "sampleID", "LFQ Intensity")]
  }))


} else {
  #Output from batches ends up being 1000s of files, put in a separate directory
  workingDirectory <- tempdir()

  #FST files allow for random access so we can process the data in batches in parallel
  datafilename <- "data.fst"

  #Create id of protein groups
  setorder(dataForProcessing, groupID)


  datafilepath <- file.path(workingDirectory,datafilename)
  if(file.exists(datafilepath)){
    file.remove(datafilepath)
  }
  write_fst(dataForProcessing, datafilepath)

  #Calculate indexes of where start and end of data is for each protein group
  #in the table for random access with fst
  dataindex <- dataForProcessing[, .(
    groupID = unique(groupID), from = .I[!duplicated(groupID)], to = .I[rev(!duplicated(rev(groupID)))]
  )]
  dataindex[,workingDirectory := workingDirectory]
  dataindex[,file := file.path(datafilename),]
  dataindex[,n := to - from + 1,]
  dataindex[,index := seq(1,.N),] #Use row of index file to avoid having to write ProteinID as a factor in each batch (factors mean that every file stores a list of all levels which wastes space)
  if (file.exists("dataindex.fst")){
    file.remove("dataindex.fst")
  }
  write_fst(dataindex, file.path("dataindex.fst"))

  nGroups <- nrow(dataindex)

  #Clear these big tables from memory for sanity
  data <- NULL
  dataindex <- NULL

  message(paste0("Starting parallel processing of ", nGroups, " proteins in ", nBatches, " batches over ", nWorkers, " workers.\n"))

  t <- system.time({
    #Create workers for parallel processing
    cl <- makeCluster(nWorkers)

    #Load packages on workers
    clusterEvalQ(cl, {
      library(data.table)
      library(fst)
      library(diann)
    })

    #Export nPatients and nBatches to workers
    clusterExport(cl, "nGroups")
    clusterExport(cl, "nBatches")

    #Export function to workers
    clusterExport(cl, "processBatch")

    #List of batches
    x.list <- sapply(seq(1,nBatches), list)

    #Process list of batches on workers
    #pblapply(x.list, processBatch, cl = cl)
    lapply(x.list, processBatch)

    #parLapply(cl = cl, x.list, processBatch)
    parallel::stopCluster(cl)
  })

  message(paste0("Elapsed: ", as.numeric(t[3]), "s"))

  file.remove("dataindex.fst")

  #Combine outputs from batches
  output <- rbindlist(lapply(seq(1,nBatches), function(i){
    tmp <- read_fst(file.path(workingDirectory,paste0(i,".fst")), as.data.table = T)
    tmp
  }))

}

output <- merge(output, groupTable, by = "groupID")
output[,sampleID := strtoi(sampleID)]
output <- merge(output, sampleTable, by = "sampleID")

output[, sampleID := NULL]
output[, groupID := NULL]

setnames(output, c("Protein group IDs"), c("id"))
output[, id := strtoi(id)]
output[, `LFQ Intensity` := as.double(`LFQ Intensity`)]

output <- dcast(output, value.var = "LFQ Intensity", id ~ Experiment)
output <- merge(proteinGroups[, c("id", "Protein IDs", "Majority protein IDs", "Protein names", "Gene names", "Fasta headers")], output, by= "id")
setorder(output, id)

#Write to tsv
fwrite(output, file = args[3], sep = "\t")
