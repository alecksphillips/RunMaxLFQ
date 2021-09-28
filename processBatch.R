processBatch <- function(batchNumber) {

  #Batch 1 includes index 1, 1 + nBatches, 1 + 2*nBatches ...
  indexes <- seq(batchNumber, nGroups, by=nBatches)

  dataindex <- fst::read_fst(
    file.path("dataindex.fst"),
    as.data.table = T
  )

  #For each protein group
  tmp <- lapply(indexes, function(i) {

    #Read indexes of protein data from the index file
    item <- dataindex[index == i,]

    #Read corresponding protein data
    dt <- fst::read_fst(
      file.path(item[,workingDirectory], item[,file]),
      columns = NULL,
      from = item[1,from],
      to = item[1,to],
      as.data.table = T
    )

    doMaxLFQ(dt, useDIANN)
  })


  #bind tables of proteins in this batch to one table
  tmp <- data.table::rbindlist(tmp, use.names = T)

  #save this batch to a file
  fstFile <- file.path(dataindex[1,workingDirectory], paste0(batchNumber,".fst"))
  if(file.exists(fstFile)){
    file.remove(fstFile)
  }
  fst::write_fst(tmp, fstFile)
  tmp <- NULL

  return(invisible(NULL))

}
