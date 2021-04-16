saveDir <- '../out/Rdata'
filename <- c('study-1', 'study-2', 'study-3')

lapply(filename, function(x){
  fileDir <- paste0(saveDir, '/', x, '.Rdata')
  scriptDir <- paste0(x, '.R')
  if(!file.exists(fileDir)){
    rm(list=setdiff(ls(), c("saveDir", "filename", "scriptDir", "fileDir")))
    source(file = scriptDir)
    save.image(file = fileDir)
  }
})
