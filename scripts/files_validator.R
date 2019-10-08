# returns error messages
files_validator <- function() {
  outcome <- c()
  exptsheet<- tryCatch({
    read.csv('./data/exptsheet.csv', header=T, stringsAsFactors = F)
  }, error = function(err){
    exptnotfound_Error <- 'Experiment sheet not found, please make sure to have one .csv file in the data directory'
    #print(exptnotfound_Error)
    #print(err)
    outcome <- c(outcome, exptnotfound_Error)
  })
  
  ## Check experiment sheet 
  exptsheetnames <- names(exptsheet)
  requiredcolumns <- c('Name','Experiment', 'DataFile', 'MetadataFile', 'Strain','Time')
  missingcols <- c()
  for (col in requiredcolumns){
    if ((col %in% exptsheetnames)==FALSE){
      missingcols <- c(missingcols, col)
    }
  }
  if (length(missingcols)>0){
    exptsheet_columns_Error <- paste(c('ERROR: Exptsheet has the following columns missing: ', missingcols), sep=' ')
    #print(exptsheet_columns_Error)
    outcome <- c(outcome, exptsheet_columns_Error)
  }else{
    #print('exptsheet.csv has all necessary columns!')
    }
  
  # is each experiment uniquely named? 
  dupnames <- duplicated(exptsheet$Name)
  if (sum(dupnames)>0){
    exptsheet_uniquenames_Error <- paste(c('ERROR: Experiment names should be unique - looks like the following names are repeated', c(exptsheet$Name[dupnames])), sep=" ")
    #print(exptsheet_uniquenames_Error)
    outcome <- c(outcome, exptsheet_uniquenames_Error)
  }else{
    #print('Experiment names are unique!')
    }
  
  ## Check data files
  datafiles <- unique(exptsheet$DataFile)
  for (datafile in datafiles){
    myfile <- read.csv(datafile, header=T, stringsAsFactors = F)
    if (('Gene' %in% names(myfile)) & ('Value' %in% names(myfile)) &
        (('padj' %in% names(myfile))|('Sig' %in% names(myfile)))){
      #print(paste0('Data File ',datafile, ' looks good!'))
    }else{
      datafile_Error <- paste0('ERROR: Please ensure the following file has at least a "Gene" and a "Value" column, and one of "padj" or "Sig" column: ', datafile)
      #print(datafile_Error)
      outcome <- c(outcome, datafile_Error)
    }
  }
  
  
  ## Check metadata files
  metadatafiles <- unique(exptsheet$MetadataFile)
  for (metadatafile in metadatafiles){
    myfile <- read.csv(metadatafile, header=T, stringsAsFactors = F)
    if (('Gene' %in% names(myfile)) ){
      #print(paste0('Metadata file ',metadatafile, ' looks good!'))
    }else{
      metadatafile_Error <- paste0('ERROR: Please ensure the following file has at least a "Gene" column: ', metadatafile)
      #print(metadatafile_Error)
      outcome <- c(outcome, metadatafile_Error)
    }
  }
  
  
  ## Check networks
  netfiles <- Sys.glob('./data/network/*_Edges.csv')
  for (netfile in netfiles){
    edges <- read.csv(netfile, header=T, stringsAsFactors = F)
    if (('source' %in% names(edges)) && ('target' %in% names(edges))){
      #print(paste0('Network file ',netfile, ' looks good!'))
    }else{
      networkfile_Error <- paste0('ERROR: Please ensure the following file has "source" and "target" columns: ', netfile)
      #print(networkfile_Error)
      outcome <- c(outcome, networkfile_Error)
    }
  }
 
  return(outcome) 
}