make_RNAseq_longtable <- function(exptsheet, thisexpt){
  exptsubsheet <- exptsheet[exptsheet$Experiment==thisexpt,]
  
  RNAdf <- data.frame()
  for (i in c(1:nrow(exptsubsheet))){
    thisfile <- exptsubsheet$DataFile[i]
    thistime <- exptsubsheet$Time[i]
    f<-read.csv(thisfile, header=T, stringsAsFactors = F)
    f$Value <- as.numeric(as.character(f$Value))
    f$Time <- thistime
    
    if ('Sig' %in% names(f)){
      f$Sig <- as.logical(f$Sig)
      RNAdf <- rbind(RNAdf, f[,c('Gene', 'Value', 'Sig', 'Time')])
      
    }else{
      f$Sig <- !is.na(f$Value) & abs(f$Value)>1 & 
        !is.na(f$padj) & f$padj < 0.05
      RNAdf <- rbind(RNAdf, f[,c('Gene', 'Value', 'padj', 'Sig', 'Time')])
    }

    
  }
  return(RNAdf)
  
  
}

