make_RNAseq_longtable <- function(exptsheet, thisexpt){
  exptsubsheet <- exptsheet[exptsheet$Experiment==thisexpt,]
  
  RNAdf <- data.frame()
  for (i in c(1:nrow(exptsubsheet))){
    thisfile <- exptsubsheet$RNAseqFile[i]
    thistime <- exptsubsheet$Time[i]
    f<-read.csv(thisfile, header=T, stringsAsFactors = F)
    f$log2FoldChange <- as.numeric(as.character(f$log2FoldChange))
    f$Sig <- !is.na(f$log2FoldChange) & abs(f$log2FoldChange)>1 & 
      !is.na(f$padj) & f$padj < 0.05
    f$Time <- thistime
    RNAdf <- rbind(RNAdf, f[,c('Gene', 'log2FoldChange', 'padj', 'Sig', 'Time')])
    
  }
  return(RNAdf)
  
  
}

