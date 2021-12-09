### 2021.10.19
nsc.df.bySample <- rbind(data.frame('NSC'=NA, 'sample'=NA, 'factor'=NA, 'type'=NA, 'type_index'=NA)) 
for (TF in c('ATF1','ELF1')){
  for (cell.num in c('FitCR', 'FlagCR', 'wtCR', 'ChIP')){
    cell.num.2 <- paste0(cell.num, '_1e5')
    if (cell.num == 'FitCR') {
      type <- 'FitCUT&RUN'
      type_index <- 4
    }else if (cell.num == 'FlagCR'){
      type <- 'Flag-CUT&RUN'
      type_index <- 2
    }else if (cell.num == 'wtCR'){
      type <- 'CUT&RUN'
      type_index <- 3
    }else {
      type <- 'ChIP-seq'
      type_index <- 1
      cell.num.2 <- 'ChIP'
    }
    
    load(paste0('/mnt/Storage/home/wangyiman/fitCUTRUN/revision_1/analysis/TF/comparison/NSC_bySample/',TF,'_',cell.num.2,'/',TF,'_', cell.num,'.filt.sample.15.SE.tagAlign.gz.cc.RData'),
         ex <- new.env())
    nsc <- ex$crosscorr$peak$y / ex$crosscorr$min.cc$y
    print(nsc)
    nsc.df.newline <- data.frame(nsc, paste0(TF, '_', cell.num), TF, type, type_index)
    colnames(nsc.df.newline) <- c('NSC', 'sample', 'factor', 'type', 'type_index')
    nsc.df.bySample <- rbind(nsc.df.bySample, nsc.df.newline)
  }
}
write.csv(nsc.df.bySample, row.names = FALSE, quote = FALSE, file = '~/fitCUTRUN/revision_1/analysis/TF/comparison/NSC_bySample/NSC_df.csv')



