performSingleSequencePredictions=function(file,allele,peptidelength,affcutoff,proccutoff,exprcutoff){
  
}

performPairedSequencePredictions=function(file,allele,peptidelength,affcutoff,proccutoff,exprcutoff){
  # prepare empty table
  epitopePredictions=data.table()
  epitopePredictionsWithFiltersApplied=data.table()
  
  # get some info on dataset
  filePath=file
  fileName = gsub(pattern = "\\..+$",
                  replacement = "",
                  x = basename(filePath))
  dirPath = dirname(filePath)
  
  # import kitchensink data & clean up
  kitchensink=unique(fread(input = filePath))
  
  sampleId=toupper(gsub(pattern = "_.+$|-.+$",
                        replacement = "",
                        x = fileName))
  variantInfo=returnProcessedVariants(id = sampleId,
                                      variants = kitchensink)
  
  progressBar = txtProgressBar(min = 0,
                               max = nrow(variantInfo),
                               width = 100,
                               style = 3)
  
  for(i in 1:nrow(variantInfo)){
    # for each variant line, make list tumor peptides which are different from normal (and corresponding normal peptides)
    # and make vector containing normal and tumor peptide stretches
    peptideList=buildPeptideList(variant = variantInfo[i,],
                                 peptidelength = peptideLength)
    peptideStretchVector=c(variantInfo[i,]$peptidecontextnormal,variantInfo[i,]$peptidecontexttumor)
    
    # if no tumor peptides found, move to next line
    if(nrow(peptideList[[2]])<1){
      setTxtProgressBar(progressBar, i)
      next
    }
    
    # perform normal and tumor affinity & processing predictions
    normalAndTumorPredictions=foreach(k=1:2) %dopar% {
      performParallelPredictions(peptides = peptideList[[k]],
                                 peptidestretch = peptideStretchVector[k],
                                 allele = hlaType,
                                 peptidelength = peptideLength)
    }
    
    # rename columns
    setnames(x = normalAndTumorPredictions[[1]],
             old = c("peptide",paste0(hlaType,"affinity"),"c_term_aa","processing_score"),
             new = c("normal_peptide",paste0("normal_",hlaType,"affinity"),"normal_c_term_aa","normal_processing_score"))
    
    setnames(x = normalAndTumorPredictions[[2]],
             old = c("peptide",paste0(hlaType,"affinity"),"c_term_aa","processing_score"),
             new = c("tumor_peptide",paste0("tumor_",hlaType,"affinity"),"tumor_c_term_aa","tumor_processing_score"))
    
    # apply various cutoffs
    normalPredictionsWithFiltersApplied=subset(x = normalAndTumorPredictions[[1]],
                                               subset = normalAndTumorPredictions[[1]][[paste0("normal_",hlaType,"affinity")]] <= affinityCutoff &
                                                 normal_processing_score >= processingCutoff &
                                                 rna_expression_fpkm > expressionCutoff)
    
    tumorPredictionsWithFiltersApplied=subset(x = normalAndTumorPredictions[[2]],
                                              subset = normalAndTumorPredictions[[2]][[paste0("tumor_",hlaType,"affinity")]] <= affinityCutoff &
                                                tumor_processing_score >= processingCutoff &
                                                rna_expression_fpkm > expressionCutoff)
    
    # determine self-sim
    
    ### use self-sim which supports 9, 10 and 11mers; get from Marit
    
    # merge all info
    if (nrow(tumorPredictionsWithFiltersApplied)>0){
      mergedTumorPredictionsWithFiltersApplied=merge(x = tumorPredictionsWithFiltersApplied,
                                                     y = normalAndTumorPredictions[[1]],
                                                     by = c("variant_id","gene_symbol","rna_expression_fpkm","c_term_pos"),
                                                     all.x = TRUE)
      
      epitopePredictionsWithFiltersApplied=rbindlist(list(epitopePredictionsWithFiltersApplied,mergedTumorPredictionsWithFiltersApplied))
    }
    
    if (nrow(normalAndTumorPredictions[[2]])>0){
      mergedPredictions=merge(x = normalAndTumorPredictions[[2]],
                              y = normalAndTumorPredictions[[1]],
                              by = c("variant_id","gene_symbol","rna_expression_fpkm","c_term_pos"),
                              all.x = TRUE)
      
      epitopePredictions=rbindlist(list(epitopePredictions,mergedPredictions))
    }
    
    setTxtProgressBar(progressBar, i)
  }
  close(progressBar)
  
  # sort tables & set new order
  setorderv(x = epitopePredictions,
            cols = paste0("tumor_",hlaType,"affinity"))
  setcolorder(x=epitopePredictions,
              neworder = c("variant_id","gene_symbol",
                           "tumor_peptide","tumor_c_term_aa",paste0("tumor_",hlaType,"affinity"),"tumor_processing_score",
                           "normal_peptide","normal_c_term_aa",paste0("normal_",hlaType,"affinity"),"normal_processing_score","c_term_pos","rna_expression_fpkm"))
  
  setorderv(x = epitopePredictionsWithFiltersApplied,
            cols = paste0("tumor_",hlaTypes,"affinity"))
  setcolorder(x = epitopePredictionsWithFiltersApplied,
              neworder = c("variant_id","gene_symbol",
                           "tumor_peptide","tumor_c_term_aa",paste0("tumor_",hlaType,"affinity"),"tumor_processing_score",
                           "normal_peptide","normal_c_term_aa",paste0("normal_",hlaType,"affinity"),"normal_processing_score","c_term_pos","rna_expression_fpkm"))
  
  # calculate percentile rank
  epitopePredictionsWithFiltersApplied[,percentile_rank:=returnPercentileRank(epitopePredictionsWithFiltersApplied[[paste0("tumor_",hlaType,"affinity")]])]
  
  # write predictions to disk
  if(nrow(epitopePredictionsWithFiltersApplied)>0){
    write.csv(x = unique(x = epitopePredictionsWithFiltersApplied,
                         by = names(epitopePredictionsWithFiltersApplied)[-match(x = c("c_term_pos","variant_id"),
                                                                                 table = names(epitopePredictionsWithFiltersApplied))]),
              file = paste0(dirPath,"/output/",paste(format(Sys.Date(),"%Y%m%d"),sampleId,hlaType,peptideLength,sep="_"),"mer_epitopes.csv"),
              row.names = FALSE)  
  } else {
    write.csv(x = "No epitopes predicted",
              file = paste0(dirPath,"/output/",paste(format(Sys.Date(),"%Y%m%d"),sampleId,hlaType,peptideLength,sep="_"),"mer_epitopes.csv"),
              row.names = FALSE)  
  }
  
  if(nrow(epitopePredictions)>0){
    write.csv(x = unique(x = epitopePredictions,
                         by = names(epitopePredictions)[-match(x = c("c_term_pos","variant_id"),
                                                               table = names(epitopePredictions))]),
              file = paste0(dirPath,"/output/",paste(format(Sys.Date(),"%Y%m%d"),sampleId,hlaType,peptideLength,sep="_"),"mer_epitopes_unfiltered.csv"),
              row.names = FALSE)  
  } else {
    write.csv(x = "No epitopes predicted",
              file = paste0(dirPath,"/output/",paste(format(Sys.Date(),"%Y%m%d"),sampleId,hlaType,peptideLength,sep="_"),"mer_epitopes_unfiltered.csv"),
              row.names = FALSE)  
  }
}