performSingleSequencePredictions=function(file,allele,peptidelength,affcutoff,proccutoff,exprcutoff){
  
}

performPairedSequencePredictions=function(file,allele,peptidelength,affcutoff,proccutoff,exprcutoff){
  # prepare empty table
  # epitopePredictions=data.table()
  # epitopePredictionsWithFiltersApplied=data.table()
  
  # get some info on dataset
  fileName = gsub(pattern = "\\..+$",
                  replacement = "",
                  x = basename(file))
  dirPath = dirname(file)
  
  # import kitchensink data & clean up
  kitchensink=unique(fread(input = file))
  
  sampleId=toupper(gsub(pattern = "_.+$|-.+$",
                        replacement = "",
                        x = fileName))
  variantInfo=returnProcessedVariants(id = sampleId,
                                      variants = kitchensink)
  
  # progressBar = txtProgressBar(min = 0,
  #                              max = nrow(variantInfo),
  #                              width = 100,
  #                              style = 3)
  
  epitopePredictions=foreach(i=1:nrow(variantInfo)) %dopar% {
    # for each variant line, make list tumor peptides which are different from normal (and corresponding normal peptides)
    # and make vector containing normal and tumor peptide stretches
    peptideList=buildPeptideList(variant = variantInfo[i,],
                                 peptidelength = peptidelength)
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
                                 allele = allele,
                                 peptidelength = peptidelength)
    }
    
    # rename columns
    setnames(x = normalAndTumorPredictions[[1]],
             old = c("peptide",paste0(allele,"affinity"),"c_term_aa","processing_score"),
             new = c("normal_peptide",paste0("normal_",allele,"affinity"),"normal_c_term_aa","normal_processing_score"))
    
    setnames(x = normalAndTumorPredictions[[2]],
             old = c("peptide",paste0(allele,"affinity"),"c_term_aa","processing_score"),
             new = c("tumor_peptide",paste0("tumor_",allele,"affinity"),"tumor_c_term_aa","tumor_processing_score"))
    
    # apply various cutoffs
    normalPredictionsWithFiltersApplied=subset(x = normalAndTumorPredictions[[1]],
                                               subset = normalAndTumorPredictions[[1]][[paste0("normal_",allele,"affinity")]] <= affcutoff &
                                                 normal_processing_score >= proccutoff &
                                                 rna_expression_fpkm > exprcutoff)
    
    tumorPredictionsWithFiltersApplied=subset(x = normalAndTumorPredictions[[2]],
                                              subset = normalAndTumorPredictions[[2]][[paste0("tumor_",allele,"affinity")]] <= affcutoff &
                                                tumor_processing_score >= proccutoff &
                                                rna_expression_fpkm > exprcutoff)
    
    # determine self-sim
    
    ### use self-sim which supports 9, 10 and 11mers; get from Marit
    
    # merge all info
    if (nrow(tumorPredictionsWithFiltersApplied)>0){
      mergedTumorPredictionsWithFiltersApplied=merge(x = tumorPredictionsWithFiltersApplied,
                                                     y = normalAndTumorPredictions[[1]],
                                                     by = c("variant_id","gene_symbol","rna_expression_fpkm","c_term_pos"),
                                                     all.x = TRUE)
    } else{
      mergedTumorPredictionsWithFiltersApplied=data.table()
    }
    
    if (nrow(normalAndTumorPredictions[[2]])>0){
      mergedPredictions=merge(x = normalAndTumorPredictions[[2]],
                              y = normalAndTumorPredictions[[1]],
                              by = c("variant_id","gene_symbol","rna_expression_fpkm","c_term_pos"),
                              all.x = TRUE)
    } else{
      mergedPredictions=data.table()
    }
    
    return(list(mergedTumorPredictionsWithFiltersApplied,mergedPredictions))
    
    # setTxtProgressBar(progressBar, i)
  }
  # close(progressBar)
  
  epitopePredictionsTumorWithFiltersApplied=sapply(seq(1,length(epitopePredictions),1), function(x) rbindlist(list(epitopePredictions)))
  
  epitopePredictionsAll=rbindlist(list(epitopePredictions))
  
  # sort tables & set new order
  setorderv(x = epitopePredictions,
            cols = paste0("tumor_",allele,"affinity"))
  setcolorder(x=epitopePredictions,
              neworder = c("variant_id","gene_symbol",
                           "tumor_peptide","tumor_c_term_aa",paste0("tumor_",allele,"affinity"),"tumor_processing_score",
                           "normal_peptide","normal_c_term_aa",paste0("normal_",allele,"affinity"),"normal_processing_score","c_term_pos","rna_expression_fpkm"))
  
  setorderv(x = epitopePredictionsWithFiltersApplied,
            cols = paste0("tumor_",allele,"affinity"))
  setcolorder(x = epitopePredictionsWithFiltersApplied,
              neworder = c("variant_id","gene_symbol",
                           "tumor_peptide","tumor_c_term_aa",paste0("tumor_",allele,"affinity"),"tumor_processing_score",
                           "normal_peptide","normal_c_term_aa",paste0("normal_",allele,"affinity"),"normal_processing_score","c_term_pos","rna_expression_fpkm"))
  
  # calculate percentile rank
  epitopePredictionsWithFiltersApplied[,percentile_rank:=returnPercentileRank(epitopePredictionsWithFiltersApplied[[paste0("tumor_",allele,"affinity")]])]
  
  # write predictions to disk
  if(nrow(epitopePredictionsWithFiltersApplied)>0){
    write.csv(x = unique(x = epitopePredictionsWithFiltersApplied,
                         by = names(epitopePredictionsWithFiltersApplied)[-match(x = c("c_term_pos","variant_id"),
                                                                                 table = names(epitopePredictionsWithFiltersApplied))]),
              file = paste0(dirPath,"/output/",paste(format(Sys.Date(),"%Y%m%d"),sampleId,allele,peptidelength,sep="_"),"mer_epitopes.csv"),
              row.names = FALSE)  
  } else {
    write.csv(x = "No epitopes predicted",
              file = paste0(dirPath,"/output/",paste(format(Sys.Date(),"%Y%m%d"),sampleId,allele,peptidelength,sep="_"),"mer_epitopes.csv"),
              row.names = FALSE)  
  }
  
  if(nrow(epitopePredictions)>0){
    write.csv(x = unique(x = epitopePredictions,
                         by = names(epitopePredictions)[-match(x = c("c_term_pos","variant_id"),
                                                               table = names(epitopePredictions))]),
              file = paste0(dirPath,"/output/",paste(format(Sys.Date(),"%Y%m%d"),sampleId,allele,peptidelength,sep="_"),"mer_epitopes_unfiltered.csv"),
              row.names = FALSE)  
  } else {
    write.csv(x = "No epitopes predicted",
              file = paste0(dirPath,"/output/",paste(format(Sys.Date(),"%Y%m%d"),sampleId,allele,peptidelength,sep="_"),"mer_epitopes_unfiltered.csv"),
              row.names = FALSE)  
  }
}