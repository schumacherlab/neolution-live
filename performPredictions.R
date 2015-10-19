# load external scripts
source("./support_scripts/supportFunctions.R")
source("./support_scripts/runConfig.R")
source("./support_scripts/parseLogic.R")
source("./support_scripts/predictionLogic.R")
source("./support_scripts/peptideConstructionLogic.R")

# register parallel back-end
registerDoMC(numberOfWorkers)

#scriptPath=thisDirectory()
#setwd(scriptPath)

runStart=format(Sys.Date(),"%Y%m%d")

message(paste0("Input file:\t\t",filePath,"\n",
               "MHC/HLA allele:\t\t",hlaType,"\n",
               "Peptide length:\t\t",peptideLength,"\n",
               "Affinity cutoff:\t",affinityCutoff," nM\n",
               "Chop score cutoff:\t",processingCutoff,"\n",
               "Expression cutoff:\t",expressionCutoff))

# prepare empty list
epitopePredictions=data.table()

# check availability of predictors
checkPredictorPaths(predictorPaths)

# create directory to hold logs/output, if necessary
dir.create(path = "./tmp",
           showWarnings=FALSE)
dir.create(path = paste0(dirPath,"/output"),
           showWarnings=FALSE)

# write run info to log
write(x = paste0(Sys.time()," - Neolution run start\n\n",
                 "branch:\t\t\t\t",system("git symbolic-ref --short -q HEAD",intern = TRUE),"\n",
                 "commit hash:\t\t",system("git rev-parse HEAD",intern = TRUE),"\n"),
      file = paste0(dirPath,"/output/",runStart,"_",fileName,"_runInfo.txt"),
      append = FALSE)

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
  normalAndTumorPredictions=foreach(i=1:2) %dopar% {
    performParallelPredictions(peptides = peptideList[[i]],
                               peptidestretch = peptideStretchVector[i],
                               allele = hlaType,
                               peptidelength = peptideLength)
  }
  
  setnames(x = normalAndTumorPredictions[[1]],
           old = c("peptide",paste0(hlaType,"affinity"),"c_term_aa","processing_score"),
           new = c("normal_peptide",paste0("normal_",hlaType,"affinity"),"normal_c_term_aa","normal_processing_score"))
  
  setnames(x = normalAndTumorPredictions[[2]],
           old = c("peptide",paste0(hlaType,"affinity"),"c_term_aa","processing_score"),
           new = c("tumor_peptide",paste0("tumor_",hlaType,"affinity"),"tumor_c_term_aa","tumor_processing_score"))
  
  # apply various cutoffs
  normalPredictionsWithFiltersApplied=subset(x = normalAndTumorPredictions[[1]],
                                             subset = paste0("normal_",hlaType,"affinity") <= affinityCutoff &
                                               normal_processing_score >= processingCutoff &
                                               rna_expression_fpkm > expressionCutoff)
  
  tumorPredictionsWithFiltersApplied=subset(x = normalAndTumorPredictions[[2]],
                                            subset = paste0("tumor_",hlaType,"affinity") <= affinityCutoff &
                                              tumor_processing_score >= processingCutoff &
                                              rna_expression_fpkm > expressionCutoff)
  
  # determine self-sim, only if peptide length == 9
  
  ### do we want to do self-sim only on 9-mers? or not do for any? ###
  
  # merge all info
  if (nrow(tumorPredictionsWithFiltersApplied)>0){
    mergedPredictions=merge(x = tumorPredictionsWithFiltersApplied,
                            y = normalPredictionsWithFiltersApplied,
                            by = c("variant_id","gene_symbol","rna_expression_fpkm","c_term_pos"),
                            all.x = TRUE)
    
    print(mergedPredictions)
    
    epitopePredictions=rbindlist(list(epitopePredictions,mergedPredictions))
  }
  
  setTxtProgressBar(progressBar, i)
}

close(progressBar)

file.remove("./tmp")

if(nrow(epitopePredictions)>0){
  write.csv(x = epitopePredictions,
            file = paste0(dirPath,"/output/",paste(format(Sys.Date(),"%Y%m%d"),sampleId,hlaType,peptideLength,sep="_"),"mer_epitopes.csv"),
            row.names = FALSE)  
} else {
  write.csv(x = "No epitopes predicted",
            file = paste0(dirPath,"/output/",paste(format(Sys.Date(),"%Y%m%d"),sampleId,hlaType,peptideLength,sep="_"),"mer_epitopes.csv"),
            row.names = FALSE)  
}


# write run info to log
write(x = paste0(Sys.time()," - Neolution run end\n\n",
                 "comments:"),
      file = paste0(dirPath,"/output/",runStart,"_",fileName,"_runInfo.txt"),
      append = TRUE)
