# load external scripts
source("./runConfig.R")
source("./support_scripts/parseLogic.R")
source("./support_scripts/predictionLogic.R")
source("./support_scripts/supportFunctions.R")
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
               "Chop score cutoff:\t",chopCutoff,"\n",
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
                 "branch:\t\t\t",system("git symbolic-ref --short -q HEAD",intern = TRUE),"\n",
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

for(i in 1:nrow(variantInfo)){
  # for each variant line, make list tumor peptides which are different from normal (and corresponding normal peptides)
  peptideList=buildPeptideList(variant = variantInfo[i,],
                               peptidelength = peptideLength)
  
  # if no tumor peptides found, move to next line
  if(nrow(peptideList[[2]])<1){
    next
  }
  
  if(nrow(peptideList[[1]])>0){
    # perform affinity predictions
    normalAffinityPredictions=performAffinityPredictions(peptides = peptideList[[1]]$normal_peptide,
                                                         allele = hlaType,
                                                         peptidelength = peptideLength)
    setnames(x = normalAffinityPredictions,
             old = c("peptide",paste0(hlaType,"affinity")),
             new = c("normal_peptide",paste0("normal_",hlaType,"affinity")))
    
    # perform processing predictions
    normalProcessingPredictions=performProcessingPredictions(peptidestretch = variantInfo[i,]$peptidecontextnormal)
    
    setnames(x = normalProcessingPredictions,
             old = c("c_term_pos","c_term_aa","processing_score"),
             new = c("c_term_pos","normal_c_term_aa","normal_processing_score"))
    
    # merge prediction info
    normalPredictions=merge(x = peptideList[[1]],
                            y = normalAffinityPredictions,
                            by = "normal_peptide")
    
    normalPredictions=merge(x = normalPredictions,
                            y = normalProcessingPredictions,
                            by = "c_term_pos")
  } else {
    normalPredictions=as.data.table(setNames(replicate(n = 8,
                                                       expr = numeric(0),
                                                       simplify = FALSE),
                                             c("c_term_pos","normal_peptide","variant_id","gene_symbol","rna_expression_fpkm",
                                               paste0("normal_",hlaType,"affinity"),"normal_c_term_aa","normal_processing_score")))
  }
  
  if(nrow(peptideList[[2]])>0){
    # perform affinity predictions
    tumorAffinityPredictions=performAffinityPredictions(peptides = peptideList[[2]]$tumor_peptide,
                                                        allele = hlaType,
                                                        peptidelength = peptideLength)
    setnames(x = tumorAffinityPredictions,
             old = c("peptide",paste0(hlaType,"affinity")),
             new = c("tumor_peptide",paste0("tumor_",hlaType,"affinity")))
    
    # perform processing predictions
    tumorProcessingPredictions=performProcessingPredictions(peptidestretch = variantInfo[i,]$peptidecontexttumor)
    
    setnames(x = tumorProcessingPredictions,
             old = c("c_term_pos","c_term_aa","processing_score"),
             new = c("c_term_pos","tumor_c_term_aa","tumor_processing_score"))
    
    # merge prediction info
    tumorPredictions=merge(x = peptideList[[2]],
                           y = tumorAffinityPredictions,
                           by = "tumor_peptide")
    
    tumorPredictions=merge(x = tumorPredictions, 
                           y = tumorProcessingPredictions, 
                           by = "c_term_pos")
  }
  
  # apply various cutoffs
  normalPredictionsWithFiltersApplied=subset(x = normalPredictions,
                                             subset = paste0("normal_",hlaType,"affinity") <= affinityCutoff &
                                               normal_processing_score >= processingCutoff &
                                               rna_expression_fpkm > expressionCutoff)
  
  tumorPredictionsWithFiltersApplied=subset(x = tumorPredictions,
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
    
    epitopePredictions=rbindlist(list(epitopePredictions,mergedPredictions))
  }
}

write.csv(x = epitopePredictions,
          file = paste0(dirPath,"/",paste(format(Sys.Date(),"%Y%m%d"),sampleId,hlaType,peptideLength,sep="_"),"mer_epitopes.csv"),
          row.names = FALSE)

# write run info to log
write(x = paste0(Sys.time()," - Neolution run end\n\n",
                 "comments:"),
      file = paste0(dirPath,"/output/",runStart,"_",fileName,"_runInfo.txt"),
      append = TRUE)
