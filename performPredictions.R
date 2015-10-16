suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(gtools))
suppressPackageStartupMessages(library(utils))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(doMC))
suppressPackageStartupMessages(library(optparse))

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

# # prepare some empty lists
# normalEpitopes_HLA=vector("list",length(hlaTypes))
# normalEpitopes_filtered=vector("list",length(hlaTypes))
# tumorEpitopes_HLA=vector("list",length(hlaTypes))
# tumorEpitopes_filtered=vector("list",length(hlaTypes))

# check availability of predictors
checkPredictorPaths(predictorPaths)

# create directory to hold logs/output, if necessary
dir.create(path = paste0(dirPath,"/tmp"),
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
  peptideList=buildPeptideList(variant = variantInfo[i],
                               peptidelength = peptideLength)
  
  # if no tumor peptides found, move to next line
  if(nrow(peptideList[[2]])==0){
    next
  }
  
  if(nrow(peptideList[[1]])>0){
    invisible(sapply(seq(1,nrow(peptideList[[1]]),1), function(x)
      write(x=sprintf(">%i\n%s",x,peptideList[[1]]$normal_peptides[x]),
            file=paste0(dirPath,"/tmp/normal_peps.fas"),
            append=TRUE,
            sep="\n")))
    
    normalPredictions=performAffinityPredictions(fastafile = paste0(dirPath,"/tmp/normal_peps.fas"),
                                                 allele = hlaType,
                                                 peptidelength = peptideLength)
    
    file.remove(paste0(dirPath,"/tmp/normal_peps.fas"))
  }
  
  if(nrow(peptideList[[2]])>0){
    invisible(sapply(seq(1,nrow(peptideList[[2]]),1), function(x)
      write(x=sprintf(">%i\n%s",x,peptideList[[2]]$tumor_peptides[x]),
            file=paste0(dirPath,"/tmp/tumor_peps.fas"),
            append=TRUE,
            sep="\n")))
    
    tumorPredictions=performAffinityPredictions(fastafile = paste0(dirPath,"/tmp/tumor_peps.fas"),
                                                 allele = hlaType,
                                                 peptidelength = peptideLength)
    
    file.remove(paste0(dirPath,"/tmp/tumor_peps.fas"))
  }
  
}

# write run info to log
write(x = paste0(Sys.time()," - Neolution run end\n\n",
                 "comments:"),
      file = paste0(dirPath,"/output/",runStart,"_",fileName,"_runInfo.txt"),
      append = TRUE)
