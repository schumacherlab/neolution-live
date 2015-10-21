# load external scripts
source("./support_scripts/supportFunctions.R")
source("./support_scripts/runConfig.R")
source("./support_scripts/parseLogic.R")
source("./support_scripts/workflowLogic.R")
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
               "Processsing cutoff:\t",processingCutoff,"\n",
               "Expression cutoff:\t",expressionCutoff))

# prepare empty table
epitopePredictions=data.table()
epitopePredictionsWithFiltersApplied=data.table()

# check availability of predictors
checkPredictorPaths(predictorPaths)

# create directory to hold logs/output, if necessary
dir.create(path = "./tmp",
           showWarnings=FALSE)
dir.create(path = paste0(dirPath,"/output"),
           showWarnings=FALSE)

# write run info to log
write(x = paste0(Sys.time()," - Neolution run start\n\n",
                 "branch:\t\t\t\t\t",system("git symbolic-ref --short -q HEAD",intern = TRUE),"\n",
                 "commit hash:\t\t",system("git rev-parse HEAD",intern = TRUE),"\n\n",
                 "Input file:\t\t\t\t\t",filePath,"\n",
                 "MHC/HLA allele:\t\t\t",hlaType,"\n",
                 "Peptide length:\t\t\t",peptideLength,"\n",
                 "Affinity cutoff:\t\t",affinityCutoff," nM\n",
                 "Processing cutoff:\t",processingCutoff,"\n",
                 "Expression cutoff:\t",expressionCutoff,"\n\n",
                 "Affinity predictor:\t\t",predictorPaths$netMHCpan,"\n",
                 "Processing predictor:\t",predictorPaths$netChop),
      file = paste0(dirPath,"/output/",paste(runStart,fileName,hlaType,peptideLength,sep="_"),"mer_runInfo.txt"),
      append = FALSE)


# start branching here (use switch)

# -s == TRUE (single seq input)

# -s == FALSE (normal-tumor input)
performPairedSequencePredictions(file = filePath,
                                 allele = hlaType,
                                 peptidelength = peptideLength,
                                 affcutoff = affinityCutoff,
                                 proccutoff = processingCutoff,
                                 exprcutoff = expressionCutoff)

#====================================================================================================================================#

# remove temp dir
file.remove("./tmp")

# write run info to log
write(x = paste0(Sys.time()," - Neolution run end\n\n",
                 "comments:"),
      file = paste0(dirPath,"/output/",paste(runStart,fileName,hlaType,peptideLength,sep="_"),"mer_runInfo.txt"),
      append = TRUE)
