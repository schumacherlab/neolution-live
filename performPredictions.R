# load external scripts
source("./support_scripts/supportFunctions.R")
source("./support_scripts/runConfig.R")
source("./support_scripts/parseLogic.R")
source("./support_scripts/predictionLogic.R")
source("./support_scripts/selfSimilarityLogic.R")
source("./support_scripts/peptideConstructionLogic.R")
source("./support_scripts/workflowSinglePredictions.R")
source("./support_scripts/workflowPairedPredictions.R")

# register parallel back-end
registerDoMC(cores = numberOfWorkers)

#scriptPath = thisDirectory()

# check availability of predictors
checkPredictorPaths(predictorPaths)

# create directory to hold logs/output, if necessary
dir.create(path = paste0(dirPath, "/output"),
           showWarnings = FALSE)

# collect information on run; print to console and write to log
runStart = format(Sys.time(),"%Y%m%d-%H%M")
message(paste0("Input file:\t\t", filePath, "\n",
               "MHC/HLA allele:\t\t", hlaType, "\n",
               "Peptide length:\t\t", peptideLength, "\n",
               "Affinity cutoff:\t", affinityCutoff, " nM\n",
               "Processsing cutoff:\t", processingCutoff, "\n",
               "Expression cutoff:\t", expressionCutoff, "\n",
               "Single seq predictions:\t\t", doSingleSequencePrediction, "\n",
               "Simple self-similarity:\t\t", doSimpleSelfSimilarity, "\n",
               "Extended self-similarity:\t", doExtendedSelfSimilarity, "\n",
               "Use self-epitope list:\t\t", addSelfEpitopes, "\n"))

# write run info to log
write(x = paste0(Sys.time()," - Neolution run start\n\n",
                 "branch:\t\t\t\t\t", system("git symbolic-ref --short -q HEAD", intern = TRUE),"\n",
                 "commit hash:\t\t", system("git rev-parse HEAD", intern = TRUE),"\n\n",
                 "Input file:\t\t\t\t\t", filePath, "\n",
                 "MHC/HLA allele:\t\t\t", hlaType, "\n",
                 "Peptide length:\t\t\t", peptideLength, "\n",
                 "Affinity cutoff:\t\t", affinityCutoff, " nM\n",
                 "Processing cutoff:\t", processingCutoff, "\n",
                 "Expression cutoff:\t", expressionCutoff, "\n\n",
                 "Single seq predictions:\t\t", doSingleSequencePrediction, "\n",
                 "Simple self-similarity:\t\t", doSimpleSelfSimilarity, "\n",
                 "Extended self-similarity:\t", doExtendedSelfSimilarity, "\n", 
                 "Use self-epitope list:\t\t", addSelfEpitopes, "\n\n",
                 "Affinity predictor:\t\t", predictorPaths$netMHCpan, "\n",
                 "Processing predictor:\t", predictorPaths$netChop, "\n"),
      file = paste0(dirPath, "/output/", paste(runStart, fileName, hlaType, peptideLength, sep = "_"), "mer_runInfo.txt"),
      append = FALSE)

# re-direct output to log file
sink(file = file(description = paste0(dirPath, "/output/", paste(runStart, fileName, hlaType, peptideLength, sep = "_"), "mer_runErrors.txt"),
                 open = "wt"),
     type = "message",
     append = TRUE)

# check if we want single seq predictions or paired normal-tumor predictions
switch(EXPR = as.character(doSingleSequencePrediction),
       # --single == TRUE (single seq input)
       "TRUE" = performSingleSequencePredictions(file = filePath,
                                                 allele = hlaType,
                                                 peptidelength = peptideLength,
                                                 affcutoff = affinityCutoff,
                                                 proccutoff = processingCutoff,
                                                 exprcutoff = expressionCutoff), 
       # --single == FALSE (normal-tumor input)
       "FALSE" = performPairedSequencePredictions(file = filePath,
                                                  allele = hlaType,
                                                  peptidelength = peptideLength,
                                                  affcutoff = affinityCutoff,
                                                  proccutoff = processingCutoff,
                                                  exprcutoff = expressionCutoff)
)

#====================================================================================================================================#
# return output to console
sink()

# write run info to log
write(x = paste0(Sys.time()," - Neolution run end\n\n",
                 "comments:"),
      file = paste0(dirPath, "/output/", paste(runStart, fileName, hlaType, peptideLength, sep = "_"), "mer_runInfo.txt"),
      append = TRUE)
