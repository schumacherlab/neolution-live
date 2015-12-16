# load external scripts
source("./support_scripts/supportFunctions.R")
source("./support_scripts/runConfig.R")
source("./support_scripts/parseLogic.R")
source("./support_scripts/predictionLogic.R")
source("./support_scripts/selfSimilarityLogic.R")
source("./support_scripts/peptideConstructionLogic.R")
source("./support_scripts/workflowSinglePredictions.R")
source("./support_scripts/workflowPairedPredictions.R")
source("./support_scripts/fasDatabaseLogic.R")

# register parallel back-end
registerDoMC(cores = numberOfWorkers)

#scriptPath = thisDirectory()

# create directory to hold logs/output, if necessary
dir.create(path = paste0(runParameters$filepath, "/output"),
           showWarnings = FALSE)

# collect information on run; print to console and write to log
runStart = format(Sys.time(),"%Y%m%d-%H%M")
message(paste0("Input file:\t\t", paste(runParameters$filepath, runParameters$filename, sep = "/"), "\n",
               "MHC/HLA allele:\t\t", runParameters$allele, "\n",
               "Peptide length:\t\t", runParameters$peptidelength, "\n",
               "Affinity cutoff:\t", runParameters$affinity, " nM\n",
               "Processsing cutoff:\t", runParameters$processing, "\n",
               "Expression cutoff:\t", runParameters$expression, "\n",
               "Single seq predictions:\t\t", runParameters$single_sequence, "\n",
               "Simple self-similarity:\t\t", runParameters$simple_selfsim, "\n",
               "Extended self-similarity:\t", runParameters$extended_selfsim, "\n",
               "Use self-epitope list:\t\t", runParameters$use_selflist, "\n",
               "Use FASdb lookups:\t\t", runParameters$use_fasdb, "\n"))

# write run info to log
write(x = paste0(Sys.time()," - Neolution run start\n\n",
                 "Branch:\t\t\t\t\t", system("git symbolic-ref --short -q HEAD", intern = TRUE),"\n",
                 "Commit hash:\t\t", system("git rev-parse HEAD", intern = TRUE),"\n\n",
                 "Affinity predictor:\t\t", predictorPaths$netMHCpan, "\n",
                 "Processing predictor:\t", predictorPaths$netChop, "\n\n",
                 "Input file:\t\t\t\t\t", paste(runParameters$filepath, runParameters$filename, sep = "/"), "\n",
                 "MHC/HLA allele:\t\t\t", runParameters$allele, "\n",
                 "Peptide length:\t\t\t", runParameters$peptidelength, "\n",
                 "Affinity cutoff:\t\t", runParameters$affinity, " nM\n",
                 "Processing cutoff:\t", runParameters$processing, "\n",
                 "Expression cutoff:\t", runParameters$expression, "\n\n",
                 "Use FASdb lookups:\t\t\t\t", runParameters$use_fasdb, "\n",
                 "Single seq predictions:\t\t", runParameters$single_sequence, "\n",
                 "Simple self-similarity:\t\t", runParameters$simple_selfsim, "\n",
                 "Extended self-similarity:\t", runParameters$extended_selfsim, "\n",
                 "Use self-epitope list:\t\t", runParameters$use_selflist, "\n"),
      file = paste0(runParameters$filepath,
                    "/output/",
                    paste(runStart,
                          runParameters$filename_no_ext,
                          runParameters$allele, 
                          runParameters$peptidelength,
                          sep = "_"),
                    "mer_runInfo.txt"),
      append = FALSE)



# re-direct output to log file
sink(file = file(description = paste0(runParameters$filepath,
                                      "/output/",
                                      paste(runStart,
                                            runParameters$filename_no_ext,
                                            runParameters$allele,
                                            runParameters$peptidelength, sep = "_"),
                                      "mer_runLog.txt"),
                 open = "wt"),
     type = "message",
     append = TRUE)

# check if we want single seq predictions or paired normal-tumor predictions
switch(EXPR = as.character(runParameters$single_sequence),
       # --single == TRUE (single seq input)
       "TRUE" = performSingleSequencePredictions(config = runParameters), 
       # --single == FALSE (normal-tumor input)
       "FALSE" = performPairedSequencePredictions(config = runParameters)
)

#====================================================================================================================================#

# write run info to log
write(x = paste0(Sys.time()," - Neolution run end\n\n",
                 "comments:"),
      file = paste0(runParameters$filepath,
                    "/output/",
                    paste(runStart,
                          runParameters$filename_no_ext,
                          runParameters$allele,
                          runParameters$peptidelength, sep = "_"),
                    "mer_runInfo.txt"),
      append = TRUE)
