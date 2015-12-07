#### Lookup parameters and other configuration options can be set here ####

# set number of cores available for parallel processing
numberOfWorkers = 32

# set predictor paths
predictorPaths = data.table(netMHCpan = "/home/NFS/users/l.fanchi/netMHCpan-2.4/netMHCpan",
                            netChop = "/home/NFS/users/l.fanchi/netchop-3.1/bin/netChop")

# set path of self-epitope lists
selfEpitopeListPath = "./selflists_v2.4"

# set path of temporary directory (set chmod 1777 on tempdir)
temporaryDirectoryPath = "~/scratch"

### ============================================================================================================================================##

# commandline option specification, do not change
optionList = list(make_option(opt_str = c("-f", "--file"),
                              action = "store",
                              type = "character",
                              default = NULL,
                              help = "Full path to file containing variant calls (required)"),
                  make_option(opt_str = c("-m", "--mhc"),
                              action = "store",
                              type = "character",
                              default = NULL,
                              help = "MHC/HLA allele, formatted as follows: A0201 (required)"),
                  make_option(opt_str = c("-l", "--length"),
                              action = "store",
                              type = "integer",
                              default = NULL,
                              help = "Peptide length (required)"),
                  make_option(opt_str = c("-a", "--affinity"),
                              action = "store",
                              type = "double",
                              default = 500,
                              help = "netMHCpan affinity cutoff (optional, default: <= %default nM)"),
                  # make_option(opt_str = c("-r", "--rank"),
                  #             action = "store",
                  #             type = "integer",
                  #             default = -1,
                  #             help = "netMHCpan affinity percentile rank cutoff (optional, default cutoff is nM affinity)"),
                  make_option(opt_str = c("-p", "--processing"),
                              action = "store",
                              type = "double",
                              default = 0.5,
                              help = "netChop processing score cutoff (optional, default: >= %default)"),
                  make_option(opt_str = c("-e", "--expression"),
                              action = "store",
                              type = "double",
                              default = 0,
                              help = "RNA expression cutoff (optional, default: > %default)"),
                  make_option(opt_str = c("--single"),
                              action = "store_true",
                              type = "logical",
                              default = FALSE,
                              help = "Single sequence predictions (not paired normal-tumor) (optional, default: %default)"),
                  make_option(opt_str = c("--selfsim"),
                              action = "store_true",
                              type = "logical",
                              default = FALSE,
                              help = "Perform simple self-similarity check; compatible with 9-, 10-, 11-mers (optional, default: %default)"),
                  make_option(opt_str = c("--extselfsim"),
                              action = "store_true",
                              type = "logical",
                              default = FALSE,
                              help = "Perform extended self-similarity check; only compatible with 9-mers (optional, default: %default)"),
                  make_option(opt_str = c("--selflist"),
                              action = "store_true",
                              type = "logical",
                              default = FALSE,
                              help = "Add predicted self-epitopes to self-similarity check, requires length- & HLA-matched selflist (optional, default: %default)"),
                  make_option(opt_str = c("--fasdb"),
                              action = "store_true",
                              type = "logical",
                              default = FALSE,
                              help = "Look up peptide affinity in FASdb, predict if not found; only compatible with 9-mers & netMHCpan-2.4 (optional, default: %default)")
)

# check availability of predictors
checkPredictorPaths(predictorPaths)

# find out if netMHCpan version is 2.4 (bit hackish)
isFasDbVersion = any(grepl(pattern = "netmhcpan 2\\.4",
                           x = readLines(con = predictorPaths$netMHCpan),
                           ignore.case = TRUE))

# parse commandline arguments
commandlineArguments = parse_args(OptionParser(option_list = optionList))

# parse other arguments
if (is.null(commandlineArguments$file)) {
  message("File input (-f or --file) is required argument, use -h for help")
  q(status = 1)
} else if (file.exists(commandlineArguments$file)) {
  filePath = commandlineArguments$file
  fileName = substring(text = basename(filePath),
                       first =  1, 
                       last = max(unlist(gregexpr(pattern = ".", 
                                                  text = basename(filePath),
                                                  fixed = TRUE)))-1
  )
  dirPath = dirname(filePath)
} else {
  message("Can't find file, make sure to provide full path to file")
  q(status = 1)
}

if (is.null(commandlineArguments$mhc)) {
  message("MHC/HLA input (-m or --mhc) is required argument, use -h for help")
  q(status = 1)
} else if (nchar(commandlineArguments$mhc) != 5) {
  message("MHC/HLA input should be formatted as follows: A0201")
  q(status = 1)
} else {
  hlaType = toupper(commandlineArguments$mhc)
}

if (is.null(commandlineArguments$length)) {
  message("Peptide length input (-l or --length) is required argument, use -h for help")
  q(status = 1)
} else if (commandlineArguments$length >=8 & commandlineArguments$length <= 11) {
  peptideLength = commandlineArguments$length
} else {
  message("Peptide length input (-l or --length) should be >=8 and <= 11, use -h for help")
  q(status = 1)
}

if (is.numeric(commandlineArguments$affinity)) {
  affinityCutoff = commandlineArguments$affinity
} else {
  message("Affinity cutoff input (-a or --affinity) should be numeric, use -h for help")
  q(status = 1)
}

if (is.numeric(commandlineArguments$processing)) {
  processingCutoff = commandlineArguments$processing
} else {
  message("Processing cutoff input (-p or --processing) should be numeric, use -h for help")
  q(status = 1)
}

if (is.numeric(commandlineArguments$expression)) {
  expressionCutoff = commandlineArguments$expression
} else {
  message("Expression cutoff input (-e or --expression) should be numeric, use -h for help")
  q(status = 1)
}

doSingleSequencePrediction = commandlineArguments$single

if (commandlineArguments$selfsim & commandlineArguments$extselfsim) {
  message("Please choose ONE type of self-similarity check, use -h for help")
  q(status = 1)
} else if (commandlineArguments$extselfsim & peptideLength != 9) {
  message("Extended selfsim can only be used for 9-mers, use -h for help")
  q(status = 1)
} else {
  doSimpleSelfSimilarity = commandlineArguments$selfsim
  doExtendedSelfSimilarity = commandlineArguments$extselfsim
}

if (doSingleSequencePrediction & (doSimpleSelfSimilarity | doExtendedSelfSimilarity) & commandlineArguments$selflist == FALSE) {
  message("Self-similarity check on single sequences can only be performed with a self-epitope list!")
  q(status = 1)
} else {
  addSelfEpitopes = commandlineArguments$selflist
}

if (commandlineArguments$fasdb & isFasDbVersion & peptideLength == 9) {
  useFasDb = commandlineArguments$fasdb
} else if (commandlineArguments$fasdb & (isFasDbVersion == FALSE | peptideLength != 9)) {
  message("FASdb peptide affinity lookups are only supported for 9-mers & when using netMHCpan-2.4, use -h for help")
  q(status = 1)
} else {
  useFasDb = commandlineArguments$fasdb
}
