#### Lookup parameters and other configuration options can be set here ####

# set number of cores available for parallel processing
numberOfWorkers=10

# set predictor paths
predictorPaths=data.table(netMHCpan="/home/NFS/users/l.fanchi/netMHCpan-2.4/netMHCpan",
                          netChop="/home/NFS/users/l.fanchi/netchop-3.1/bin/netChop")

# commandline option specification, do not change
optionList = list(make_option(opt_str = c("-a", "--affinity"),
                              action = "store",
                              type = "double",
                              default = 500,
                              help = "netMHCpan affinity cutoff (optional, default: <= %default nM)"),
                  make_option(opt_str = c("-r", "--rank"),
                              action = "store",
                              type = "integer",
                              default = -1,
                              help = "netMHCpan affinity percentile rank cutoff (optional, default cutoff is nM affinity)"),
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
                  make_option(opt_str = c("-f", "--file"),
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
                              help = "Peptide length (required)")
)

# parse commandline arguments
commandlineArguments=parse_args(OptionParser(option_list=optionList))

# parse other arguments
if (is.null(commandlineArguments$file)) {
  message("File input (-f or --file) is required argument, use -h for help")
  q(status=1)
} else if (file.exists(commandlineArguments$file)) {
  filePath = commandlineArguments$file
  fileName = gsub(pattern = "\\..+$",
                  replacement = "",
                  x = basename(filePath))
  dirPath = dirname(filePath)
} else {
  message("Can't find file, make sure to provide full path to file")
  q(status=1)
}

if (is.null(commandlineArguments$mhc)) {
  message("MHC/HLA input (-m or --mhc) is required argument, use -h for help")
  q(status=1)
} else if (nchar(commandlineArguments$mhc)!=5) {
  message("MHC/HLA input should be formatted as follows: A0201")
  q(status=1)
} else {
  hlaType = toupper(commandlineArguments$mhc)
}

if (is.null(commandlineArguments$length)) {
  message("Peptide length input (-l or --length) is required argument, use -h for help")
  q(status=1)
} else if (commandlineArguments$length >=8 & commandlineArguments$length <=11) {
  peptideLength = commandlineArguments$length
} else {
  message("Peptide length input (-l or --length) should be >=8 and <= 11, use -h for help")
  q(status=1)
}

if (is.numeric(commandlineArguments$affinity)){
  affinityCutoff = commandlineArguments$affinity
} else {
  message("Affinity cutoff input (-a or --affinity) should be numeric, use -h for help")
  q(status=1)
}

if (is.numeric(commandlineArguments$processing)){
  processingCutoff = commandlineArguments$processing
} else {
  message("Processing cutoff input (-p or --processing) should be numeric, use -h for help")
  q(status=1)
}

if (is.numeric(commandlineArguments$expression)){
  expressionCutoff = commandlineArguments$expression
} else {
  message("Expression cutoff input (-e or --expression) should be numeric, use -h for help")
  q(status=1)
}

