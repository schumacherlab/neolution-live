# # commandline option specification, do not change
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
                  make_option(opt_str = c("-r", "--rank"),
                              action = "store",
                              type = "double",
                              default = NULL,
                              help = "netMHCpan percentile rank cutoff (optional, default cutoff is affinity)"),
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
                  make_option(opt_str = c("--structural"),
                              action = "store_true",
                              type = "logical",
                              default = FALSE,
                              help = "Structural variant predictions (optional, default: %default)"),
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
                              help = "Look up peptide affinity in FASdb, predict if not found; only compatible with 9-mers & netMHCpan-2.4 (optional, default: %default)"),
                  make_option(opt_str = c("--panversion"),
                              action = "store",
                              type = "double",
                              default = 3.0,
                              help = "netMHCpan version (optional, default: %default)")
)

# parse commandline arguments
commandlineArguments = parse_args(OptionParser(option_list = optionList))

# prepare table for holding run configuration
runParameters = vector(mode = "list", length = 16)
runParameters = setNames(object = runParameters, nm = c("filename", "filename_no_ext", "filepath",
                                                        "allele", "peptidelength", "affinity", "rank", "processing", "expression",
                                                        "single_sequence", "structural_variants", "simple_selfsim", "extended_selfsim", "use_selflist", "use_fasdb", "panversion"))

# parse other arguments
if (is.null(commandlineArguments$file)) {
  message("File input (-f or --file) is required argument, use -h for help")
  q(status = 1)
} else if (file.exists(commandlineArguments$file)) {
  runParameters$filename = basename(commandlineArguments$file)
  runParameters$filename_no_ext = substring(text = runParameters$filename,
                                            first = 1,
                                            last = max(unlist(gregexpr(pattern = ".",
                                                                       text = runParameters$filename,
                                                                       fixed = TRUE))) - 1)
  runParameters$filepath = dirname(commandlineArguments$file)
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
  runParameters$allele = toupper(commandlineArguments$mhc)
}

if (is.null(commandlineArguments$length)) {
  message("Peptide length input (-l or --length) is required argument, use -h for help")
  q(status = 1)
} else if (commandlineArguments$length >= 8 & commandlineArguments$length <= 11) {
  runParameters$peptidelength = commandlineArguments$length
} else {
  message("Peptide length input (-l or --length) should be >=8 and <= 11, use -h for help")
  q(status = 1)
}

if (is.numeric(commandlineArguments$affinity)) {
  runParameters$affinity = commandlineArguments$affinity
} else {
  message("Affinity cutoff input (-a or --affinity) should be numeric, use -h for help")
  q(status = 1)
}

if (is.numeric(commandlineArguments$rank) | is.null(commandlineArguments$rank)) {
  runParameters$rank = commandlineArguments$rank
} else {
  message("Rank cutoff input (-r or --rank) should be numeric, use -h for help")
  q(status = 1)
}

if (is.numeric(commandlineArguments$processing)) {
  runParameters$processing = commandlineArguments$processing
} else {
  message("Processing cutoff input (-p or --processing) should be numeric, use -h for help")
  q(status = 1)
}

if (is.numeric(commandlineArguments$expression)) {
  runParameters$expression = commandlineArguments$expression
} else {
  message("Expression cutoff input (-e or --expression) should be numeric, use -h for help")
  q(status = 1)
}

runParameters$single_sequence = commandlineArguments$single

runParameters$structural_variants = commandlineArguments$structural

if (commandlineArguments$selfsim & commandlineArguments$extselfsim) {
  message("Please choose ONE type of self-similarity check, use -h for help")
  q(status = 1)
} else if (commandlineArguments$extselfsim & runParameters$peptidelength != 9) {
  message("Extended selfsim can only be used for 9-mers, use -h for help")
  q(status = 1)
} else {
  runParameters$simple_selfsim = commandlineArguments$selfsim
  runParameters$extended_selfsim = commandlineArguments$extselfsim
}

if (runParameters$single_sequence & (runParameters$simple_selfsim | runParameters$extended_selfsim) & commandlineArguments$selflist == FALSE) {
  message("Self-similarity check on single sequences can only be performed with a self-epitope list!")
  q(status = 1)
} else {
  runParameters$use_selflist = commandlineArguments$selflist
}

if (is.numeric(commandlineArguments$panversion)) {
  runParameters$panversion = format(round(commandlineArguments$panversion, 2), nsmall = 1)
} else {
  message("Please input version number for netMHCpan (e.g. 2.4, 2.8 or 3.0")
  q(status = 1)
}

# if (commandlineArguments$fasdb & runParameters$panversion == "2.4" & runParameters$peptidelength == 9) {
#   runParameters$use_fasdb = commandlineArguments$fasdb
# } else if (commandlineArguments$fasdb & (runParameters$panversion != "2.4" | runParameters$peptidelength != 9)) {
#   message("FASdb peptide affinity lookups are only supported for 9-mers & when using netMHCpan-2.4, use -h for help")
#   q(status = 1)
# } else {
#   runParameters$use_fasdb = commandlineArguments$fasdb
# }

if (commandlineArguments$fasdb & runParameters$panversion == "2.4" & is.numeric(runParameters$rank)) {
  message("Using FASdb with netMHCpan-2.4 and rank cutoffs is not supported")
  q(status = 1)
} else {
  runParameters$use_fasdb = commandlineArguments$fasdb
}
