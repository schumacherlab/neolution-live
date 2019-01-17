regexPatterns <- list(
  file_extension = '\\.[^.]+$', # match file extension (everything after last dot, inclusive)
  snp_identifier = '[gr]s\\d+', # for matching SNPs
  gs_identifier = 'gs\\d+', # for matching snps not found in dbSNP, keep boundless (no '^' or '$')
  rs_identifier = 'rs\\d+', # for matching snps found in dbSNP, keep boundless (no '^' or '$')
  cosmic_identifier = 'COSM\\d+', # for matching variants found in COSMIC coding muts database, keep boundless
  seqdata_prefix = '(^.+[ATCG]{7,8})(.+)', # for isolating GCF prefix
  allele_exclusion = 'C[0-9]{4}') # for excluding particular alleles from analysis


processDocOpt <- function(commandlineArguments) {
  runParameters <- list()
  runParameters$file <- commandlineArguments$file
  ## TODO make output file configurable
  runParameters$allele <- toupper(commandlineArguments$mhc)
  runParameters$logfile <- commandlineArguments$logfile
  runParameters$peptidelength <-
    suppressWarnings(as.integer(commandlineArguments$length))
  runParameters$debug <-
    suppressWarnings(as.logical(commandlineArguments$debug))
  runParameters$copyinput <- commandlineArguments$copyinput %||% F
  runParameters$copyinput <-
    suppressWarnings(as.logical(runParameters$copyinput))
  runParameters$affinity <-
    suppressWarnings(as.numeric(commandlineArguments$affinity))
  runParameters$rank <-
    suppressWarnings(as.numeric(commandlineArguments$rank))
  runParameters$ncores <-
    suppressWarnings(as.integer(commandlineArguments$ncores))
  runParameters$run_mode <- match.arg(commandlineArguments$run_mode, 
    choices = c('single', 'structural', 'paired'), several.ok = F)
  runParameters$model <-
    suppressWarnings(as.numeric(commandlineArguments$model))
  runParameters$use_rfModel <- !is.na(runParameters$model)
  runParameters$processing <-
    suppressWarnings(as.numeric(commandlineArguments$processing))
  runParameters$expression <-
    suppressWarnings(as.numeric(commandlineArguments$expression))
  runParameters$simple_selfsim <-
    suppressWarnings(as.logical(commandlineArguments$selfsim))
  runParameters$extended_selfsim <-
    suppressWarnings(as.logical(commandlineArguments$extselfsim))
  runParameters$use_selflist <-
    suppressWarnings(as.logical(commandlineArguments$selflist))
  runParameters$panversion <-
    suppressWarnings(as.numeric(commandlineArguments$panversion))
  runParameters$verbose <-
    suppressWarnings(as.logical(commandlineArguments$verbose))

  ## List all potential config dirs
  potential_config_locs <- 
    c(getwd(),
      '~/.config',
      file.path('~/libs', 'neolution-live'),
      file.path('~/libs', 'neolution'),
      sapply(.libPaths(), function(p_dir) file.path(p_dir, 'neolution-live')),
      sapply(.libPaths(), function(p_dir) file.path(p_dir, 'neolution')))

  config_fn <-
    sapply(potential_config_locs, list.files, pattern = 'neolution_config.yaml',
      full.names = T) %>%
    .[sapply(., length) == 1] %>%
    .[[1]]

  if (is.null(config_fn) || is.na(config_fn) || length(config_fn) == 0) {
    stop('Could not find a neolution_config.yaml anywhere')
  }

  config <- yaml::read_yaml(config_fn)
  runParameters$netChop <- path.expand(config[['netChop']])
  runParameters$netMHCpan <- path.expand(config[['netMHCpan']])
  stopifnot(file.exists(runParameters$netMHCpan))
  runParameters$selfEpitopeListPath <- 
    path.expand(config[['selfEpitopeListPath']])
  stopifnot(dir.exists(runParameters$selfEpitopeListPath))
  runParameters$randomForestModelPath <- 
    path.expand(config[['randomForestModelPath']])
  stopifnot(file.exists(runParameters$randomForestModelPath))
  runParameters$temporaryDirectoryPath <- 
    path.expand(config[['temporaryDirectoryPath']])
  stopifnot(dir.exists(runParameters$temporaryDirectoryPath))

  runParameters <- lapply(runParameters,
    function(x) {
      if (!is.na(x) && class(x) == 'character' && x == 'NA') return(NA)
      else return(x)
    })

  if (file.exists(runParameters$file)) {
    runParameters$filename <- basename(runParameters$file)
    runParameters$filename_no_ext <-
      substring(text = runParameters$filename,
        first = 1, last = max(unlist(gregexpr(pattern = ".",
            text = runParameters$filename, fixed = TRUE))) - 1)
    if (is.null(runParameters$filepath))
      runParameters$filepath <- dirname(runParameters$file)
  } else {
    message('Cannot find file, make sure to provide full path to file')
  }

  if (!is.na(runParameters$allele) && nchar(runParameters$allele) != 5) {
    message('MHC/HLA input should be formatted as follows: A0201')
  }

  if (runParameters$peptidelength < 8 || runParameters$peptidelength > 11) {
    message('Peptide length input (-l or --length) should be >=8 and <= 11')
  }

  if (runParameters$simple_selfsim && runParameters$extended_selfsim) {
    message("Please choose ONE type of self-similarity check, use -h for help")
  }

  if (runParameters$extended_selfsim && runParameters$peptidelength != 9) {
    message("Extended selfsim can only be used for 9-mers, use -h for help")
  }

  if (runParameters$run_mode == 'single' &&
      (runParameters$simple_selfsim || runParameters$extended_selfsim) &&
      !runParameters$selflist) {
    message(paste('Self-similarity check on single sequences can only be',
        'performed with a self-epitope list!'))
  }
  return(runParameters)
}
