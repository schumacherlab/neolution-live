processDocOpt <- function(cla = list()) {
  if (!is.null(cla) && length(cla) > 0) {
    cla <- lapply(cla,
      function(x) {
        if (!is.na(x) && class(x) == 'character' && x == 'NA') return(NA)
        else return(x)
      })
  }

  rp <- list()
  rp$file <- cla$file
  rp$allele <- cla$mhc %||% NA
  rp$logfile <- cla$logfile %||% NA
  rp$peptidelength <- suppressWarnings(as.integer(cla$length)) %||% 9
  rp$debug <- suppressWarnings(as.logical(cla$debug)) %||% F
  rp$copyinput <- cla$copyinput %||% F
  rp$copyinput <- suppressWarnings(as.logical(rp$copyinput)) %||% F
  rp$affinity <- suppressWarnings(as.numeric(cla$affinity)) %||% NA
  rp$rank <- suppressWarnings(as.numeric(cla$rank)) %||% 1.9
  rp$ncores <- suppressWarnings(as.integer(cla$ncores)) %||% 1
  rp$run_mode <- cla$run_mode %||% 'paired'
  rp$run_mode <- match.arg(rp$run_mode, 
    choices = c('single', 'structural', 'paired'), several.ok = F)
  rp$model <- suppressWarnings(as.numeric(cla$model)) %||% NA
  rp$use_rfModel <- !is.na(rp$model) %||% F
  rp$processing <- suppressWarnings(as.numeric(cla$processing)) %||% .5
  rp$expression <- suppressWarnings(as.numeric(cla$expression)) %||% 0
  rp$simple_selfsim <- suppressWarnings(as.logical(cla$selfsim)) %||% F
  rp$extended_selfsim <- suppressWarnings(as.logical(cla$extselfsim)) %||% F
  rp$use_selflist <- suppressWarnings(as.logical(cla$selflist)) %||% F
  rp$panversion <- suppressWarnings(as.numeric(cla$panversion)) %||% 3.0
  rp$verbose <- suppressWarnings(as.logical(cla$verbose)) %||% F

  ## List all potential config dirs
  potential_config_locs <- 
    c(getwd(),
      '~/.config',
      file.path('~/libs', 'neolution-live'),
      file.path('~/libs', 'neolution'),
      sapply(.libPaths(), function(p_dir) file.path(p_dir, 'neolution-live')),
      sapply(.libPaths(), function(p_dir) file.path(p_dir, 'neolution')))

  config_fn <-
    sapply(potential_config_locs, list.files, pattern = '.*neolution.*\\.yaml',
      full.names = T) %>%
    .[sapply(., length) == 1] %>%
    .[[1]]

  if (is.null(config_fn) || is.na(config_fn) || length(config_fn) == 0) {
    stop('Could not find a neolution_config.yaml anywhere')
  }

  config <- yaml::read_yaml(config_fn)
  rp$netChop <- path.expand(config[['netChop']])
  rp$netMHCpan <- path.expand(config[['netMHCpan']])
  stopifnot(file.exists(rp$netMHCpan))
  rp$selfEpitopeListPath <- 
    path.expand(config[['selfEpitopeListPath']])
  stopifnot(dir.exists(rp$selfEpitopeListPath))
  rp$randomForestModelPath <- 
    path.expand(config[['randomForestModelPath']])
  stopifnot(file.exists(rp$randomForestModelPath))
  rp$temporaryDirectoryPath <- 
    path.expand(config[['temporaryDirectoryPath']])
  stopifnot(dir.exists(rp$temporaryDirectoryPath))

  if (file.exists(rp$file)) {
    rp$filename <- basename(rp$file)
    rp$filename_no_ext <-
      substring(text = rp$filename,
        first = 1, last = max(unlist(gregexpr(pattern = ".",
            text = rp$filename, fixed = TRUE))) - 1)
    if (is.null(rp$filepath))
      rp$filepath <- dirname(rp$file)
  } else {
    message('Cannot find file, make sure to provide full path to file')
  }

  if (!is.na(rp$allele) && nchar(rp$allele) != 5) {
    message('MHC/HLA input should be formatted as follows: A0201')
  } else {
    rp$allele <- toupper(rp$allele)
  }

  if (rp$peptidelength < 8 || rp$peptidelength > 11) {
    message('Peptide length input (-l or --length) should be >=8 and <= 11')
  }

  if (rp$simple_selfsim && rp$extended_selfsim) {
    message("Please choose ONE type of self-similarity check, use -h for help")
  }

  if (rp$extended_selfsim && rp$peptidelength != 9) {
    message("Extended selfsim can only be used for 9-mers, use -h for help")
  }

  if (rp$run_mode == 'single' &&
      (rp$simple_selfsim || rp$extended_selfsim) &&
      !rp$selflist) {
    message(paste('Self-similarity check on single sequences can only be',
        'performed with a self-epitope list!'))
  }
  return(rp)
}
