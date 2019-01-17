performParallelPredictions <- function(peptides, peptidestretch, allele,
  peptidelength, runParameters) {

  if (runParameters$structural_variants) {
    affinityPredictions <- lapply(peptides,
      function(x) { performAffinityPredictions(
        peptides = x$peptide,
        allele = allele,
        peptidelength = peptidelength,
        runParameters = runParameters)
      })

    # perform processing predictions
    processingPredictions <- lapply(peptidestretch,
      function(x) { performProcessingPredictions(
        peptidestretch = x,
        runParameters = runParameters)
      })

    # merge prediction info
    predictions <- lapply(seq(1, length(affinityPredictions)),
      function(x) {
        if (nrow(peptides[[x]]) < 1) {
          return(data.table())
        }

        dtf <- merge(x = peptides[[x]],
          y = unique(x = affinityPredictions[[x]],
            by = 'peptide'),
          by = 'peptide')

        dtf <- merge(x = dtf,
          y = processingPredictions[[x]],
          by = 'c_term_pos')
        return(dtf)
      })

    predictions <- rbindlist(predictions)
    predictions[, xmer := peptidelength]
  } else {
    if (nrow(peptides) > 0) {
      # perform affinity predictions
      affinityPredictions <- performAffinityPredictions(
        peptides = peptides$peptide,
        allele = allele,
        runParameters = runParameters,
        peptidelength = peptidelength)

      # perform processing predictions
      processingPredictions <- performProcessingPredictions(
        peptidestretch = peptidestretch,
        runParameters = runParameters)

      # merge prediction info
      predictions <- merge(
        x = peptides,
        y = unique(x = affinityPredictions, by = 'peptide'),
        by = 'peptide')

      predictions <- merge(
        x = predictions,
        y = processingPredictions,
        by = 'c_term_pos')

      predictions[, xmer := peptidelength]
    } else {
      predictions <- emptyTableWithColumnNamesAndColumnClasses(
        colnames = c(names(peptides), "xmer", "hla_allele",
          paste0(allele, 'affinity'),
          paste0(allele, 'percentile_rank'),
          'c_term_aa', 'processing_score'),
        colclasses = c(unlist(lapply(peptides, class),
            use.names = FALSE),
          'numeric', 'character', 'numeric', 'numeric', 'character', 'numeric'))
    }
  }
  return(predictions)
}


performAffinityPredictions <- function(
  peptides, allele, peptidelength, runParameters) {

  if (length(peptides) == 0) {
    dtf <- processPredictionOutput('', allele, runParameters = runParameters)
  } else {
    ## Skip binding affinity predictions if allele is NA
    if (is.na(allele)) {
      # hla_allele   peptide A0201affinity A0201percentile_rank
      # HLA-A*02:01 CCGSSSGGS         44368                   85
      dtf <- data.table('hla_allele' = as.character(NA), 'peptide' = peptides,
        'NAaffinity' = as.numeric(NA), 'NApercentile_rank' = as.numeric(NA))
      return(dtf)
    }
    randomString <- sample(c(letters, LETTERS), size = 8, replace = T) %>%
      paste(collapse = '')

    dir.create(file.path(runParameters$temporaryDirectoryPath,
        randomString))

    # write peptides to disk in temp dir
    invisible(sapply(seq(1, length(peptides), 1),
        function(x)
          write(x = sprintf('%s', peptides[x]),
            file = file.path(runParameters$temporaryDirectoryPath,
              randomString, paste0(randomString, '_peps.fas')),
            append = TRUE,
            sep = "\n")))

    # perform predictions
    if (runParameters$panversion >= 3) {
      command <- paste0(# "nice -n 9 ",
        runParameters$netMHCpan,
        ' -a HLA-', gsub(pattern = '^([A-Z]{1}[0-9]{2})([0-9]{2})$',
                         replacement = '\\1:\\2',
                         x = allele),
        ' -l ', peptidelength,
        ' -p ',
        ' -f ', file.path(runParameters$temporaryDirectoryPath,
          randomString, paste0(randomString, '_peps.fas')),
        if (runParameters$panversion >= 4) paste(' -BA'))
    } else {
      command <- paste0(# "nice -n 9 ",
        runParameters$netMHCpan,
        ' -a HLA-', gsub(pattern = '^([A-Z]{1}[0-9]{2})([0-9]{2})$',
                         replacement = '\\1:\\2',
                         x = allele),
        ' -l ', peptidelength,
        ' -p ',
        ' -f ', file.path(runParameters$temporaryDirectoryPath,
          randomString, paste0(randomString, '_peps.fas')),
        ' -tdir ', file.path(runParameters$temporaryDirectoryPath,
          randomString),
        ' -ic50')
    }

    output <- system(command = command, intern = TRUE)

    file.remove(file.path(runParameters$temporaryDirectoryPath,
        randomString, paste0(randomString, '_peps.fas')))
    file.remove(file.path(runParameters$temporaryDirectoryPath,
        randomString))

    dtf <- processPredictionOutput(output, allele,
      runParameters = runParameters)
  }

  dtf <- subset(x = dtf,
    select = colnames(dtf[, -dropNa(match(
          x = c('position', 'variant_id', 'peptide_score_log50k',
            'peptide_core', 'Of', 'Gp', 'Gl', 'Ip', 'Il', 'Icore'),
          table = names(dtf))),
      with = FALSE]))
  return(dtf)
}


processPredictionOutput <- function(raw_output, allele,
  runParameters) {
  ## Regex the raw netMHC output
  raw_output <- raw_output[-grep(pattern = "^\\#.+|^\\-.+|^Protein.+|pos.+|^HLA.+|^$",
                                ignore.case = TRUE,
                                x = raw_output)]
  raw_output <- gsub(pattern = "^[[:blank:]]+| <= WB| <= SB",
                    replacement = "",
                    x = raw_output)
  clean_output <- gsub(pattern = "[[:blank:]]+",
                      replacement = "\t",
                      x = raw_output)

  ## Parse cleaned output
  if (length(clean_output) > 0 & runParameters$panversion >= 3) {
    dtf <- as.data.table(read.table(text = clean_output,
                                    stringsAsFactors = FALSE))
    setnames(x = dtf, old = names(dtf),
      new = c("position", "hla_allele", "peptide", "peptide_core", "Of", "Gp",
        "Gl", "Ip", "Il", "Icore", "variant_id", "peptide_score_log50k",
        paste0(allele, "affinity"), paste0(allele, "percentile_rank")))
  } else if (length(clean_output) > 0) {
    dtf <- as.data.table(read.table(text = clean_output,
        stringsAsFactors = FALSE))
    setnames(x = dtf, old = names(dtf),
      new = c("position", "hla_allele", "peptide", "variant_id",
        "peptide_score_log50k", paste0(allele, "affinity"),
        paste0(allele, "percentile_rank")))
  } else {
    dtf <- emptyTableWithColumnNamesAndColumnClasses(
      colnames = c("position", "hla_allele", "peptide", "variant_id",
        "peptide_score_log50k", paste0(allele, "affinity"),
        paste0(allele, "percentile_rank")),
      colclasses = c("numeric", "character", "character", "numeric", "numeric",
        "numeric", "numeric")
      )
  }
  return(dtf)
}


performProcessingPredictions <- function(
  peptidestretch,
  processing_threshold = runParameters$processing %||% .5,
  temp_root = runParameters$temporaryDirectoryPath %||% tempdir(),
  runParameters = list()) {

  random_string <- c(letters, LETTERS, seq(0, 9)) %>%
    sample(size = 8, replace = T) %>%
    paste(collapse = '')
  tmp_dir <- file.path(temp_root, random_string)
  dir.create(tmp_dir, showWarnings = T)
  stopifnot(dir.exists(tmp_dir))
  system(sprintf('chmod 777 %s', tmp_dir), wait = T)

  fasta_fn <- tempfile(fileext = '.fasta', tmpdir = tmp_dir)
  write(x = sprintf(">1\n%s", peptidestretch),
    file = fasta_fn, append = F, sep = "\n")
  output <- paste(runParameters$netChop, '--threshold',
    processing_threshold, '--method netchop --noplot',
    fasta_fn) %>%
    system(intern = T)
  unlink(tmp_dir, recursive = T)

  if (length(output) > 17) {
    output <- output[-grep(
      pattern = '^\\#.+|^\\-.+|^Number.+|^NetChop.+|pos.+|^$|^sh:|^1:|plot.png',
      x = output)]
    output <- gsub(pattern = "^[[:blank:]]+", replacement = "", x = output)
    output <- gsub(pattern = "[[:blank:]]+", replacement = "\t", x = output)

    ## Parse cleaned output
    dtf <- as.data.table(read.table(text = output, stringsAsFactors = FALSE))
    setnames(x = dtf, old = names(dtf),
      new = c('c_term_pos', 'c_term_aa', 'processing_score', 'extra'))
    dtf[, extra := NULL]
  } else {
    dtf <- emptyTableWithColumnNamesAndColumnClasses(
      colnames = c("c_term_pos", "c_term_aa", "processing_score"),
      colclasses = c("numeric", "character", "numeric"))
  }

  return(dtf)
}
