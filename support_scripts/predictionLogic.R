## this support script contains functions for peptide affinity and processing score predictions

performParallelPredictions = function(peptides, peptidestretch, allele, peptidelength) {
  if (runParameters$structural_variants) {
    # perform affinity predictions
    affinityPredictions = lapply(peptides,
                                 function(x) {
                                   performAffinityPredictions(peptides = x$peptide,
                                                              allele = allele,
                                                              peptidelength = peptidelength)
                                 })

    # perform processing predictions
    processingPredictions = lapply(peptidestretch,
                                   function(x) {
                                     performProcessingPredictions(peptidestretch = x)
                                   })

    # merge prediction info
    predictions = lapply(seq(1, length(affinityPredictions)),
                         function(x) {
                           if (nrow(peptides[[x]]) < 1) {
                             return(data.table())
                           }
                           data = merge(x = peptides[[x]],
                                        y = unique(x = affinityPredictions[[x]],
                                                   by = "peptide"),
                                        by = "peptide")

                           data = merge(x = data,
                                        y = processingPredictions[[x]],
                                        by = "c_term_pos")
                           return(data)
                         })

    predictions = rbindlist(predictions)

    predictions[, xmer := peptidelength]
  } else {
    if (nrow(peptides) > 0) {
      # perform affinity predictions
      affinityPredictions = performAffinityPredictions(peptides = peptides$peptide,
                                                       allele = allele,
                                                       peptidelength = peptidelength)

      # perform processing predictions
      processingPredictions = performProcessingPredictions(peptidestretch = peptidestretch)

      # merge prediction info
      predictions = merge(x = peptides,
                          y = unique(x = affinityPredictions,
                                     by = "peptide"),
                          by = "peptide")

      predictions = merge(x = predictions,
                          y = processingPredictions,
                          by = "c_term_pos")

      predictions[, xmer := peptidelength]
    } else {
      predictions = emptyTableWithColumnNamesAndColumnClasses(colnames = c(names(peptides),
                                                                           "xmer", "hla_allele", paste0(allele, "affinity"), paste0(allele, "percentile_rank"), "c_term_aa", "processing_score"),
                                                              colclasses = c(unlist(lapply(peptides, class),
                                                                                    use.names = FALSE),
                                                                             "numeric", "character", "numeric", "numeric", "character", "numeric"))
    }
  }

  return(predictions)
}

performAffinityPredictions = function(peptides, allele, peptidelength) {
  if (length(peptides) < 1 ) {
    data = processPredictionOutput('', allele)
  } else {
    # generate random string to give to temp dir and temp file
    randomString = paste(sample(c(letters, LETTERS),
                                size = 8,
                                replace = T),
                         collapse = '')

    dir.create(file.path(runOptions$general$temporaryDirectoryPath, randomString))

    # write peptides to disk in temp dir
    invisible(sapply(seq(1, length(peptides), 1),
                     function(x)
                       write(x = sprintf("%s", peptides[x]),
                             file = file.path(runOptions$general$temporaryDirectoryPath, randomString, paste0(randomString, '_peps.fas')),
                             append = TRUE,
                             sep = "\n")))

    # perform predictions
    if (as.numeric(runParameters$panversion) >= 3) {
      command = paste0(# "nice -n 9 ",
        runOptions$predictors$netMHCpan,
        ' -a HLA-', gsub(pattern = '^([A-Z]{1}[0-9]{2})([0-9]{2})$',
                         replacement = '\\1:\\2',
                         x = allele),
        ' -l ', peptidelength,
        ' -p ',
        ' -f ', file.path(runOptions$general$temporaryDirectoryPath, randomString, paste0(randomString, '_peps.fas')),
        if (as.numeric(runParameters$panversion) == 4) paste(' -BA'))
    } else {
      command = paste0(# "nice -n 9 ",
        runOptions$predictors$netMHCpan,
        ' -a HLA-', gsub(pattern = '^([A-Z]{1}[0-9]{2})([0-9]{2})$',
                         replacement = '\\1:\\2',
                         x = allele),
        ' -l ', peptidelength,
        ' -p ',
        ' -f ', file.path(runOptions$general$temporaryDirectoryPath, randomString, paste0(randomString, '_peps.fas')),
        ' -tdir ', file.path(runOptions$general$temporaryDirectoryPath, randomString),
        ' -ic50')
    }

    output = system(command = command,
                    intern = TRUE)

    file.remove(file.path(runOptions$general$temporaryDirectoryPath, randomString, paste0(randomString, '_peps.fas')))
    file.remove(file.path(runOptions$general$temporaryDirectoryPath, randomString))

    data = processPredictionOutput(output, allele)
  }

  data = subset(x = data,
                select = colnames(data[, -dropNa(match(x = c("position", "variant_id", "peptide_score_log50k", "peptide_core", "Of", "Gp", "Gl", "Ip", "Il", "Icore"),
                                                       table = names(data))),
                                       with = FALSE]))
  return(data)
}


processPredictionOutput = function(raw_output, allele) {
  # regex the raw netMHC output
  raw_output = raw_output[-grep(pattern = "^\\#.+|^\\-.+|^Protein.+|pos.+|^HLA.+|^$",
                                ignore.case = TRUE,
                                x = raw_output)]
  raw_output = gsub(pattern = "^[[:blank:]]+| <= WB| <= SB",
                    replacement = "",
                    x = raw_output)
  clean_output = gsub(pattern = "[[:blank:]]+",
                      replacement = "\t",
                      x = raw_output)

  # parse cleaned output
  if (length(clean_output) > 0 & as.numeric(runParameters$panversion) >= 3) {
    data = as.data.table(read.table(text = clean_output,
                                    stringsAsFactors = FALSE))
    setnames(x = data, old = names(data), new = c("position", "hla_allele", "peptide", "peptide_core", "Of", "Gp", "Gl", "Ip", "Il", "Icore", "variant_id",
                                                  "peptide_score_log50k", paste0(allele, "affinity"), paste0(allele, "percentile_rank")))
  } else if (length(clean_output) > 0) {
    data = as.data.table(read.table(text = clean_output,
                                    stringsAsFactors = FALSE))
    setnames(x = data, old = names(data), new = c("position", "hla_allele", "peptide", "variant_id",
                                                  "peptide_score_log50k", paste0(allele, "affinity"), paste0(allele, "percentile_rank")))
  } else {
    data = emptyTableWithColumnNamesAndColumnClasses(colnames = c("position", "hla_allele", "peptide", "variant_id", "peptide_score_log50k", paste0(allele, "affinity"), paste0(allele, "percentile_rank")),
                                                     colclasses = c("numeric", "character", "character", "numeric", "numeric", "numeric", "numeric"))
  }
  return(data)
}


performProcessingPredictions = function(peptidestretch) {
  # generate random string to give to temp dir and temp file
  randomString = paste(sample(c(letters, LETTERS),
                              size = 8,
                              replace = T),
                       collapse = '')

  dir.create(file.path(runOptions$general$temporaryDirectoryPath, randomString))

  # write peptidestretch to disk in temp dir
  write(x = sprintf(">1\n%s", peptidestretch),
        file = file.path(runOptions$general$temporaryDirectoryPath, randomString, paste0(randomString, "_peptidestretch.fas")),
        append = TRUE,
        sep = "\n")

  # perform predictions
  output = system(command = paste(# "nice -n 9",
  	runOptions$predictors$netChop,
    '-tdir', file.path(runOptions$general$temporaryDirectoryPath, randomString),
  	file.path(runOptions$general$temporaryDirectoryPath, randomString, paste0(randomString, "_peptidestretch.fas")),
    sep = " "),
    intern = TRUE)

  file.remove(file.path(runOptions$general$temporaryDirectoryPath, randomString, paste0(randomString, "_peptidestretch.fas")))
  file.remove(file.path(runOptions$general$temporaryDirectoryPath, randomString))

  # perform regex on netChop output
  if (length(output) > 17) {
    output = output[-grep(pattern = "^\\#.+|^\\-.+|^Number.+|^NetChop.+|pos.+|^$",
                          x = output)]
    output = gsub(pattern = "^[[:blank:]]+",
                  replacement = "",
                  x = output)
    output = gsub(pattern = "[[:blank:]]+",
                  replacement = "\t",
                  x = output)

    # parse cleaned output
    data = as.data.table(read.table(text = output,
                                    stringsAsFactors = FALSE))
    data = data[, -match(x = c("V3", "V5"), table = names(data)), with = FALSE]
  } else {
    data = emptyTableWithColumnNamesAndColumnClasses(colnames = c("c_term_pos", "c_term_aa", "processing_score"),
                                                     colclasses = c("numeric", "character", "numeric"))
  }

  setnames(x = data,
           old = names(data),
           new = c("c_term_pos", "c_term_aa", "processing_score"))

  return(data)
}
