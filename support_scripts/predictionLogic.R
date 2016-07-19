## this support script contains functions for peptide affinity and processing score predictions

performParallelPredictions = function(peptides, peptidestretch, allele, peptidelength) {
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

  return(predictions)
}

performAffinityPredictions = function(peptides, allele, peptidelength) {
  # generate random number to give to temp dir and temp file
  randomNumber = ceiling(runif(n = 1,
                               min = 0,
                               max = 10 ^ 10))

  dir.create(paste0(temporaryDirectoryPath, "/", randomNumber))

  # write peptides to disk in temp dir
  invisible(sapply(seq(1, length(peptides), 1), function(x)
    write(x = sprintf("%s", peptides[x]),
          file = paste0(temporaryDirectoryPath, "/", randomNumber, "/", randomNumber, "_peps.fas"),
          append = TRUE,
          sep = "\n")))

  ## use on local Mac
  # perform predictions on HPC from local Mac
  # sshconn = pipe(paste0('ssh -l l.fanchi paranoid "',
  #                    paste0("nice -n 9 ",
  #                           predictorPaths$netMHCpan,
  #                           ' -a HLA-',gsub('^([A-Z]{1}[0-9]{2})([0-9]{2})$', '\\1:\\2', allele),
  #                           ' -l ',peptidelength,
  #                           ' -f ',fastafile,'"'))) # start pipe to paranoid
  # output = readLines(sshconn) # start job on paranoid and read terminal output
  # close(sshconn) # close pipe

  ## use on HPC
  # perform predictions
  command = ifelse(test = runParameters$panversion == "3.0",
                   yes = paste0(# "nice -n 9 ",
                     predictorPaths$netMHCpan,
                     ' -a HLA-', gsub(pattern = '^([A-Z]{1}[0-9]{2})([0-9]{2})$',
                                      replacement = '\\1:\\2',
                                      x = allele),
                     ' -l ', peptidelength,
                     ' -p ',
                     ' -f ', paste0(temporaryDirectoryPath, "/", randomNumber, "/", randomNumber, "_peps.fas")),
                   no = paste0(# "nice -n 9 ",
                     predictorPaths$netMHCpan,
                     ' -a HLA-', gsub(pattern = '^([A-Z]{1}[0-9]{2})([0-9]{2})$',
                                      replacement = '\\1:\\2',
                                      x = allele),
                     ' -l ', peptidelength,
                     ' -p ',
                     ' -f ', paste0(temporaryDirectoryPath, "/", randomNumber, "/", randomNumber, "_peps.fas"),
                     ' -tdir ', paste0(temporaryDirectoryPath,"/",randomNumber),
                     ' -ic50')
  )
  output = system(command = command,
                  intern = TRUE)

  file.remove(paste0(temporaryDirectoryPath, "/", randomNumber, "/", randomNumber, "_peps.fas"))
  file.remove(paste0(temporaryDirectoryPath, "/", randomNumber))

  # perform regex on netMHC output
  processPreditionOutput = function(raw_output) {
    raw_output = raw_output[-grep(pattern = "^\\#.+|^\\-.+|^Protein.+|pos.+|^HLA.+|^$",
                                  ignore.case = TRUE,
                                  x = raw_output)]
    raw_output = gsub(pattern = "^[[:blank:]]+| <= WB| <= SB",
                      replacement = "",
                      x = raw_output)
    clean_output = gsub(pattern = "[[:blank:]]+",
                        replacement = "\t",
                        x = raw_output)
    if (length(clean_output) > 0 & runParameters$panversion == "3.0") {
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

  data = processPreditionOutput(output)

  data = subset(x = data,
                select = colnames(data[, -dropNa(match(x = c("position", "variant_id", "peptide_score_log50k", "peptide_core", "Of", "Gp", "Gl", "Ip", "Il", "Icore"),
                                                        table = names(data))),
                                       with = FALSE]))
  return(data)
}

performProcessingPredictions = function(peptidestretch) {
  # generate random number to give to temp dir and temp file
  randomNumber = ceiling(runif(n = 1,
                               min = 0,
                               max = 10 ^ 10))

  dir.create(paste0(temporaryDirectoryPath, "/", randomNumber))

  # write peptidestretch to disk in temp dir
  write(x = sprintf(">1\n%s", peptidestretch),
        file = paste0(temporaryDirectoryPath, "/", randomNumber, "/", randomNumber, "_peptidestretch.fas"),
        append = TRUE,
        sep = "\n")

  #   sshconn = pipe(paste('ssh -l l.fanchi paranoid "',
  #                      paste(netChoppath,
  #                            ' ',peptidestretch,'"',
  #                            sep = ""),
  #                      sep = "")) # start pipe to paranoid
  #   output = readLines(sshconn) # start job on paranoid and read terminal output
  #   close(sshconn)

  ## use on HPC
  # perform predictions
  output = system(command = paste(# "nice -n 9",
    predictorPaths$netChop,
    '-tdir', paste0(temporaryDirectoryPath, "/", randomNumber),
    paste0(temporaryDirectoryPath, "/", randomNumber, "/", randomNumber, "_peptidestretch.fas"),
    sep = " "),
    intern = TRUE)

  file.remove(paste0(temporaryDirectoryPath, "/", randomNumber, "/", randomNumber, "_peptidestretch.fas"))
  file.remove(paste0(temporaryDirectoryPath, "/", randomNumber))

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
    # read data into table
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
