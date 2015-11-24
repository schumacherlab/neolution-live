## this support script allows live predictions through netMHCpan v2.4 for peptides that have proximal mutations (and there not likely to be in the FASdb) or are not in the FASdb for some other reason
## it will also allow us to, in the future, either (1) switch completely to 'live' predictions or (2) perform 'live' predictions for in-dels

performParallelPredictions = function(peptides, peptidestretch, allele, peptidelength) {
  if(nrow(peptides) > 0) {
    # perform affinity predictions
    affinityPredictions = performAffinityPredictions(peptides = peptides$peptide,
                                                     allele = allele,
                                                     peptidelength = peptidelength)
    
    # perform processing predictions
    processingPredictions = performProcessingPredictions(peptidestretch = peptidestretch)
    
    # merge prediction info
    predictions = merge(x = peptides,
                        y = affinityPredictions,
                        by = "peptide")
    
    predictions = merge(x = predictions,
                        y = processingPredictions,
                        by = "c_term_pos")
    
    predictions[, xmer := peptidelength]
  } else {
    predictions = as.data.table(setNames(replicate(n = length(c(names(peptides),
                                                                "xmer", "hla_allele", paste0(allele, "affinity"), "c_term_aa", "processing_score")),
                                                   expr = numeric(0),
                                                   simplify = FALSE),
                                         c(names(peptides),
                                           "xmer", "hla_allele", paste0(allele, "affinity"), "c_term_aa", "processing_score")))
  }
  
  return(predictions)
}

performAffinityPredictions = function(peptides, allele, peptidelength) {
  # generate random number to give to temp dir and temp file
  randomNumber = sample(x = 1:1000000000,
                        replace = FALSE,
                        size = 1)
  
  dir.create(paste0(temporaryDirectoryPath, "/", randomNumber))
  
  # write peptides to disk in temp dir
  invisible(sapply(seq(1, length(peptides), 1), function(x)
    write(x = sprintf(">%i\n%s", x, peptides[x]),
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
  output = system(command = paste0(# "nice -n 9 ",
    predictorPaths$netMHCpan,
    ' -a HLA-', gsub(pattern = '^([A-Z]{1}[0-9]{2})([0-9]{2})$',
                     replacement = '\\1:\\2',
                     x = allele),
    ' -l ', peptidelength,
    ' -f ', paste0(temporaryDirectoryPath, "/", randomNumber, "/", randomNumber, "_peps.fas"),
    ' -tdir ', paste0(temporaryDirectoryPath,"/",randomNumber)),
    intern = TRUE)
  
  file.remove(paste0(temporaryDirectoryPath, "/", randomNumber, "/", randomNumber, "_peps.fas"))
  file.remove(paste0(temporaryDirectoryPath, "/", randomNumber))
  
  # perform regex on netMHC output
  if (length(output) > 57) {
    output = output[-c(1:57)]
    output = output[-grep(pattern = "^\\#.+|^\\-.+|^Protein.+|pos.+|^HLA.+|^$",
                          x = output)]
    output = gsub(pattern = "^[[:blank:]]+| <= WB| <= SB",
                  replacement = "",
                  x = output)
    output = gsub(pattern = "[[:blank:]]+",
                  replacement = "\t",
                  x = output)
    # read data into table
    data = as.data.table(read.table(text = output,
                                    stringsAsFactors = FALSE));data$V7 = NULL
  } else {
    data = as.data.table(setNames(replicate(n = length(c("position", "hla_allele", "peptide", "variant_id", "pept_score", paste0(allele, "affinity"))),
                                            expr = numeric(0),
                                            simplify = FALSE),
                                  c("position", "hla_allele", "peptide", "variant_id", "pept_score", paste0(allele, "affinity"))))
  }
  
  setnames(x = data,
           old = names(data),
           new = c("position", "hla_allele", "peptide", "variant_id", "pept_score", paste0(allele, "affinity")))
  data = subset(x = data,
                select = colnames(data[, -match(x = c("position", "variant_id", "pept_score"), table = names(data)), with = FALSE]))
  return(data)
}

performProcessingPredictions = function(peptidestretch) {
  # generate random number to give to temp dir and temp file
  randomNumber = sample(x = 1:1000000000,
                        replace = FALSE,
                        size = 1)
  
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
    output = output[-c(1:17)]
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
    data = as.data.table(setNames(replicate(n = length(c("c_term_pos", "c_term_aa", "processing_score")),
                                            expr = numeric(0),
                                            simplify = FALSE),
                                  c("c_term_pos", "c_term_aa", "processing_score")))
  }
  
  setnames(x = data,
           old = names(data),
           new = c("c_term_pos", "c_term_aa", "processing_score"))
  
  return(data)
}