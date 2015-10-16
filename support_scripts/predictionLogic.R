## this support script allows live predictions through netMHCpan v2.4 for peptides that have proximal mutations (and there not likely to be in the FASdb) or are not in the FASdb for some other reason
## it will also allow us to, in the future, either (1) switch completely to 'live' predictions or (2) perform 'live' predictions for in-dels

performAffinityPredictions=function(fastafile,allele,peptidelength){
  ## use on local Mac
  # perform predictions on HPC from local Mac
#   sshconn=pipe(paste('ssh -l l.fanchi paranoid "',
#                       paste(netMHCpanpath,
#                             ' -a HLA-',gsub('^([A-Z]{1}[0-9]{2})([0-9]{2})$', '\\1:\\2', allele),
#                             ' -l ',xmer,
#                             ' -f ',peptides,'"',
#                             sep=""),
#                       sep="")) # start pipe to paranoid
#   output=readLines(sshconn) # start job on paranoid and read terminal output
#   close(sshconn) # close pipe
  
  ## use on HPC
  # perform predictions
  output=system(command=paste0("nice -n 9 ",
                               predictorPaths$netMHCpan,
                               ' -a HLA-',gsub('^([A-Z]{1}[0-9]{2})([0-9]{2})$', '\\1:\\2', allele),
                               ' -l ',peptidelength,
                               ' -f ',fastafile),
                intern=TRUE)
  
  # perform regex on netMHC output
  if (length(output)>57){
    output=output[-c(1:57)]
    output=output[-grep(pattern="^\\#.+|^\\-.+|^Protein.+|pos.+|^HLA.+|^$",
                        x=output)]
    output=gsub(pattern="^[[:blank:]]+| <= WB| <= SB",
                replacement="",
                x=output)
    output=gsub(pattern="[[:blank:]]+",
                replacement="\t",
                x=output)
    # read data into table
    data=as.data.table(read.table(text=output,stringsAsFactors = FALSE));data$V7=NULL
  } else{
    data=as.data.table(setNames(replicate(6,numeric(0), simplify = F), c("position","allele","peptide","mutation_id","pept_score",paste(allele,"affinity",sep=""))))
  }
  
  setnames(data,
            colnames(data),
            c("position","allele","peptide","mutation_id","pept_score",paste(allele,"affinity",sep="")))
  data=subset(x = data,
              select=colnames(data[,-match(c("position","allele","mutation_id","pept_score"),colnames(data)),with=FALSE]))
  return(data)
}

performProcessingPredictions=function(fastafile){
#   sshconn=pipe(paste('ssh -l l.fanchi paranoid "',
#                      paste(netChoppath,
#                            ' ',peptidestretch,'"',
#                            sep=""),
#                      sep="")) # start pipe to paranoid
#   output=readLines(sshconn) # start job on paranoid and read terminal output
#   close(sshconn)
  
  ## use on HPC
  # perform predictions
  output=system(command=paste("nice -n 9",
                              predictorPaths$netChop,
                              fastafile,
                              sep=" "),
                intern=TRUE)
  
  # perform regex on netChop output
  if (length(output)>17){
    output=output[-c(1:17)]
    output=output[-grep(pattern="^\\#.+|^\\-.+|^Number.+|^NetChop.+|pos.+|^$",
                        x=output)]
    output=gsub(pattern="^[[:blank:]]+",
                replacement="",
                x=output)
    output=gsub(pattern="[[:blank:]]+",
                replacement="\t",
                x=output)
    # read data into table
    data=as.data.table(read.table(text=output,stringsAsFactors = FALSE));data=data[,-match(c("V3","V5"),colnames(data)),with=FALSE]
  } else{
    data=as.data.table(setNames(replicate(3,numeric(0), simplify = F), c("c_term_pos","aa","processing_score")))
  }
  
  setnames(data,
           colnames(data),
           c("c_term_pos","aa","processing_score"))
  
  return(data)
}