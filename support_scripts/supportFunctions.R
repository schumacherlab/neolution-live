## Support functions

# define negation of %in% operator, returns logical vector for 'NOT IN'
'%ni%'=Negate('%in%')


# return directory where script is executed
thisDirectory=function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    # Rscript
    abspath=normalizePath(sub(needle, "", cmdArgs[match]))
    absdir=dirname(abspath)
    return(absdir)
  } else {
    # 'source'd via R console
    abspath=normalizePath(sys.frames()[[1]]$ofile)
    absdir=dirname(abspath)
    return(absdir)
  }
}

# allow natural sorting on multiple columns
multiMixedOrder=function(..., na.last = TRUE, decreasing = FALSE){
  do.call(order, c(
    lapply(list(...), function(l){
      if(is.character(l)){
        factor(l, levels=mixedsort(unique(l)))
      } else {
        l
      }
    }),
    list(na.last = na.last, decreasing = decreasing)
  ))
}

# check for presence of necessary tools
checkPredictorPaths=function(paths){
  if (!all(sapply(paths,file.exists))){
    stop(paste0("One or more predictors not found, please check runConfig.R. Expected paths are:\n",
                paste0(paths,
                       collapse = "\n")))
  }
}

# check if rna expression data is available for input tissue type
checkTissueTypeInput=function(inputtissue,availabletissues){
  if (!any(inputtissue==availabletissues)){
    stop(paste("Input tissue type",
                inputtissue,
                "not supported.\nAvailable tissues:",
                paste(availabletissues,
                      collapse=",")))
  }
}