## Support functions

# define negation of %in% operator, returns logical vector for 'NOT IN'
'%ni%'=Negate('%in%')

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
checkPredictorPaths=function(){
  predictors_list=c(netMHCpanpath,netChoppath)
  print("Checking for presence of predictors...")
  for(i in 1:length(predictors_list)){
    if (!file.exists(predictors_list[i])){
      stop(paste("Missing: ",predictors_list[i],", stopping script execution",sep=""))
    } else {
      print(paste("Present: ",predictors_list[i],", continuing...",sep=""))
    }
  }
}

# check if rna expression data is available for input tissue type
checkTissueTypeInput=function(inputtissue,availabletissues){
  if (any(inputtissue==availabletissues)){
    #message(paste("Input tissuetype: ",tissueTypeInput,sep=""))
  } else{
    #message(paste("Input tissue type (",tissueTypeInput,") has no available RNAseq data",sep=""))
    #message(paste("Available tissues: ", paste(availableTissueTypes,collapse=","),sep=""))
    stop("Input tissue type not supported")
  }
}