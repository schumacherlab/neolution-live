## functions used to report statistics on generated epitopes

## get some stats
for(i in 1:length(hlaTypes)){
  if (nrow(tumorEpitopes_filtered[[i]])>0){
    numberOfForeignEpitopes[i]=as.character(nrow(tumorEpitopes_filtered[[i]]))
  }else{
    numberOfForeignEpitopes[i]="0"
  }
}

## print summary of all predictions
message("----------------------------------------------------------------------------------------------------")
for(k in 1:length(hlaTypes)){
  message(paste(numberOfForeignEpitopes[k]," epitopes (",paste(xMer,collapse="-"),"mers) predicted for dataset: ",fileName," | Allele: ",hlaTypes[k],sep=""))
  message(paste("Total epitopes: ",numberOfTotalEpitopes[k]," | After affinity & chop filter: ",numberOfFilteredEpitopes[k]," | After self-sim check: ",numberOfForeignEpitopes[k],sep=""))
  message(paste("Affinity cutoff: ",affinityLimit,"nM | Chop score cutoff: ",chopLimit," | RNA expression cutoff: ",expressionLimit,sep=""))
  message("----------------------------------------------------------------------------------------------------")
}