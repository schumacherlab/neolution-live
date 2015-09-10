## Peptide construction functions
# returns complementary base
strandOrientationSubstitutionTable=data.frame(orig=c("A","T","C","G"),
                                              subs=c("T","A","G","C"),
                                              stringsAsFactors=FALSE)

returnComplementaryBase=function(inputBase){
  return(strandOrientationSubstitutionTable[match(inputBase,strandOrientationSubstitutionTable$orig),]$subs)
}

# performs multiple base replacements in case >1 mutation in same codon occurs
returnCodonAfterMultiReplacement=function(inputMutations,duplicateProteinPosition){
  codonInfoForMultiReplacement=subset(inputMutations,subset=ProteinPosition==duplicateProteinPosition)
  newMutantCodon=codonInfoForMultiReplacement$WildtypeCodon[1]
    
  for(l in 1:nrow(codonInfoForMultiReplacement)){
    substr(x=newMutantCodon,
            start=codonInfoForMultiReplacement$CodonPosition[l],
            stop=codonInfoForMultiReplacement$CodonPosition[l])=codonInfoForMultiReplacement$MutantCode[l]
  }
  return(newMutantCodon)
}

# builds peptides from mutation context and mutant amino acid info
buildPeptides=function(mutationContext,mutationInfo){
  n=nrow(mutationContext)-8
  mutationContext=mutationContext[order(mutationContext$proteinPosition),]
  
  # make WT peptides
  WT=sapply(seq(1,n,1),function(i) paste(mutationContext$AA[i:(i+8)],collapse=""))
  
  # replace WT with MT aminoacid in mutationContext
  for(i in 1:nrow(mutationInfo)){
    mutationContext[match(mutationInfo$ProteinPosition[i],mutationContext$proteinPosition),]$AA=mutationInfo$MutantAminoAcid[i]
  }
  
  # make MT peptides
  MT=sapply(seq(1,n,1),function(i) paste(mutationContext$AA[i:(i+8)],collapse=""))
  
  # get peptideStart positions for relevant amino acids
  PeptideStart=mutationContext$proteinPosition[1:n]
  
  # make data table
  dt=data.table(peptideWT=WT,
                peptideMT=MT,
                PeptideStart=PeptideStart)
  
  # for each row, compare wt and mut and omit rows with identical peptides
  dt=dt[!sapply(seq(1,nrow(dt),1), function(x) identical(x=dt$peptideWT[x],y=dt$peptideMT[x])),]

  return(dt)
}

# returns the position of an amino acid, given a starting position of the peptide (relative to the protein) and the location of the amino acid in the protein
returnSNVAminoAcidPosition=function(peptideStart,proteinPosition){
  sapply(seq(1,length(proteinPosition),1), function(x) 
                                            {pos=proteinPosition[x]-peptideStart+1; if (pos<=9 & pos>0){return (pos)}})
}

returnChromosomalLocations=function(peptideStart, inputMutations){
  peptideEnd=peptideStart+8
  return(subset(x=inputMutations, subset=ProteinPosition %in% c(peptideStart:peptideEnd))$ChromosomeLocation)
}

returnMutantCode=function(peptideStart, inputMutations){
  peptideEnd=peptideStart+8
  return(subset(x=inputMutations, subset=ProteinPosition %in% c(peptideStart:peptideEnd))$MutantCode)
}

returnCodonPosition=function(peptideStart, inputMutations){
  peptideEnd=peptideStart+8
  return(subset(x=inputMutations, subset=ProteinPosition %in% c(peptideStart:peptideEnd))$CodonPosition)
}

returnCodonChange=function(peptideStart,proteinPosition,wildtypeCodon,mutantCodon){
  sapply(seq(1,length(proteinPosition),1), function(x) 
                                            {pos=proteinPosition[x]-peptideStart+1; if (pos<=9 & pos>0){return (paste(wildtypeCodon[x],mutantCodon[x],sep="->"))}})
}

returnAminoAcidChange=function(peptideStart,proteinPosition,wildtypeAA,mutantAA){
  sapply(seq(1,length(proteinPosition),1), function(x) 
                                            {pos=proteinPosition[x]-peptideStart+1; if (pos<=9 & pos>0){return (paste(wildtypeAA[x],mutantAA[x],sep="->"))}})
}

returnPeptidesWithProximalMutations=function(wtAndMtPeptides,proteinPositions){
  if (length(proteinPositions)>1){
    sapply(seq(1,nrow(wtAndMtPeptides),1), function(x)
      if (length(grep(x=seq(wtAndMtPeptides$PeptideStart[x],wtAndMtPeptides$PeptideStart[x]+8,1),pattern=paste(proteinPositions,collapse="|")))>1){
        return (wtAndMtPeptides[x,])
      }
      ,simplify = FALSE)
  } else if (length(proteinPositions)==1){
    sapply(seq(1,nrow(wtAndMtPeptides),1), function(x)
      return (wtAndMtPeptides[x,])
    ,simplify = FALSE)
  }
}


getWildtypeAndMutatedEpitopes=function(sampleID,query1InputMutation,mutantCode,index){
  epitopesWT=data.table()
  epitopesMT=data.table()
  entryWithProximalMutations=data.table()
  
  for(j in 1:nrow(query1InputMutation)){
    affinityPredictions=vector("list", 5)
    # find proximal mutations
    # make sure to include sample ID when analysing input files wih multiple patient identifiers!!!!! for debugging remove sample statement
    if (query1InputMutation$Strand[j]==1){
#       inputAndProximalMutations=subset(input,
#                                         subset=(ChromosomeID==ChromosomeID[index] & 
#                                                 (ChromosomeLocation>=(ChromosomeLocation[index]-(23+query1InputMutation$CodonPosition[j])) & 
#                                                   ChromosomeLocation<=(ChromosomeLocation[index]+(27-query1InputMutation$CodonPosition[j])))
#                                                 ))
      inputAndProximalMutations=subset(input,
                                       subset=(SAMPLE_ID==sampleID & 
                                                 (ChromosomeID==ChromosomeID[index] & 
                                                    (ChromosomeLocation>=(ChromosomeLocation[index]-(23+query1InputMutation$CodonPosition[j])) & 
                                                       ChromosomeLocation<=(ChromosomeLocation[index]+(27-query1InputMutation$CodonPosition[j]))))
                                               )
                                       )
    }else{
#       inputAndProximalMutations=subset(input,
#                                         subset=(ChromosomeID==ChromosomeID[index] & 
#                                                 (ChromosomeLocation>=(ChromosomeLocation[index]-(27-query1InputMutation$CodonPosition[j])) & 
#                                                   ChromosomeLocation<=(ChromosomeLocation[index]+(23+query1InputMutation$CodonPosition[j])))
#                                                 ))
      inputAndProximalMutations=subset(input,
                                       subset=(SAMPLE_ID==sampleID & 
                                                 (ChromosomeID==ChromosomeID[index] & 
                                                    (ChromosomeLocation>=(ChromosomeLocation[index]-(27-query1InputMutation$CodonPosition[j])) & 
                                                       ChromosomeLocation<=(ChromosomeLocation[index]+(23+query1InputMutation$CodonPosition[j]))))
                                               )
                                       )
    }
    
#     inputAndProximalMutations=unique(inputAndProximalMutations,
#                                       by=colnames(inputAndProximalMutations)[-1]) # add sample identifier when abovementioned is applicable (remove '[-1]')
     inputAndProximalMutations=unique(inputAndProximalMutations,
                                       by=colnames(inputAndProximalMutations)) # add sample identifier when abovementioned is applicable

    # get mutation data for input and proximal mutations
    if(nrow(inputAndProximalMutations)>1){
      # check orientation of Strand and, if necessary, replace with complementary base
      inputAndProximalMutations$MutantCode=sapply(seq(1,nrow(inputAndProximalMutations),1),function(i) 
                                                                                            if(query1InputMutation$Strand[j]==-1){
                                                                                              returnComplementaryBase(inputAndProximalMutations$MutantCode[i])
                                                                                            } else {
                                                                                              inputAndProximalMutations$MutantCode[i]
                                                                                            })
      
      # get mutation info (codonPositions and proteinPositions etc) for input and proximal mutations (doing input mutation query again, unnecessary, but easier, change in future?)
      query1InputAndProximalMutations=lapply(seq(1,nrow(inputAndProximalMutations),1),function(i) 
                                                                                        queryDatabaseWithChromIdAndLocation(inputAndProximalMutations$ChromosomeID[i],
                                                                                                                            inputAndProximalMutations$ChromosomeLocation[i]))
      query1InputAndProximalMutations=unique(do.call(rbind, query1InputAndProximalMutations))
      
      setnames(query1InputAndProximalMutations,
                colnames(query1InputAndProximalMutations),
                c("ChromosomeID","ChromosomeLocation","ENSG","CodonPosition","ProteinPosition"))
      inputAndProximalMutations=unique(subset(merge(inputAndProximalMutations,query1InputAndProximalMutations,by=c("ChromosomeID","ChromosomeLocation")),
                                        subset=ENSG==query1InputMutation$ENSG[j]))
      inputAndProximalMutations$HAS_PROX_MUT=TRUE
      
      # store entry that has proximal mutations
      entryWithProximalMutations=unique(rbindlist(list(entryWithProximalMutations,query1InputMutation)))
    }else{
      #setnames(query1InputMutation,colnames(query1InputMutation),c("ChromosomeID","ChromosomeLocation","ENSG","CodonPosition","ProteinPosition"))
      
      # no proximal mutations found, continue with single input mutation
      # check orientation of Strand and, if necessary, replace with complementary base
      if(query1InputMutation$Strand[j]==-1){
        inputAndProximalMutations$MutantCode=returnComplementaryBase(inputAndProximalMutations$MutantCode)
      }
      
      inputAndProximalMutations=unique(merge(inputAndProximalMutations,
                                              query1InputMutation[j,],
                                              by=c("ChromosomeID","ChromosomeLocation")))
      inputAndProximalMutations$HAS_PROX_MUT=FALSE
    }
    
    # find context of input mutation (aminoacids and corresponding codons -8 & +8 positions around mutated protein position)
    query2InputMutation=unique(queryDatabaseWithENSGAndProteinPositionFor9Mers(query1InputMutation$ENSG[j],
                                                                                inputAndProximalMutations[match(query1InputMutation$ChromosomeLocation[j],
                                                                                                                inputAndProximalMutations$ChromosomeLocation),]$ProteinPosition),
                        by="proteinPosition")
    if(nrow(query2InputMutation)>0){
      inputAndProximalMutations[,WildtypeCodon:=character(.N)];inputAndProximalMutations[,MutantCodon:=character(.N)]
      
      # find wildtype codon(s) for mutation(s)
      inputAndProximalMutations$WildtypeCodon=sapply(seq(1,nrow(inputAndProximalMutations),1),function(i) 
                                                                                                query2InputMutation[match(inputAndProximalMutations$ProteinPosition[i],
                                                                                                                          query2InputMutation$proteinPosition),]$codon)
      # substitute wildtype base for mutant base
      inputAndProximalMutations$MutantCodon=sapply(seq(1,nrow(inputAndProximalMutations),1),function(i) {
                                                                                              substr(inputAndProximalMutations$WildtypeCodon[i],
                                                                                                      inputAndProximalMutations$CodonPosition[i],
                                                                                                      inputAndProximalMutations$CodonPosition[i])=inputAndProximalMutations$MutantCode[i]
                                                                                              return(inputAndProximalMutations$WildtypeCodon[i])
                                                                                              },
                                                    USE.NAMES=FALSE)
      
      # if >1 mutation falls into the same codon (or protein position), perform multireplacement (replace multiple nucleotides in same codon)
      duplicatedProteinPositions=duplicated(inputAndProximalMutations$ProteinPosition)
      
      if(any(duplicatedProteinPositions)==TRUE){
        for(k in 1:length(which(duplicatedProteinPositions))){
          inputAndProximalMutations[ProteinPosition==inputAndProximalMutations[which(duplicatedProteinPositions)[k],]$ProteinPosition,]$MutantCodon=returnCodonAfterMultiReplacement(inputAndProximalMutations,
                                                                                                                                                                                      inputAndProximalMutations[which(duplicatedProteinPositions)[k],]$ProteinPosition)
        }
      }
    }else{
      return("Context of input mutation not found")
    }
    
    # find aminoacid(s) for mutant codon(s)
    query3InputAndProximalMutations=rbindlist(lapply(inputAndProximalMutations$MutantCodon,function(x) queryDatabaseWithCodon(x)))
    if(nrow(query3InputAndProximalMutations)>0){
      inputAndProximalMutations[,WildtypeAminoAcid:=character(.N)];inputAndProximalMutations[,MutantAminoAcid:=character(.N)]
      inputAndProximalMutations$WildtypeAminoAcid=sapply(seq(1,nrow(inputAndProximalMutations),1),function(i) 
                                                                                                    query2InputMutation[match(inputAndProximalMutations$ProteinPosition[i],query2InputMutation$proteinPosition),]$AA)
      inputAndProximalMutations$MutantAminoAcid=sapply(seq(1,nrow(inputAndProximalMutations),1),function(i) 
                                                                                                    query3InputAndProximalMutations[match(inputAndProximalMutations$MutantCodon[i],query3InputAndProximalMutations$codon),]$AA)
    }else{
      return("Mutant amino acid not found")
    }
    
    # find geneID & ENST for corresponding ENSG
    query4Output=queryDatabaseWithENSGForGeneIdAndENST(query1InputMutation$ENSG[j])
    if(nrow(query4Output)>0){
      inputAndProximalMutations$GeneID=as.character(query4Output$geneID[1])
      inputAndProximalMutations$ENST=as.character(query4Output$ENST[1])
    }else{
      inputAndProximalMutations$GeneID="N.A."
      inputAndProximalMutations$ENST="N.A."
    }
    
    # if 9 or more amino acids are present, build peptide(s)
    if(nrow(query2InputMutation)>=9){
      wtAndMtPeptides=buildPeptides(query2InputMutation,
                                    inputAndProximalMutations)
    }else{
      return("Less than 9 amino acids in context")
    }
    
    # if all WT peptides are equal to MT, skip
    if(all(wtAndMtPeptides$peptideWT==wtAndMtPeptides$peptideMT)){
      return("Wildtype and mutant peptides are equal")
    }
    
    outputAffinityWT=list()
    outputAffinityMT=list()
    outputProcessingWT=list()
    outputProcessingMT=list()
    
    for (k in 1:nrow(wtAndMtPeptides)){
      outputAffinityWT[[k]]=queryDatabaseWithPeptideForAffinityScore(wtAndMtPeptides$peptideWT[k])
      outputAffinityMT[[k]]=queryDatabaseWithPeptideForAffinityScore(wtAndMtPeptides$peptideMT[k])
      
      outputProcessingWT[[k]]=queryDatabaseWithENSGPeptideAndPeptideStartForProcessingScore(inputAndProximalMutations$ENSG[1],
                                                                                            wtAndMtPeptides$peptideWT[k],
                                                                                            wtAndMtPeptides$PeptideStart[k])
      outputProcessingMT[[k]]=queryDatabaseWithENSGPeptideAndPeptideStartForProcessingScore(inputAndProximalMutations$ENSG[1],
                                                                                            wtAndMtPeptides$peptideMT[k],
                                                                                            wtAndMtPeptides$PeptideStart[k])
    }
    
    outputAffinityWT=unique(rbindlist(outputAffinityWT));outputProcessingWT=unique(rbindlist(outputProcessingWT))
    outputAffinityMT=unique(rbindlist(outputAffinityMT));outputProcessingMT=unique(rbindlist(outputProcessingMT))
    
    # if peptide has proximal mutations, assume we're missing them and generate affinity data for all peptides for which we don't have data
    # note: some peptides with proximal mutations might have affinity data, since they could be in other transcripts
    # note 2: since we are working with the slimmed-down FASdb, we can't discriminate between: not-in-db due to prox mut and not-in-db due to >500nM, therefore we have predict affinity for all peptides with prox mut
    if(inputAndProximalMutations$HAS_PROX_MUT[1]==TRUE){
      peptidesWithProximalMutations=rbindlist(returnPeptidesWithProximalMutations(wtAndMtPeptides,unique(inputAndProximalMutations,by="ProteinPosition")$ProteinPosition))
      allPeptidesWithProximalMutations=data.table(peptide=c(peptidesWithProximalMutations$peptideWT,peptidesWithProximalMutations$peptideMT))
      peptidesWithProximalMutationsWithoutAffinity=unique(subset(x=allPeptidesWithProximalMutations,subset=peptide %ni% c(outputAffinityWT$peptide,outputAffinityMT$peptide)))
      
      randomNumbers=sample(1:1000000,3,replace=FALSE)
      
      if(nrow(peptidesWithProximalMutationsWithoutAffinity)>0){
        invisible(sapply(seq(1,nrow(peptidesWithProximalMutationsWithoutAffinity),1), function(x) 
                                                                              write(x=sprintf(">%i\n%s",x,peptidesWithProximalMutationsWithoutAffinity$peptide[x]),
                                                                                    file=paste0("./tmp/",randomNumbers[1],"_peps.fas"),
                                                                                    append=TRUE,
                                                                                    sep="\n")))
        
        # perform affinity predictions - change peps.fas path for production code
        affinityPredictions=foreach(k=1:length(hlaTypes)) %dopar% {
          affinityPredictions[[k]]=performAffinityPredictions(paste0("./tmp/",randomNumbers[1],"_peps.fas"),hlaTypes[k],peptideLength)
          return(affinityPredictions[[k]])
        }
        affinityPredictions=Reduce(function(x,y) merge(x,y,by="peptide"), affinityPredictions)
        
        affinityPredictionsWT=merge(data.table(peptide=peptidesWithProximalMutations$peptideWT),affinityPredictions,by = "peptide")
        affinityPredictionsMT=merge(data.table(peptide=peptidesWithProximalMutations$peptideMT),affinityPredictions,by = "peptide")
        
        outputAffinityWT=unique(rbindlist(list(outputAffinityWT,affinityPredictionsWT)),by="peptide")
        outputAffinityMT=unique(rbindlist(list(outputAffinityMT,affinityPredictionsMT)),by="peptide")
        
        file.remove(paste0("./tmp/",randomNumbers[1],"_peps.fas"))
      }
      
      # check if we need to generate processing data for some peptides
      if (all(allPeptidesWithProximalMutations$peptide %in% c(outputProcessingWT$peptide,outputProcessingMT$peptide))==FALSE){
        # if so, which peptides do we need processing scores for
        wtPeptidesWithoutProcessing=subset(x=peptidesWithProximalMutations,subset=peptideWT %ni% outputProcessingWT$peptide,select=c("peptideWT","PeptideStart"))
        mtPeptidesWithoutProcessing=subset(x=peptidesWithProximalMutations,subset=peptideMT %ni% outputProcessingMT$peptide,select=c("peptideMT","PeptideStart"))
        
        # create AA context of SNV (15 left and 15 right of SNV AA pos)
        wtAminoAcidContext=unique(queryDatabaseWithENSGAndProteinPositionForProcessing(query1InputMutation$ENSG[j],
                                                                                       inputAndProximalMutations[match(query1InputMutation$ChromosomeLocation[j],
                                                                                                                       inputAndProximalMutations$ChromosomeLocation),]$ProteinPosition),
                                  by="proteinPosition")
        wtAminoAcidContext=wtAminoAcidContext[order(wtAminoAcidContext$proteinPosition,decreasing = FALSE)];wtAminoAcidContext$codon=NULL
        mtAminoAcidContext=wtAminoAcidContext
        
        for(k in 1:nrow(inputAndProximalMutations)){
          mtAminoAcidContext[match(inputAndProximalMutations$ProteinPosition[k],mtAminoAcidContext$proteinPosition),]$AA=inputAndProximalMutations$MutantAminoAcid[k]
        }
        write(x=sprintf(">1\n%s",paste(wtAminoAcidContext$AA,collapse="")),
              file=paste0("./tmp/",randomNumbers[2],"_wtPeptidestretch.fas"),
              append=FALSE,
              sep="\n")
        
        write(x=sprintf(">1\n%s",paste(mtAminoAcidContext$AA,collapse="")),
              file=paste0("./tmp/",randomNumbers[3],"_mtPeptidestretch.fas"),
              append=FALSE,
              sep="\n")
        
        # perform processing predictions & change 'position' so positions reflect n-termini of 9mer peptides (should be changed in future if 10 & 11-mers are to be included!!)
        processingPredictionsWT=performProcessingPredictions(paste0("./tmp/",randomNumbers[2],"_wtPeptidestretch.fas"))
        processingPredictionsMT=performProcessingPredictions(paste0("./tmp/",randomNumbers[3],"_mtPeptidestretch.fas"))
        
        processingPredictionsWT=merge(data.table(peptide=wtPeptidesWithoutProcessing$peptideWT,peptideStart=wtPeptidesWithoutProcessing$PeptideStart),
                                      data.table(peptideStart=wtAminoAcidContext$proteinPosition-8,processingScore=processingPredictionsWT$processingScore),
                                      by="peptideStart")
        processingPredictionsMT=merge(data.table(peptide=mtPeptidesWithoutProcessing$peptideMT,peptideStart=mtPeptidesWithoutProcessing$PeptideStart),
                                      data.table(peptideStart=mtAminoAcidContext$proteinPosition-8,processingScore=processingPredictionsMT$processingScore),
                                      by="peptideStart")
        
        processingPredictionsWT[,"ENSG":=as.list(query1InputMutation$ENSG[j])]
        processingPredictionsMT[,"ENSG":=as.list(query1InputMutation$ENSG[j])]
        
        outputProcessingWT=unique(rbindlist(list(outputProcessingWT,processingPredictionsWT),use.names=TRUE),by=c("peptide","peptideStart"))
        outputProcessingMT=unique(rbindlist(list(outputProcessingMT,processingPredictionsMT),use.names=TRUE),by=c("peptide","peptideStart"))
        
        file.remove(paste0("./tmp/",randomNumbers[2],"_wtPeptidestretch.fas"))
        file.remove(paste0("./tmp/",randomNumbers[3],"_mtPeptidestretch.fas"))
      }
    }
    
    if (nrow(outputAffinityWT)>0 && any(outputAffinityWT$peptide %in% outputProcessingWT$peptide)){ # need to test for empty tables, as slimmed-down db can return NULL (good safety check anyway), also check if we have processing data for peptides
      outputWT=as.data.table(merge(outputAffinityWT,outputProcessingWT,by="peptide"));setnames(outputWT,"peptide","peptideWT")
      outputWT[,"SAMPLE_ID":=as.list(sampleID)]
      outputWT[,"CHR_ID":=as.list(query1InputMutation$ChromosomeID[j])]
      outputWT[,"CHR_LOC":=sapply(seq(1,nrow(outputWT),1), function(x) paste(returnChromosomalLocations(peptideStart=outputWT$peptideStart[x],
                                                                                                        inputMutations=inputAndProximalMutations),
                                                                              collapse="|"))]
      outputWT[,"MUT_CODE":=sapply(seq(1,nrow(outputWT),1), function(x) paste(returnMutantCode(peptideStart=outputWT$peptideStart[x],
                                                                                                inputMutations=inputAndProximalMutations),
                                                                              collapse="|"))]
      outputWT[,"POS_IN_CODON":=sapply(seq(1,nrow(outputWT),1), function(x) paste(returnCodonPosition(peptideStart=outputWT$peptideStart[x],
                                                                                                      inputMutations=inputAndProximalMutations),
                                                                                  collapse="|"))]
      outputWT[,"GENE_SYMBOL":=as.list(inputAndProximalMutations$GeneID[1])]
      outputWT[,"ENST":=as.list(inputAndProximalMutations$ENST[1])]
      outputWT[,"CODON_CHANGE":=sapply(seq(1,nrow(outputWT),1), function(x) 
                                                                  paste(unlist(sapply(outputWT$peptideStart[x],
                                                                                      returnCodonChange,
                                                                                      proteinPosition=unique(inputAndProximalMutations,by="ProteinPosition")$ProteinPosition,
                                                                                      wildtypeCodon=unique(inputAndProximalMutations,by="ProteinPosition")$WildtypeCodon,
                                                                                      mutantCodon=unique(inputAndProximalMutations,by="ProteinPosition")$MutantCodon)),
                                                                  collapse="|"))]
      outputWT[,"AA_CHANGE":=sapply(seq(1,nrow(outputWT),1), function(x) 
                                                                  paste(unlist(sapply(outputWT$peptideStart[x],
                                                                                      returnAminoAcidChange,
                                                                                      proteinPosition=unique(inputAndProximalMutations,by="ProteinPosition")$ProteinPosition,
                                                                                      wildtypeAA=unique(inputAndProximalMutations,by="ProteinPosition")$WildtypeAminoAcid,
                                                                                      mutantAA=unique(inputAndProximalMutations,by="ProteinPosition")$MutantAminoAcid)),
                                                                  collapse="|"))]
      
      outputWT[,"AA_POSITION":=sapply(seq(1,nrow(outputWT),1), function(x) 
                                                                  paste(unlist(sapply(outputWT$peptideStart[x],
                                                                                      returnSNVAminoAcidPosition,
                                                                                      proteinPosition=unique(inputAndProximalMutations,by="ProteinPosition")$ProteinPosition)),
                                                                  collapse="|"))]
      outputWT[,"HAS_PROX_MUT":=as.list(inputAndProximalMutations$HAS_PROX_MUT[1])]
    }else{
      outputWT=data.table()
    }
    
    if (nrow(outputAffinityMT)>0 && any(outputAffinityMT$peptide %in% outputProcessingMT$peptide)){ # need to test for empty tables, as slimmed-down db can return NULL (good safety check anyway), also check if we have processing data for peptides
      outputMT=as.data.table(merge(outputAffinityMT,outputProcessingMT,by="peptide"));setnames(outputMT,"peptide","peptideMT")
      outputMT[,"SAMPLE_ID":=as.list(sampleID)]
      outputMT[,"CHR_ID":=as.list(query1InputMutation$ChromosomeID[j])]
      outputMT[,"CHR_LOC":=sapply(seq(1,nrow(outputMT),1), function(x) paste(returnChromosomalLocations(peptideStart=outputMT$peptideStart[x],
                                                                                                        inputMutations=inputAndProximalMutations),
                                                                              collapse="|"))]
      outputMT[,"MUT_CODE":=sapply(seq(1,nrow(outputMT),1), function(x) paste(returnMutantCode(peptideStart=outputMT$peptideStart[x],
                                                                                                inputMutations=inputAndProximalMutations),
                                                                              collapse="|"))]
      outputMT[,"POS_IN_CODON":=sapply(seq(1,nrow(outputMT),1), function(x) paste(returnCodonPosition(peptideStart=outputMT$peptideStart[x],
                                                                                                      inputMutations=inputAndProximalMutations),
                                                                              collapse="|"))]
      outputMT[,"GENE_SYMBOL":=as.list(inputAndProximalMutations$GeneID[1])]
      outputMT[,"ENST":=as.list(inputAndProximalMutations$ENST[1])]
      outputMT[,"CODON_CHANGE":=sapply(seq(1,nrow(outputMT),1), function(x) 
                                                                    paste(unlist(sapply(outputMT$peptideStart[x],
                                                                                        returnCodonChange,proteinPosition=unique(inputAndProximalMutations,by="ProteinPosition")$ProteinPosition,
                                                                                        wildtypeCodon=unique(inputAndProximalMutations,by="ProteinPosition")$WildtypeCodon,
                                                                                        mutantCodon=unique(inputAndProximalMutations,by="ProteinPosition")$MutantCodon)),
                                                                    collapse="|"))]
      outputMT[,"AA_CHANGE":=sapply(seq(1,nrow(outputMT),1), function(x) 
                                                                    paste(unlist(sapply(outputMT$peptideStart[x],
                                                                                        returnAminoAcidChange,
                                                                                        proteinPosition=unique(inputAndProximalMutations,by="ProteinPosition")$ProteinPosition,
                                                                                        wildtypeAA=unique(inputAndProximalMutations,by="ProteinPosition")$WildtypeAminoAcid,
                                                                                        mutantAA=unique(inputAndProximalMutations,by="ProteinPosition")$MutantAminoAcid)),
                                                                    collapse="|"))]
      outputMT[,"AA_POSITION":=sapply(seq(1,nrow(outputMT),1), function(x) 
                                                                    paste(unlist(sapply(outputMT$peptideStart[x],
                                                                                        returnSNVAminoAcidPosition,
                                                                                        proteinPosition=unique(inputAndProximalMutations,by="ProteinPosition")$ProteinPosition)),
                                                                    collapse="|"))]
      outputMT[,"HAS_PROX_MUT":=as.list(inputAndProximalMutations$HAS_PROX_MUT[1])]
    }else{
      outputMT=data.table()
    }
    
    epitopesWT=rbind(epitopesWT,outputWT)
    epitopesMT=rbind(epitopesMT,outputMT)
    #output=merge(outputWT,outputMT,by="peptide") # can only merge WT and MT peptides when using full FASdb
  }
  return(list(epitopesWT,epitopesMT,entryWithProximalMutations))
}

