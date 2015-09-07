######################################################################################################################################
### IMPORTANT NOTE BEFORE RUNNING: IF INPUT IS ALEXANDROV LISTS, MAKE SURE INCLUDE SAMPLE_ID IN CODE FOR PROXIMAL LOCATION FINDING ###
######################################################################################################################################
# NOT IMPLEMENTED:  PBMEC score retrieval & PMBEC filtering

suppressMessages(library(RMySQL))
suppressMessages(library(data.table))
suppressMessages(library(gtools))
suppressMessages(library(utils))
suppressMessages(library(parallel))
suppressMessages(library(foreach))
suppressMessages(library(doMC))

#scriptPath="/home/NKI/l.fanchi/dev_environments/fasdb"  # use when debugging dev
scriptPath=thisDirectory()
setwd(scriptPath)

# load external scripts
source(paste(scriptPath,"supportFunctions.R",sep="/"))
source(paste(scriptPath,"runConfig.R",sep="/"))
source(paste(scriptPath,"queryLogic.R",sep="/"))
source(paste(scriptPath,"peptideConstructionLogic.R",sep="/"))
source(paste(scriptPath,"predictionLogic.R",sep="/"))
source(paste(scriptPath,"selfSimilarityLogic.R",sep="/"))
source(paste(scriptPath,"recordProgress.R",sep="/"))

# get commandline arguments and process them
argsFalse=commandArgs(FALSE)
argsTrue=commandArgs(TRUE)

fileInput=argsTrue[1] # uncomment for production code
fileName=strsplit(fileInput,"/");fileName=gsub(".csv","",fileName[[1]][length(fileName[[1]])],fixed=TRUE)
tissueTypeInput=argsTrue[2]
lookupProgress=as.numeric(argsTrue[3])+1

# register parallel back-end
registerDoMC(numberOfWorkers)

# prepare some empty lists
affinityPredictions=vector("list",length(hlaTypes))
normalEpitopes=data.table()
normalEpitopes_HLA=vector("list",length(hlaTypes))
normalEpitopes_filtered=vector("list",length(hlaTypes))
tumorEpitopes=data.table()
tumorEpitopes_HLA=vector("list",length(hlaTypes))
tumorEpitopes_filtered=vector("list",length(hlaTypes))
selfEpitopes=vector("list", length(hlaTypes))
entriesWithProximalMutations=data.table()

# read RNA expression data and make a list of which tissue types are available
rnaExpressionData=fread(paste(scriptPath,"rna_expr_data",rnaExpressionFile,sep="/"),sep="auto",skip=0,header=TRUE)
availableTissueTypes=gsub("RNA_LKLIHD_",
                          "",
                          colnames(rnaExpressionData[,-c(which(colnames(rnaExpressionData) %in% c("ENSG_v58","ENSG","ENTREZ_ID","GENE_SYMBOL"))),with=FALSE]))

checkTissueTypeInput(tissueTypeInput,availableTissueTypes)

# check availability of predictors
#checkPredictorPaths()

## import list of data for query building
input=fread(fileInput)
# remove exact duplicate entries
#input=unique(input,by=colnames(input)[-match("SAMPLE_ID",colnames(input))]) # we want this to only happen within a SAMPLE ID, remove this statement for production code
input=unique(input,by=colnames(input)) 

# create directory to hold logs/output, if necessary
dir.create(paste(scriptPath,"logs",sep="/"),showWarnings=FALSE)
dir.create(paste(scriptPath,"output",sep="/"),showWarnings=FALSE)

# start looping over every row in input data and perform queries
for (i in lookupProgress:nrow(input)){
  sampleID=as.character(input$SAMPLE_ID[i])
  chromID=input$ChromosomeID[i]
  chromLoc=input$ChromosomeLocation[i]
  mutantCode=input$MutantCode[i]
  
  ## start database queries and build output list
  # get mutation data for input mutation (need codonPosition to find proximal mutations and possibly multiple ENSGs on that genomic location)
  query1InputMutation=queryDatabaseWithChromIdAndLocation(chromID,
                                                          chromLoc)    
  
  if (nrow(query1InputMutation)>0){
    # get transcribed strand for input SNV
    query1InputMutation$Strand=rbindlist(lapply(query1InputMutation$ENSG,queryDatabaseWithENSGForStrand))$Strand
    setnames(query1InputMutation,
              colnames(query1InputMutation),
              c("ChromosomeID","ChromosomeLocation","ENSG","CodonPosition","ProteinPosition","Strand"))
  }else{
    output=data.table(SampleID=sampleID,
                      ChromosomeID=chromID,
                      ChromosomeLocation=chromLoc,
                      MutantCode=mutantCode,
                      StopReason="Genomic location not found, possibly transcript is missing")
    
    if(!file.exists(paste(scriptPath,"/output/",fileName,"_prematureStops.csv",sep=""))){
      write(paste(colnames(output),collapse=",",sep=""),
                  paste(scriptPath,"/output/",fileName,"_prematureStops.csv",sep=""),
                  append = TRUE)
    }
    write(paste(output,collapse=",",sep=""),
                paste(scriptPath,"/output/",fileName,"_prematureStops.csv",sep=""),
                append = TRUE)
    # send progress to "lookupProgress" table in FASdb
    dbUpdateRecord(dbtable="lookupProgress",data=data.table(dataset=fileName,progress=i,total=nrow(input)),primary="dataset",vars=c("progress","total"))
    next
  }

  epitopeOutput=getWildtypeAndMutatedEpitopes(sampleID,
                                              query1InputMutation,
                                              mutantCode,
                                              i)
  
  if (class(epitopeOutput)=="character"){
    output=data.table(SampleID=sampleID,
                      ChromosomeID=chromID,
                      ChromosomeLocation=chromLoc,
                      MutantCode=input$MutantCode[i],
                      StopReason=epitopeOutput)
    # send progress to "lookupProgress" table in FASdb
    dbUpdateRecord(dbtable="lookupProgress",data=data.table(dataset=fileName,progress=i,total=nrow(input)),primary="dataset",vars=c("progress","total"))
    next
  }else if(length(epitopeOutput)==3){
    normalEpitopes=epitopeOutput[[1]]
    tumorEpitopes=epitopeOutput[[2]]
    entriesWithProximalMutations=epitopeOutput[[3]]
  }else{
    # send progress to "lookupProgress" table in FASdb
    dbUpdateRecord(dbtable="lookupProgress",data=data.table(dataset=fileName,progress=i,total=nrow(input)),primary="dataset",vars=c("progress","total"))
    next
  }
  
  ## split code for processing normal, tumor and entries with prox muts
  # code for handling processing & logging of entries with proximal mutations
  if (nrow(entriesWithProximalMutations)>0){
    entriesWithProximalMutations=merge(entriesWithProximalMutations,
                                        subset(input,select=c("SAMPLE_ID","ChromosomeID","ChromosomeLocation")),by=c("ChromosomeID","ChromosomeLocation"))
    setcolorder(entriesWithProximalMutations,
                c("SAMPLE_ID","ChromosomeID","ChromosomeLocation","ENSG","CodonPosition","ProteinPosition","Strand"))
    
    if(!file.exists(paste(scriptPath,"/output/",fileName,"_entriesWithProximalMutations.csv",sep=""))){
      write(paste(colnames(entriesWithProximalMutations),collapse=",",sep=""),
            paste(scriptPath,"/output/",fileName,"_entriesWithProximalMutations.csv",sep=""),
            append = TRUE)
    }

    for (j in nrow(entriesWithProximalMutations)){
      write(paste(entriesWithProximalMutations[j],collapse=",",sep=""),
            paste(scriptPath,"/output/",fileName,"_entriesWithProximalMutations.csv",sep=""),
            append = TRUE) 
    }
  }
  
  # code for handling processing & logging of normal epitopes
  if (nrow(normalEpitopes)>0){
    # change column names to match input of self-sim script
    setnames(normalEpitopes,
              colnames(normalEpitopes),
              c("PEPTIDE_NORMAL",paste("AFFINITY_Kd.nM_NORMAL_",hlaTypes,sep=""),"ENSG_v58","PEPTIDE_START","CHOP_SCORE_NORMAL","SAMPLE_ID","CHR_ID","CHR_LOC",
                "MUT_CODE","GENE_SYMBOL_FASDB","ENST","CODON_CHANGE","AA_CHANGE","AA_POSITION","HAS_PROX_MUT"))
    setcolorder(normalEpitopes,
              c("SAMPLE_ID","CHR_ID","CHR_LOC","MUT_CODE","GENE_SYMBOL_FASDB","ENSG_v58","ENST","CODON_CHANGE","AA_CHANGE","AA_POSITION","PEPTIDE_START",
                "PEPTIDE_NORMAL",paste("AFFINITY_Kd.nM_NORMAL_",hlaTypes,sep=""),"CHOP_SCORE_NORMAL","HAS_PROX_MUT"))
    
    # merge RNAseq data for wildtype epitopes (only since dataWT and dataMT cannot be merged using slimmed-down FASdb)
    normalEpitopes_RNA=merge(normalEpitopes,
                              subset(rnaExpressionData, select=c("GENE_SYMBOL","ENSG_v58","ENSG",paste("RNA_LKLIHD_",tissueTypeInput,sep=""))),
                              by="ENSG_v58")
    normalEpitopes_RNA[,(paste("RNA_LKLIHD_",tissueTypeInput,sep="")):=lapply(.SD,as.numeric),.SDcols=paste("RNA_LKLIHD_",tissueTypeInput,sep="")]
    
    # split processing of different HLA types
    for(j in 1:length(hlaTypes)){
      # split predictions by HLA for wildtype epitopes (only since dataWT and dataMT cannot be merged using slimmed-down FASdb)
      normalEpitopes_HLA[[j]]=subset(normalEpitopes_RNA,
                                      select=c("SAMPLE_ID","CHR_ID","CHR_LOC","MUT_CODE","GENE_SYMBOL","GENE_SYMBOL_FASDB","ENSG","ENSG_v58","ENST","CODON_CHANGE",
                                               "AA_CHANGE","AA_POSITION","PEPTIDE_NORMAL",paste("AFFINITY_Kd.nM_NORMAL_",hlaTypes[j],sep=""),"CHOP_SCORE_NORMAL",
                                               paste("RNA_LKLIHD_",tissueTypeInput,sep=""),"HAS_PROX_MUT"))
    }
    
    # remove duplicates
    for (j in 1:length(hlaTypes)){
      # remove duplicates for wildtype epitopes (only since dataWT and dataMT cannot be merged using slimmed-down FASdb)
      normalEpitopes_HLA[[j]]=normalEpitopes_HLA[[j]][complete.cases(normalEpitopes_HLA[[j]]),]
      normalEpitopes_HLA[[j]]=unique(normalEpitopes_HLA[[j]],
                                      by=c('SAMPLE_ID','ENSG_v58','PEPTIDE_NORMAL','CHOP_SCORE_NORMAL'))
    }
    
    # filter predictions based on predicted binding affinity, chop score & RNA expression
    for (j in 1:length(hlaTypes)){
      # filtering for wildtype epitopes (only since dataWT and dataMT cannot be merged using slimmed-down FASdb)
      normalEpitopes_filtered[[j]]=normalEpitopes_HLA[[j]][normalEpitopes_HLA[[j]][[paste("AFFINITY_Kd.nM_NORMAL_",hlaTypes[j],sep="")]]<=affinityLimit & 
                                                            normalEpitopes_HLA[[j]][["CHOP_SCORE_NORMAL"]]>=chopLimit,] #& 
                                                            #normalEpitopes_HLA[[j]][[paste("RNA_LKLIHD_",tissueTypeInput,sep="")]]>expressionLimit,]
      row.names(normalEpitopes_filtered[[j]])=NULL
    }
    
    for (j in 1:length(hlaTypes)) {
      if (nrow(normalEpitopes_filtered[[j]])>0){
        setcolorder(normalEpitopes_filtered[[j]],
                    c("SAMPLE_ID","CHR_ID","CHR_LOC","MUT_CODE","GENE_SYMBOL","ENSG","GENE_SYMBOL_FASDB","ENSG_v58","ENST","CODON_CHANGE","AA_CHANGE",
                      "AA_POSITION","PEPTIDE_NORMAL",paste("AFFINITY_Kd.nM_NORMAL_",hlaTypes[j],sep=""),"CHOP_SCORE_NORMAL",
                      paste("RNA_LKLIHD_",tissueTypeInput,sep=""),"HAS_PROX_MUT"))
        if (!file.exists(paste(scriptPath,"/output/",fileName,"_",tissueTypeInput,"_",hlaTypes[j],"_normalEpitopes_fasdb.csv", sep=""))){
          write(paste(colnames(normalEpitopes_filtered[[j]]),collapse=",",sep=""),
                paste(scriptPath,"/output/",fileName,"_",tissueTypeInput,"_",hlaTypes[j],"_normalEpitopes_fasdb.csv", sep=""),
                append = TRUE)
        }
        for (k in 1:nrow(normalEpitopes_filtered[[j]])){
          write(paste(normalEpitopes_filtered[[j]][k],collapse=",",sep=""),
                paste(scriptPath,"/output/",fileName,"_",tissueTypeInput,"_",hlaTypes[j],"_normalEpitopes_fasdb.csv", sep=""),
                append = TRUE)
        }
      }
    }
  }
  
  if (nrow(tumorEpitopes)>0){
    # change column names to match input of self-sim script
    setnames(tumorEpitopes,
             colnames(tumorEpitopes),
             c("PEPTIDE_TUMOR",paste("AFFINITY_Kd.nM_TUMOR_",hlaTypes,sep=""),"ENSG_v58","PEPTIDE_START","CHOP_SCORE_TUMOR","SAMPLE_ID","CHR_ID","CHR_LOC",
               "MUT_CODE","GENE_SYMBOL_FASDB","ENST","CODON_CHANGE","AA_CHANGE","AA_POSITION","HAS_PROX_MUT"))
    setcolorder(tumorEpitopes,
                c("SAMPLE_ID","CHR_ID","CHR_LOC","MUT_CODE","GENE_SYMBOL_FASDB","ENSG_v58","ENST","CODON_CHANGE","AA_CHANGE","AA_POSITION","PEPTIDE_START",
                  "PEPTIDE_TUMOR",paste("AFFINITY_Kd.nM_TUMOR_",hlaTypes,sep=""),"CHOP_SCORE_TUMOR","HAS_PROX_MUT"))
    
    # merge epitopes with TCGA RNAseq info
    tumorEpitopes_RNA=merge(tumorEpitopes,
                            subset(rnaExpressionData,select=c("GENE_SYMBOL","ENSG_v58","ENSG",paste("RNA_LKLIHD_",tissueTypeInput,sep=""))),
                            by="ENSG_v58")
    tumorEpitopes_RNA[,(paste("RNA_LKLIHD_",tissueTypeInput,sep="")):=lapply(.SD,as.numeric),.SDcols=paste("RNA_LKLIHD_",tissueTypeInput,sep="")]
    
    # split processing of different HLA types
    for(j in 1:length(hlaTypes)){
      #tumorEpitopes_HLA[[j]]=subset(tumorEpitopes_RNA, select=c("SAMPLE_ID","GENE_SYMBOL","GENE_SYMBOL_FASDB","ENSG","ENSG_v58","ENST","CODON_CHANGE","PEPTIDE_TUMOR",
      #"PEPTIDE_NORMAL",paste("AFFINITY_Kd.nM_TUMOR_",hlaTypes[i],sep=""),paste("AFFINITY_Kd.nM_NORMAL_",hlaTypes[i],sep=""),"CHOP_SCORE_TUMOR","CHOP_SCORE_NORMAL",
      #paste("RNA_LKLIHD_",tissueTypeInput,sep=""))) # can only be used when using full FASdb
      tumorEpitopes_HLA[[j]]=subset(tumorEpitopes_RNA,
                                    select=c("SAMPLE_ID","CHR_ID","CHR_LOC","MUT_CODE","GENE_SYMBOL","GENE_SYMBOL_FASDB","ENSG","ENSG_v58","ENST","CODON_CHANGE",
                                             "AA_CHANGE","AA_POSITION","PEPTIDE_TUMOR",paste("AFFINITY_Kd.nM_TUMOR_",hlaTypes[j],sep=""),"CHOP_SCORE_TUMOR",
                                             paste("RNA_LKLIHD_",tissueTypeInput,sep=""),"HAS_PROX_MUT"))
    }
    
    # remove duplicates
    for (j in 1:length(hlaTypes)){
      # remove duplicates for tumor epitopes
      tumorEpitopes_HLA[[j]]=tumorEpitopes_HLA[[j]][complete.cases(tumorEpitopes_HLA[[j]]),]
      tumorEpitopes_HLA[[j]]=unique(tumorEpitopes_HLA[[j]],
                                    by=c('SAMPLE_ID',"ENSG_v58",'PEPTIDE_TUMOR','CHOP_SCORE_TUMOR'))
    }
    
    # filter predictions based on predicted binding affinity, chop score & RNA expression
    for (j in 1:length(hlaTypes)){
      # filtering for tumor epitopes
      tumorEpitopes_filtered[[j]]=tumorEpitopes_HLA[[j]][tumorEpitopes_HLA[[j]][[paste("AFFINITY_Kd.nM_TUMOR_",hlaTypes[j],sep="")]]<=affinityLimit & 
                                                           tumorEpitopes_HLA[[j]][["CHOP_SCORE_TUMOR"]]>=chopLimit,] #& 
                                                          #tumorEpitopes_HLA[[j]][[paste("RNA_LKLIHD_",tissueTypeInput,sep="")]]>expressionLimit,]
      row.names(tumorEpitopes_filtered[[j]])=NULL
    }
    
    # perform self-similarity check on tumor epitopes
    tumorEpitopes_filtered=foreach(j=1:length(hlaTypes)) %dopar% {
      if (nrow(tumorEpitopes_filtered[[j]])>0){
        tumorEpitopes_filtered[[j]]=performSelfSimilarityCheck(normalEpitopes_filtered[[j]],
                                                               tumorEpitopes_filtered[[j]],
                                                               selfEpitopes[[j]])
      }
      return(tumorEpitopes_filtered[[j]])
    }
    
    # write predicted epitopes to disk
    for (j in 1:length(hlaTypes)) {
      if (nrow(tumorEpitopes_filtered[[j]])>0){
        #setcolorder(epitopeInput_filtered[[i]], c("SAMPLE_ID","GENE_SYMBOL","ENSG","GENE_SYMBOL_FASDB","ENSG_v58","ENST","CODON_CHANGE","PEPTIDE_TUMOR","PEPTIDE_NORMAL",
        #paste("AFFINITY_Kd.nM_TUMOR_",hlaTypes[i],sep=""),paste("AFFINITY_Kd.nM_NORMAL_",hlaTypes[i],sep=""),"CHOP_SCORE_TUMOR","CHOP_SCORE_NORMAL",paste("RNA_LKLIHD_",tissueTypeInput,sep=""),
        #"HAS_PROX_MUT","NON_SELF")) # only in case full FASdb is used
        setcolorder(tumorEpitopes_filtered[[j]],
                    c("SAMPLE_ID","CHR_ID","CHR_LOC","MUT_CODE","GENE_SYMBOL","ENSG","GENE_SYMBOL_FASDB","ENSG_v58","ENST","CODON_CHANGE","AA_CHANGE",
                      "AA_POSITION","PEPTIDE_TUMOR",paste("AFFINITY_Kd.nM_TUMOR_",hlaTypes[j],sep=""),"CHOP_SCORE_TUMOR",
                      paste("RNA_LKLIHD_",tissueTypeInput,sep=""),"HAS_PROX_MUT","NON_SELF"))
        if (!file.exists(paste(scriptPath,"/output/",fileName,"_",tissueTypeInput,"_",hlaTypes[j],"_tumorEpitopes_fasdb.csv", sep=""))){
          write(paste(colnames(tumorEpitopes_filtered[[j]]),collapse=",",sep=""),
                paste(scriptPath,"/output/",fileName,"_",tissueTypeInput,"_",hlaTypes[j],"_tumorEpitopes_fasdb.csv", sep=""),
                append = TRUE)
        }
        for (k in 1:nrow(tumorEpitopes_filtered[[j]])){
          write(paste(tumorEpitopes_filtered[[j]][k],collapse=",",sep=""),
                paste(scriptPath,"/output/",fileName,"_",tissueTypeInput,"_",hlaTypes[j],"_tumorEpitopes_fasdb.csv", sep=""),
                append = TRUE)
        }
      }
    }
  }
  # send progress to "lookupProgress" table in FASdb
  dbUpdateRecord(dbtable="lookupProgress",data=data.table(dataset=fileName,progress=i,total=nrow(input)),primary="dataset",vars=c("progress","total"))
}
cat('\n')

detach(package:RMySQL)
detach(package:data.table)
detach(package:gtools)
detach(package:doMC)
detach(package:parallel)
detach(package:foreach)
