#####################################################################################################################################
### this script takes an input file (specified or as commandline argument) and queries the Schumi FAS db with the entries         ###
#####################################################################################################################################
### from v3 onwards, the script does not use ENSG and Strand data from input, but performs the lookup in the Schumi FASdb instead ###
#####################################################################################################################################

#####################################################################################################################################
### IMPORTANT NOTE BEFORE RUNNING: IF INPUT IS ALEXANDROV SNVS, MAKE SURE INCLUDE SAMPLE_ID IN CODE FOR PROXIMAL LOCATION FINDING ###
#####################################################################################################################################

# NOT IMPLEMENTED:  filtering on affinity/chop/pmbec (currently, database only contains rows that contain at least 1 affinity <500nM AKA 'slimmed-down' db)
# NOT IMPLEMENTED:  PBMEC score retrieval
# IMPLEMENT?:       Parallel queries? PREFERRED: run multiple tumor types at the same time by starting multiple instances

# changelog:
# v5.1 (20150820)
# - fixed a bug which prevented correct output when two transcripts were present on one chromosomal location
# - fixed a bug that caused peptides with proximal mutations to get incorrect or no processing scores
# - fixed a bug where merging two data.tables could result in cartesian joins (joining tables with duplicate 'by' values)
# - moved dbConnect calls into tryCatch blocks, so errors during connection establishment are caught as well
# v5 (20150723)
# - split main code logic into multiple R subscripts for clarity
# - refactored some code
# - changed script flow to do complete pipeline (lookup, filtering & self-sim) on one SNV before moving to next; output is written to disk after each SNV
# - added SNV info to output (ChromosomeID, ChromosomeLocation, Mutant Base)
# - (wt->mut) codon change now reported
# - (wt->mut) amino acid change now reported
# - location of mutated amino acids now reported
# - presence of proximal mutations now reported
# - epitopes not passing self-sim test are reported in tumor epitope list
# - changed perPepAffinityScores (high affinity values simulated) table to perPepAffinityScores_NEW (all affinity values calculated)
# - added filter in buildPeptides function, so now only peptides with AA mutations are returned (previously synonymous mutations could be returned)
# v4.3 (20141110)
# - RNAlikelihood scores now merge based on ENSG, not Gene Symbol
# v4.2 (20141106)
# - fixed bug where multiple ENSGs on one Chromosomal location weren't handled correctly
# - cases where database returns duplicates handled correctly (shouldn't be an issue if database is filled correctly)
# - added output of list of peptides not filtered for self-similarity
# v4.1 (20141104)
# - add wildtype epitopes to self-peptide list for self-similarity checking
# v4.0 (20141021)
# - extended self-similarity check performed on output peptides
# - RNAlikelihood scores included in output list
# v3.5 (20141016)
# - entries which have proximal mutations are now saved in separate list to allow quantification
# v3.4 (20140915)
# - optimized code and added some comments
# v3.3 (20140915)
# - added peptideStart to ProcessingScore lookup query (in case peptide is present multiple times in protein)
# v3.2 (20140915)
# - multiple mutations in same codon are now handled correctly
# v3.1 (20140912)
# - refactored code to fix scoping issues
# v3 (20140911)
# - implemented looking up ENSG and transcribed strand from new table in database (ENSGtoStrand)
# - removed dependency on input ENSG & Strand (both now obtained from FASdb)
# - script can now handle cases where one SNV returns two (or more) ENSGs/Transcribed strands
# v2.1 (20140910)
# - removed bug where, in case of one or more proximal SNVs, the mutant codons for the proximal mutations could be constructed incorrectly
# v2
# - proximal mutations are taken into account (mutations surrounding chromLoc/Pos)

suppressMessages(library(RMySQL))
suppressMessages(library(data.table))
suppressMessages(library(gtools))
suppressMessages(library(utils))
suppressMessages(library(parallel))
suppressMessages(library(foreach))
suppressMessages(library(doMC))

# get number of cores available for parallel processing
numberOfWorkers=5
registerDoMC(numberOfWorkers)
#message(paste("Using",numberOfWorkers,"cores for parallel processing",sep=" "))

# set min. affinity (in nM), min. chop score & min. RNA expression
affinityLimit=500
chopLimit=0.5
expressionLimit=0

# get commandline arguments and process them
argsFalse=commandArgs(FALSE)
argsTrue=commandArgs(TRUE)
scriptPath="/home/NKI/l.fanchi/working_environments/fasdb_run" # uncomment for production code
fileInput=argsTrue[1] # uncomment for production code
fileName=strsplit(fileInput,"/");fileName=gsub(".csv","",fileName[[1]][length(fileName[[1]])],fixed=TRUE)
tissueTypeInput=argsTrue[2]
lookupProgress=as.numeric(argsTrue[3])+1
xMer=9
setwd(scriptPath)

# load support functions
source(paste(scriptPath,"supportFunctions.R",sep="/"))

# fixed objects
#scriptPath="/home/NKI/l.fanchi/working_environments/fasdb_test"  # comment out for production code
#fileInput="/home/NKI/l.fanchi/working_environments/fas_database/input_data/20141016_rte_database_input.csv" # comment out for production code
currentDate=format(Sys.Date(), "%Y%m%d")
hlaTypes=c("A0101","A0201","A0301","B0702","B0801")
rnaExpressionFile="20150723_RNAexpression_likelihood_TCGA_IlluminaHiSeq.csv"

# prepare some empty lists
selfEpitopes=vector("list", length(hlaTypes))
affinityPredictions=vector("list",length(hlaTypes))
tumorEpitopes_HLA=vector("list",length(hlaTypes))
tumorEpitopes_filtered=vector("list",length(hlaTypes))

# create lists to hold wildtype output (only since dataWT and dataMT cannot be merged using slimmed-down FASdb)
normalEpitopes_HLA=vector("list",length(hlaTypes))
normalEpitopes_filtered=vector("list",length(hlaTypes))

# prepare list of tumor types with available RNAseq data
rnaExpressionData=fread(paste(scriptPath,"rna_expr_data",rnaExpressionFile,sep="/"),sep="auto",skip=0,header=TRUE)
availableTissueTypes=gsub("RNA_LKLIHD_",
                          "",
                          colnames(rnaExpressionData[,-c(which(colnames(rnaExpressionData) %in% c("ENSG_v58","ENSG","ENTREZ_ID","GENE_SYMBOL"))),with=FALSE]))

checkTissueTypeInput=function(){
  if (any(tissueTypeInput==availableTissueTypes)){
    #message(paste("Input tissuetype: ",tissueTypeInput,sep=""))
  } else{
    #message(paste("Input tissue type (",tissueTypeInput,") has no available RNAseq data",sep=""))
    #message(paste("Available tissues: ", paste(availableTissueTypes,collapse=","),sep=""))
    stop("Try again")
  }
}
checkTissueTypeInput()

source(paste(scriptPath,"queryLogic.R",sep="/"))
source(paste(scriptPath,"peptideConstructionLogic.R",sep="/"))
source(paste(scriptPath,"selfSimilarityLogic.R",sep="/"))
source(paste(scriptPath,"predictionLogic.R",sep="/"))
source(paste(scriptPath,"recordProgress.R",sep="/"))

# check availability of predictors
#checkNetAppPaths()

# empty tables to hold final query outputs
normalEpitopes=data.table()
tumorEpitopes=data.table()
entriesWithProximalMutations=data.table()

## import list of data for query building
input=fread(fileInput)
# remove exact duplicate entries
#input=unique(input,by=colnames(input)[-match("SAMPLE_ID",colnames(input))]) # we want this to only happen within a SAMPLE ID, remove this statement for production code
input=unique(input,by=colnames(input)) 

#message(paste("Input file:",fileName))
#message(paste("Number of input rows:",nrow(input)))

# create directory to hold output
dir.create(paste(scriptPath,currentDate,sep="/"),showWarnings=F)

# start looping over every row in input data and perform queries
#message("********** Query FASdb for SNV & peptide info..")
for (i in lookupProgress:nrow(input)){
  #cat('\rCurrently on row ',i,' of ',nrow(input),"; ",ceiling(i/nrow(input)*100),"%",sep="")
  
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
    
    if(!file.exists(paste(scriptPath,"/",currentDate,"/",currentDate,"_",fileName,"_prematureStops.csv",sep=""))){
      write(paste(colnames(output),collapse=",",sep=""),
                  paste(scriptPath,"/",currentDate,"/",currentDate,"_",fileName,"_prematureStops.csv",sep=""),
                  append = TRUE)
    }
    write(paste(output,collapse=",",sep=""),
                paste(scriptPath,"/",currentDate,"/",currentDate,"_",fileName,"_prematureStops.csv",sep=""),
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
    
    if(!file.exists(paste(scriptPath,"/",currentDate,"/",currentDate,"_",fileName,"_entriesWithProximalMutations.csv",sep=""))){
      write(paste(colnames(entriesWithProximalMutations),collapse=",",sep=""),
            paste(scriptPath,"/",currentDate,"/",currentDate,"_",fileName,"_entriesWithProximalMutations.csv",sep=""),
            append = TRUE)
    }

    for (j in nrow(entriesWithProximalMutations)){
      write(paste(entriesWithProximalMutations[j],collapse=",",sep=""),
            paste(scriptPath,"/",currentDate,"/",currentDate,"_",fileName,"_entriesWithProximalMutations.csv",sep=""),
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
        if (!file.exists(paste(scriptPath,"/",currentDate,"/",currentDate,"_",fileName,"_",tissueTypeInput,"_",hlaTypes[j],"_normalEpitopes_fasdb.csv", sep=""))){
          write(paste(colnames(normalEpitopes_filtered[[j]]),collapse=",",sep=""),
                paste(scriptPath,"/",currentDate,"/",currentDate,"_",fileName,"_",tissueTypeInput,"_",hlaTypes[j],"_normalEpitopes_fasdb.csv", sep=""),
                append = TRUE)
        }
        for (k in 1:nrow(normalEpitopes_filtered[[j]])){
          write(paste(normalEpitopes_filtered[[j]][k],collapse=",",sep=""),
                paste(scriptPath,"/",currentDate,"/",currentDate,"_",fileName,"_",tissueTypeInput,"_",hlaTypes[j],"_normalEpitopes_fasdb.csv", sep=""),
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
      #tumorEpitopes_HLA[[j]]=subset(tumorEpitopes_RNA, select=c("SAMPLE_ID","GENE_SYMBOL","GENE_SYMBOL_FASDB","ENSG","ENSG_v58","ENST","CODON_CHANGE","PEPTIDE_TUMOR","PEPTIDE_NORMAL",paste("AFFINITY_Kd.nM_TUMOR_",hlaTypes[i],sep=""),paste("AFFINITY_Kd.nM_NORMAL_",hlaTypes[i],sep=""),"CHOP_SCORE_TUMOR","CHOP_SCORE_NORMAL",paste("RNA_LKLIHD_",tissueTypeInput,sep=""))) # can only be used when using full FASdb
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
        #tumorEpitopes_filtered[[i]]=tumorEpitopes_filtered[[i]][mixedorder(tumorEpitopes_filtered[[i]]$AFFINITY_Kd.nM_TUMOR),];row.names(tumorEpitopes_filtered[[i]])=NULL # to sort lists on affinity; comment out for large lists
        #setcolorder(epitopeInput_filtered[[i]], c("SAMPLE_ID","GENE_SYMBOL","ENSG","GENE_SYMBOL_FASDB","ENSG_v58","ENST","CODON_CHANGE","PEPTIDE_TUMOR","PEPTIDE_NORMAL",paste("AFFINITY_Kd.nM_TUMOR_",hlaTypes[i],sep=""),paste("AFFINITY_Kd.nM_NORMAL_",hlaTypes[i],sep=""),"CHOP_SCORE_TUMOR","CHOP_SCORE_NORMAL",paste("RNA_LKLIHD_",tissueTypeInput,sep=""),"HAS_PROX_MUT","NON_SELF")) # only in case full FASdb is used
        setcolorder(tumorEpitopes_filtered[[j]],
                    c("SAMPLE_ID","CHR_ID","CHR_LOC","MUT_CODE","GENE_SYMBOL","ENSG","GENE_SYMBOL_FASDB","ENSG_v58","ENST","CODON_CHANGE","AA_CHANGE",
                      "AA_POSITION","PEPTIDE_TUMOR",paste("AFFINITY_Kd.nM_TUMOR_",hlaTypes[j],sep=""),"CHOP_SCORE_TUMOR",
                      paste("RNA_LKLIHD_",tissueTypeInput,sep=""),"HAS_PROX_MUT","NON_SELF"))
        if (!file.exists(paste(scriptPath,"/",currentDate,"/",currentDate,"_",fileName,"_",tissueTypeInput,"_",hlaTypes[j],"_tumorEpitopes_fasdb.csv", sep=""))){
          write(paste(colnames(tumorEpitopes_filtered[[j]]),collapse=",",sep=""),
                paste(scriptPath,"/",currentDate,"/",currentDate,"_",fileName,"_",tissueTypeInput,"_",hlaTypes[j],"_tumorEpitopes_fasdb.csv", sep=""),
                append = TRUE)
        }
        #write.csv(tumorEpitopes_filtered[[i]],paste(scriptPath,"/",currentDate,"/",currentDate,"_",fileName,"_",tissueTypeInput,"_",hlaTypes[i],"_fasdb.csv", sep=""),row.names=FALSE) # we want every epitope to be written to disk as soon as it is predicted
        for (k in 1:nrow(tumorEpitopes_filtered[[j]])){
          write(paste(tumorEpitopes_filtered[[j]][k],collapse=",",sep=""),
                paste(scriptPath,"/",currentDate,"/",currentDate,"_",fileName,"_",tissueTypeInput,"_",hlaTypes[j],"_tumorEpitopes_fasdb.csv", sep=""),
                append = TRUE)
        }
      }
    }
  }
  # send progress to "lookupProgress" table in FASdb
  dbUpdateRecord(dbtable="lookupProgress",data=data.table(dataset=fileName,progress=i,total=nrow(input)),primary="dataset",vars=c("progress","total"))
}
cat('\n')

#stop("Stopping now..")

detach(package:RMySQL)
detach(package:data.table)
detach(package:gtools)
detach(package:doMC)
detach(package:parallel)
detach(package:foreach)
