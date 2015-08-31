### script aggregates all specified logs and parses the files, finally outputs one table per tumor type
### since logs for the same tumor type can span multiple dates, this is fastest way of merging the logs

library(data.table)
library(dplyr)
library(gtools)
library(stringr)

## prepare empty list
mergedLogs=list()

## set directory where outputs can be found
scriptPath="/home/NKI/l.fanchi/working_environments/fasdb_run"
allDirectories=list.dirs(scriptPath)

## set regex pattern used to grab files from output directories
regexPattern="_tumorEpitopes"

## grab directories with output
outputDirectories=grep(pattern="2015[0-9]{4}$",
                       x=allDirectories,
                       value=TRUE)

## grab filenames based on regex pattern
logFilenames=list.files(path=outputDirectories,
                        pattern=regexPattern,
                        recursive=TRUE,
                        full.names=TRUE)

if(regexPattern=="_tumorEpitopes"){
  # list of HLA types
  hlaTypes=c("A0101","A0201","A0301","B0702","B0801")
  # make list of tumor type abbreviations
  rnaExpressionFile="20150723_RNAexpression_likelihood_TCGA_IlluminaHiSeq.csv"
  rnaExpressionData=fread(paste(scriptPath,"rna_expr_data",rnaExpressionFile,sep="/"),sep="auto",skip=0,header=TRUE)
  availableTissueTypes=gsub("RNA_LKLIHD_",
                            "",
                            colnames(rnaExpressionData[,-c(which(colnames(rnaExpressionData) %in% c("ENSG_v58","ENSG","ENTREZ_ID","GENE_SYMBOL"))),with=FALSE]))
  
  logFilenamesPerTumorType=lapply(seq(1,length(availableTissueTypes),1), function(x) subset(x=logFilenames,subset=grepl(pattern=availableTissueTypes[x],x=logFilenames)))
  
  logsPerTumorType=lapply(seq(1,length(logFilenamesPerTumorType),1), function(x) unique(rbindlist(lapply(X=logFilenamesPerTumorType[[x]],FUN=fread),use.names=TRUE,fill=TRUE)))
  
  invisible(sapply(logsPerTumorType,returnCountsPerHlaType,hlatypes=hlaTypes))
  
}

returnCountsPerHlaType=function(df,hlatypes){
  if(nrow(df)>0){
    tumorType=gsub(pattern="RNA_LKLIHD_",replacement="",x=grep(pattern="RNA_LKLIHD",x=names(df),value=TRUE))
    df_subsets=lapply(X=hlatypes, function(x)
    {
      tmp_subset=subset(df,
                        select=c("SAMPLE_ID","CHR_ID","CHR_LOC","MUT_CODE","GENE_SYMBOL","GENE_SYMBOL_FASDB","ENSG","ENSG_v58","ENST","CODON_CHANGE",
                                 "AA_CHANGE","AA_POSITION","PEPTIDE_TUMOR",paste("AFFINITY_Kd.nM_TUMOR_",x,sep=""),"CHOP_SCORE_TUMOR",
                                 paste("RNA_LKLIHD_",tumorType,sep=""),"HAS_PROX_MUT","NON_SELF"))
      tmp_complete=tmp_subset[complete.cases(tmp_subset),]
    })
    
    df_subsets_rna=lapply(X=hlatypes, function(x)
    {
      tmp_subset=subset(df,
                        select=c("SAMPLE_ID","CHR_ID","CHR_LOC","MUT_CODE","GENE_SYMBOL","GENE_SYMBOL_FASDB","ENSG","ENSG_v58","ENST","CODON_CHANGE",
                                 "AA_CHANGE","AA_POSITION","PEPTIDE_TUMOR",paste("AFFINITY_Kd.nM_TUMOR_",x,sep=""),"CHOP_SCORE_TUMOR",
                                 paste("RNA_LKLIHD_",tumorType,sep=""),"HAS_PROX_MUT","NON_SELF"))
      tmp_complete=tmp_subset[complete.cases(tmp_subset),]
      tmp_filtered=tmp_complete[tmp_complete[[paste("RNA_LKLIHD_",tumorType,sep="")]]>0,]
    })
    
    df_subsets_rna_selfsim=lapply(X=hlatypes, function(x)
    {
      tmp_subset=subset(df,
                        select=c("SAMPLE_ID","CHR_ID","CHR_LOC","MUT_CODE","GENE_SYMBOL","GENE_SYMBOL_FASDB","ENSG","ENSG_v58","ENST","CODON_CHANGE",
                                 "AA_CHANGE","AA_POSITION","PEPTIDE_TUMOR",paste("AFFINITY_Kd.nM_TUMOR_",x,sep=""),"CHOP_SCORE_TUMOR",
                                 paste("RNA_LKLIHD_",tumorType,sep=""),"HAS_PROX_MUT","NON_SELF"))
      tmp_complete=tmp_subset[complete.cases(tmp_subset),]
      tmp_filtered=tmp_complete[tmp_complete[[paste("RNA_LKLIHD_",tumorType,sep="")]]>0 &
                                  tmp_complete[["NON_SELF"]]==TRUE,]
    })
    
    cat("Unique",tumorType,"epitopes for\n",
        "HLA-A0101:",nrow(df_subsets[[1]]),"\t| non-self:",nrow(df_subsets_rna_selfsim[[1]]),"\n",
        "HLA-A0201:",nrow(df_subsets[[2]]),"\t| non-self:",nrow(df_subsets_rna_selfsim[[2]]),"\n",
        "HLA-A0301:",nrow(df_subsets[[3]]),"\t| non-self:",nrow(df_subsets_rna_selfsim[[3]]),"\n",
        "HLA-B0702:",nrow(df_subsets[[4]]),"\t| non-self:",nrow(df_subsets_rna_selfsim[[4]]),"\n",
        "HLA-B0801:",nrow(df_subsets[[5]]),"\t| non-self:",nrow(df_subsets_rna_selfsim[[5]]),"\n\n")
  }
}

## process premature stop logs
if(regexPattern=="_prematureStops"){
  ## load logs into list
  logs=lapply(logFilenames,fread)
  
  logDatasetNames=sapply(str_extract(string=logFilenames,pattern="_[A-Za-z_]+_"),
                         function(string) gsub(pattern="^_|_$",
                                               replacement="",
                                               x=string),
                         USE.NAMES=FALSE)
  
  parsedLogs=lapply(seq(1,length(logs),1), function(x) {
    dt=data.table(TumorType=rep(x=logDatasetNames[x],times=nrow(logs[[x]])),
                  SampleID=logs[[x]]$V1,
                  ChromosomeID=logs[[x]]$V2,
                  ChromosomeLocation=logs[[x]]$V3,
                  MutantCode=logs[[x]]$V4,
                  StopReason=paste0(logs[[x]]$V5,',',logs[[x]]$V6)
                  )
    dt[,ChromosomeID:=as.character(dt$ChromosomeID)]
    return(dt)
  }
  )
}

## bind logs of all tumor types
allParsedLogs=bind_rows(parsedLogs)

## take subsets per tumor type and write list to disk
for(i in 1:length(sort(unique(logDatasetNames)))){
  mergedLogs[[i]]=subset(allParsedLogs,subset=TumorType==sort(unique(logDatasetNames))[i])
  write.csv(x=mergedLogs[[i]],file=paste0("/home/NKI/l.fanchi/20150824_mergedPrematureStopsLogs/",sort(unique(logDatasetNames))[i],".csv"),row.names=FALSE)
}
