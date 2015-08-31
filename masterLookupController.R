## this script controls for which datasets the lookups are being performed
## if lookups need to be started, it also checks where lookups should start (in case of premature termination)

suppressMessages(library(RMySQL))
suppressMessages(library(gtools))
suppressMessages(library(data.table))
suppressMessages(library(parallel))
suppressMessages(library(foreach))
suppressMessages(library(doMC))

numberOfWorkers=5
registerDoMC(numberOfWorkers)

queryDatabaseForLookupProgress=function(){
  res=NULL
  attempt=1
  while(is.null(res) && attempt<=10){
    attempt=attempt+1
    dbConnection=dbConnect(MySQL(),
                           host="medoid",
                           user="l.fanchi",
                           password="MpRi1RKd",
                           dbname="SchumiDB")
    tryCatch(
{
  res<-dbGetQuery(dbConnection, paste('SELECT dataset, progress, total ',
                                      'FROM lookupProgress ',";",
                                      sep=""))
  dbDisconnect(dbConnection)
  return(as.data.table(res))
}
,
error=function(err){
  write(x=paste("Problem with query for lookupProgress. File:",fileName," | index:",i," | ",err),file=paste(scriptPath,"/",fileName,"_error.log",sep=""),append=TRUE)
  if (!is.null(res)){
    dbClearResult(res)
  }
  res=NULL
  dbDisconnect(dbConnection)
}
    )
  }
}

startDatabaseQuery=function(datasetpath,tissuetype,progressindex){
  system(paste("nice -n 9 Rscript /home/NKI/l.fanchi/working_environments/fasdb_run/queryDatabase_v5.1.R",datasetpath,tissuetype,progressindex))
}

scriptPath="/home/NKI/l.fanchi/working_environments/fasdb_run"
setwd(scriptPath)
datasetFolder="20150722_input_lists"

## make a table that contains which datasets use which tissuetype
tissueTypePerDataset=data.table(dataset=
                                  c("ALL","AML","Bladder","Breast","Cervix","CLL","Colorectum","Esophageal","Glioblastoma","Glioma_Low_Grade","Head_and_Neck",
                                    "Kidney_Chromophobe","Kidney_Clear_Cell","Kidney_Papillary","Liver","Lung_Adeno","Lung_Small_Cell","Lung_Squamous","Lymphoma_B-cell",
                                    "Medulloblastoma","Melanoma","Myeloma","Neuroblastoma","Ovary","Pancreas","Pilocytic_Astrocytoma","Prostate","Stomach","Thyroid","Uterus"),
                                tissueType=
                                  c("DLBC","LAML","BLCA","BRCA","DLBC","CESC","COAD","STAD","GBM","LGG","HNSC","KICH","KIRC","KIRP","LIHC","LUAD","LUSMC","LUSC","DLBC","GBM","SKCM","DLBC","ACC","OV","PAAD","LGG","PRAD","STAD","THCA","UCEC")
)
datasetOrder=as.integer(c(5,19,12,22,17,15,9,14,20,26,25,6,3,1,8,2,18,23,28,10,30,29,13,27,11,21,24,7,16,4))

# get all dataset files
allDatasets=list.files(path=paste(scriptPath,datasetFolder,sep="/"),full.names=FALSE)

# read progress from "lookupProgress" table in FASdb for all datasets
progress=queryDatabaseForLookupProgress()
allProgressPerDataset=merge(x = tissueTypePerDataset,y = progress,by="dataset")
allProgressPerDataset=allProgressPerDataset[c(datasetOrder),]

# exclude datasets which are complete
inProgressPerDataset=rbindlist(sapply(seq(1,nrow(allProgressPerDataset),1), function(x)
  if (!(allProgressPerDataset$progress[x]==allProgressPerDataset$total[x])){
    return (allProgressPerDataset[x,])
  },simplify = FALSE))

# create commandline calls for necessary datasets (check which are running and/or which have been completed)
#cmdlineCalls=sapply(seq(1,nrow(inProgressPerDataset),1), function(x) paste("nice -n 9 Rscript /home/NKI/l.fanchi/working_environments/fasdb_run/queryDatabase_v5.R ",datasetFolder,"/",inProgressPerDataset$dataset[x],".csv ",inProgressPerDataset$tissueType[x]," ",inProgressPerDataset$progress[x]," &",sep=""))

write(paste(Sys.time(),
            " - Job list:\n",
            paste(inProgressPerDataset$dataset,inProgressPerDataset$tissueType,inProgressPerDataset$progress,collapse="\n"),"\n",sep=""),
            file="/home/NKI/l.fanchi/masterLookupController.log",append=TRUE)

# start lookups with index and tissuetype
# invisible(x=foreach(i=1:length(cmdlineCalls)) %dopar% {
#   pid=Sys.getpid()
#   system(cmdlineCalls[i])
#   while(pid %in% list(ofPIDs)){
#     Sys.sleep(3600)
#   }
# })

invisible(x=foreach(i=1:nrow(inProgressPerDataset)) %dopar% {
  startDatabaseQuery(paste(scriptPath,"/",datasetFolder,"/",inProgressPerDataset$dataset[i],".csv",sep=""),inProgressPerDataset$tissueType[i],inProgressPerDataset$progress[i])
})