## this script controls for which datasets the lookups are to be performed
## if lookups need to be started, it also checks where lookups should start (in case of premature lookup termination)

suppressMessages(library(RMySQL))
suppressMessages(library(gtools))
suppressMessages(library(data.table))
suppressMessages(library(parallel))
suppressMessages(library(foreach))
suppressMessages(library(doMC))

numberOfWorkers=5
registerDoMC(numberOfWorkers)

# lookup the progress for all datasets
queryDatabaseForLookupProgress=function(){
	res=NULL
	attempt=1
	while(is.null(res) && attempt<=10){
		attempt=attempt+1

		tryCatch(
		{
			dbConnection=dbConnect(MySQL(),
				host="medoid",
				user="l.fanchi",
				password="MpRi1RKd",
				dbname="SchumiDB")
			res<-dbGetQuery(dbConnection, paste('SELECT dataset, progress, total ',
												'FROM lookupProgress ',";",
				sep=""))
			dbDisconnect(dbConnection)
			return(as.data.table(res))
		}
		,
		error=function(err){
			write(x=paste0("Problem with query for lookupProgress | ",err),
				file=paste(scriptPath,"/logs/","getProgress_errors.log",sep=""),
				append=TRUE)

			if (!is.null(res)){
				dbClearResult(res)
			}

			res=NULL

			if(exists("dbConnection")){
				dbDisconnect(dbConnection)
			}
		})
	}
}

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

startDatabaseQuery=function(datasetpath,tissuetype,progressindex){
  system(paste0("nice -n 9 Rscript ",scriptPath,"/queryDatabase.R ",datasetpath," ",tissuetype," ",progressindex))
}

#scriptPath="/home/NKI/l.fanchi/working_environments/fasdb"
scriptPath=thisDirectory()
datasetFolder="20150722_input_lists"
setwd(scriptPath)

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
		},
		simplify = FALSE))

dir.create(paste(scriptPath,"logs",sep="/"),showWarnings=FALSE)
write(paste(Sys.time(),
            " - Job list:\n",
            paste(inProgressPerDataset$dataset,inProgressPerDataset$tissueType,inProgressPerDataset$progress,collapse="\n"),"\n",sep=""),
            file=paste(scriptPath,"/logs/","masterLookupController.log",sep=""),
            append=TRUE)

invisible(x=foreach(i=1:nrow(inProgressPerDataset)) %dopar% {
  startDatabaseQuery(datasetpath = paste(scriptPath,"/",datasetFolder,"/",inProgressPerDataset$dataset[i],".csv",sep=""),
  					tissuetype = inProgressPerDataset$tissueType[i],
  					progressindex = inProgressPerDataset$progress[i])
})