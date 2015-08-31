## Query functions
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
        write(x=paste("Problem with query for lookupProgress. File:",fileName," | index:",i," | ",err),file=paste(scriptPath,"/logs/",fileName,"_error.log",sep=""),append=TRUE)
        if (!is.null(res)){
          dbClearResult(res)
        }
        res=NULL
        dbDisconnect(dbConnection)
      }
    )
  }
}

queryDatabaseWithChromIdAndLocation=function(chromID,chromLoc){
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
        res<-dbGetQuery(dbConnection, paste('SELECT chromosomeID, chromosomeLocation, ENSG, codonPosition, proteinPosition ',
                                            'FROM perChromosomeLocationENSG ',
                                            'WHERE chromosomeID="',chromID,'" ',
                                            'AND chromosomeLocation=',chromLoc,";",
                                            sep=""))
        dbDisconnect(dbConnection)
        return(as.data.table(res))
      }
      ,
      error=function(err){
        write(x=paste("Problem with query for SNV info. File:",fileName," | index:",i," | ",err),file=paste(scriptPath,"/logs/",fileName,"_error.log",sep=""),append=TRUE)
        if (!is.null(res)){
          dbClearResult(res)
        }
        res=NULL
        dbDisconnect(dbConnection)
      }
    )
  }
}

queryDatabaseWithENSGForStrand=function(ENSG){
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
        res<-dbGetQuery(dbConnection, paste('SELECT ENSG, Strand ',
                                            'FROM ENSGtoStrand ',
                                            'WHERE ENSG="',ENSG,'"',";",
                                            sep=""))
        dbDisconnect(dbConnection)
        return(as.data.table(res))
      }
      ,
      error=function(err){
        write(x=paste("Problem with query for strand. File:",fileName," | index:",i," | ",err),file=paste(scriptPath,"/logs/",fileName,"_error.log",sep=""),append=TRUE)
        if (!is.null(res)){
          dbClearResult(res)  
        }
        res=NULL
        dbDisconnect(dbConnection)
      }
    )
  }
}

queryDatabaseWithENSGAndProteinPositionFor9Mers=function(ENSG,protPos){
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
        res<-dbGetQuery(dbConnection, paste('SELECT ENSG, proteinPosition, codon, AA ',
                                            'FROM perProteinAAC ',
                                            'WHERE ENSG="',ENSG,'" ',
                                            'AND proteinPosition BETWEEN ',as.numeric(protPos)-8,' AND ',as.numeric(protPos)+8,";",
                                            sep=""))
        dbDisconnect(dbConnection)
        return(as.data.table(res))
      }
      ,
      error=function(err){
        write(x=paste("Problem with query for 9mers. File:",fileName," | index:",i," | ",err),file=paste(scriptPath,"/logs/",fileName,"_error.log",sep=""),append=TRUE)
        if (!is.null(res)){
          dbClearResult(res)  
        }
        res=NULL
        dbDisconnect(dbConnection)
      }
    )  
  }
}

queryDatabaseWithENSGAndProteinPositionForProcessing=function(ENSG,protPos){
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
        res<-dbGetQuery(dbConnection, paste('SELECT ENSG, proteinPosition, codon, AA ',
                                            'FROM perProteinAAC ',
                                            'WHERE ENSG="',ENSG,'" ',
                                            'AND proteinPosition BETWEEN ',as.numeric(protPos)-19,' AND ',as.numeric(protPos)+19,";",
                                            sep=""))
        dbDisconnect(dbConnection)
        return(as.data.table(res))
      }
      ,
      error=function(err){
        write(x=paste("Problem with query for peptideStretch. File:",fileName," | index:",i," | ",err),file=paste(scriptPath,"/logs/",fileName,"_error.log",sep=""),append=TRUE)
        if (!is.null(res)){
          dbClearResult(res)  
        }
        res=NULL
        dbDisconnect(dbConnection)
      }
    )  
  }
}

queryDatabaseWithCodon=function(codon){
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
        res<-dbGetQuery(dbConnection, paste('SELECT codon, AA ',
                                            'FROM codonToAAC ',
                                            'WHERE codon="',codon,'"',";",
                                            sep=""))
        dbDisconnect(dbConnection)
        return(as.data.table(res))
      }
      ,
      error=function(err){
        write(x=paste("Problem with query for aminoAcid. File:",fileName," | index:",i," | ",err),file=paste(scriptPath,"/logs/",fileName,"_error.log",sep=""),append=TRUE)
        if (!is.null(res)){
          dbClearResult(res)  
        }
        res=NULL
        dbDisconnect(dbConnection)
      }
    )  
  }
}

queryDatabaseWithENSGPeptideAndPeptideStartForProcessingScore=function(ENSG,peptide,position){
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
        res<-dbGetQuery(dbConnection, paste('SELECT ENSG, peptide, peptideStart, processingScore ',
                                            'FROM perGenePepProcessingScores ',
                                            'WHERE ENSG="',ENSG,'" ',
                                            'AND peptide="',peptide,'" ',
                                            'AND peptideStart="',position,'"',";",
                                            sep=""))
        dbDisconnect(dbConnection)
        return(as.data.table(res))
      }
      ,
      error=function(err){
        write(x=paste("Problem with query for processingScore. File:",fileName," | index:",i," | ",err),file=paste(scriptPath,"/logs/",fileName,"_error.log",sep=""),append=TRUE)
        if (!is.null(res)){
          dbClearResult(res)  
        }
        res=NULL
        dbDisconnect(dbConnection)
      }
    )  
  }
}

queryDatabaseWithPeptideForAffinityScore=function(peptide){
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
        # suppress typecast warnings (fixed in RMySQL dev v0.11, but not in stable yet)
        res<-suppressWarnings(dbGetQuery(dbConnection, paste('SELECT peptide, ',paste(hlaTypes,"affinity",collapse=", ",sep=""),' ',
                                                             'FROM perPepAffinityScores_NEW ',
                                                             'WHERE peptide="',peptide,'"',";",
                                                             sep="")))
        dbDisconnect(dbConnection)
        return(as.data.table(res))
      }
      ,
      error=function(err){
        write(x=paste("Problem with query for affinityScore. File:",fileName," | index:",i," | ",err),file=paste(scriptPath,"/logs/",fileName,"_error.log",sep=""),append=TRUE)
        if (!is.null(res)){
          dbClearResult(res)  
        }
        res=NULL
        dbDisconnect(dbConnection)
      }
    )  
  }
}

queryDatabaseWithENSGForGeneIdAndENST=function(ENSG){
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
        res<-dbGetQuery(dbConnection, paste('SELECT ENSG, ENST, geneID ',
                                            'FROM idInfo ',
                                            'WHERE ENSG="',ENSG,'"',";",
                                            sep=""))
        dbDisconnect(dbConnection)
        return(as.data.table(res))
      }
      ,
      error=function(err){
        write(x=paste("Problem with query for geneID & ENST. File:",fileName," | index:",i," | ",err),file=paste(scriptPath,"/logs/",fileName,"_error.log",sep=""),append=TRUE)
        if (!is.null(res)){
          dbClearResult(res)  
        }
        res=NULL
        dbDisconnect(dbConnection)
      }
    )  
  }
}