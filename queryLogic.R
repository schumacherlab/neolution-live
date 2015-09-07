## Error logging
logQueryErrorToDisk=function(querytype,error){
  write(x=paste0(format(Sys.Date(),"%Y%m%d")," - Problem with query for ",querytype," | file: ",fileName," | index: ",i," | error: ",error),
        file=paste(scriptPath,"/logs/",fileName,"_query_errors.log",sep=""),
        append=TRUE)
}

## Query functions
queryDatabaseWithChromIdAndLocation=function(chromID,chromLoc){
  res=NULL
  attempt=1
  while(is.null(res) && attempt<=10){
    attempt=attempt+1

    tryCatch(
      {
        dbConnection=dbConnect(MySQL(),
                               host=sqlhost,
                               user=sqluser,
                               password=sqlpass,
                               dbname=sqldbname)
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
        logQueryErrorToDisk("SNV info",err)

        if (!is.null(res)){
          dbClearResult(res)
        }

        res=NULL
        if(exists("dbConnection")){
          dbDisconnect(dbConnection)  
        }
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
                               host=sqlhost,
                               user=sqluser,
                               password=sqlpass,
                               dbname=sqldbname)
        res<-dbGetQuery(dbConnection, paste('SELECT ENSG, Strand ',
                                            'FROM ENSGtoStrand ',
                                            'WHERE ENSG="',ENSG,'"',";",
                                            sep=""))
        dbDisconnect(dbConnection)
        return(as.data.table(res))
      }
      ,
      error=function(err){
        logQueryErrorToDisk("strandInfo",err)

        if (!is.null(res)){
          dbClearResult(res)  
        }

        res=NULL
        if(exists("dbConnection")){
          dbDisconnect(dbConnection)  
        }
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
                               host=sqlhost,
                               user=sqluser,
                               password=sqlpass,
                               dbname=sqldbname)
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
        logQueryErrorToDisk("input mutation context (9mers)",err)

        if (!is.null(res)){
          dbClearResult(res)  
        }

        res=NULL
        if(exists("dbConnection")){
          dbDisconnect(dbConnection)  
        }
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
                               host=sqlhost,
                               user=sqluser,
                               password=sqlpass,
                               dbname=sqldbname)
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
        logQueryErrorToDisk("peptide context for processing predictions",err)

        if (!is.null(res)){
          dbClearResult(res)  
        }

        res=NULL
        if(exists("dbConnection")){
          dbDisconnect(dbConnection)  
        }
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
                               host=sqlhost,
                               user=sqluser,
                               password=sqlpass,
                               dbname=sqldbname)
        res<-dbGetQuery(dbConnection, paste('SELECT codon, AA ',
                                            'FROM codonToAAC ',
                                            'WHERE codon="',codon,'"',";",
                                            sep=""))
        dbDisconnect(dbConnection)
        return(as.data.table(res))
      }
      ,
      error=function(err){
        logQueryErrorToDisk("codon to aac conversion",err)

        if (!is.null(res)){
          dbClearResult(res)  
        }

        res=NULL
        if(exists("dbConnection")){
          dbDisconnect(dbConnection)  
        }
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
                               host=sqlhost,
                               user=sqluser,
                               password=sqlpass,
                               dbname=sqldbname)
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
        logQueryErrorToDisk("processingScore",err)

        if (!is.null(res)){
          dbClearResult(res)  
        }

        res=NULL
        if(exists("dbConnection")){
          dbDisconnect(dbConnection)  
        }
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
                               host=sqlhost,
                               user=sqluser,
                               password=sqlpass,
                               dbname=sqldbname)
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
        logQueryErrorToDisk("affinityScore",err)

        if (!is.null(res)){
          dbClearResult(res)  
        }

        res=NULL
        if(exists("dbConnection")){
          dbDisconnect(dbConnection)  
        }
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
                               host=sqlhost,
                               user=sqluser,
                               password=sqlpass,
                               dbname=sqldbname)
        res<-dbGetQuery(dbConnection, paste('SELECT ENSG, ENST, geneID ',
                                            'FROM idInfo ',
                                            'WHERE ENSG="',ENSG,'"',";",
                                            sep=""))
        dbDisconnect(dbConnection)
        return(as.data.table(res))
      }
      ,
      error=function(err){
        logQueryErrorToDisk("geneSymbol & ENST info",err)

        if (!is.null(res)){
          dbClearResult(res)
        }

        res=NULL
        if(exists("dbConnection")){
          dbDisconnect(dbConnection)  
        }
      }
      )
  }
}
