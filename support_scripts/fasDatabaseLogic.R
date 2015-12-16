suppressPackageStartupMessages(library(RMySQL))

# set mysql configuration
sqlConfiguration = data.table(sqlhost = "medoid",
                              sqluser = "l.fanchi",
                              sqlpass = "MpRi1RKd",
                              sqldbname = "SchumiDB")

# wrapper function for logging query errors
logQueryErrorToDisk = function(querytype, file, index, error) {
  write(x = paste0(Sys.time(), " - Problem with query for ", querytype,
                   " | file: ", file,
                   " | index: ", index,
                   " | error: ", error),
        file = paste0(config$filepath, "/output/", file, "_query_errors.log"),
        append = TRUE)
}

performFasDbPredictions = function(index, peptides, peptidestretch, allele, peptidelength) {
  if(nrow(peptides) > 0) {
    # do affinity lookups in FASdb
    affinityLookups = queryDatabaseWithPeptideForAffinityScore(index = i,
                                                               peptides = peptides$peptide,
                                                               allele = allele)
    
    predictions = merge(x = peptides,
                        y = unique(x = affinityLookups,
                                   by = "peptide"),
                        by = "peptide",
                        all.x = TRUE)
    
    # perform affinity predictions on peptides missing affinity data
    if (any(is.na(predictions[[paste0(allele, "affinity")]]))){
      affinityPredictions = performAffinityPredictions(peptides = subset(x = predictions,
                                                                         subset = is.na(predictions[[paste0(allele, "affinity")]]) == TRUE)$peptide,
                                                       allele = allele,
                                                       peptidelength = peptidelength)  
    } else {
      affinityPredictions = emptyTableWithColumnNamesAndColumnClasses(colnames = c("hla_allele", "peptide", paste0(allele, "affinity")),
                                                                      colclasses = c("character", "character", "numeric"))
    }
    
    # perform processing predictions
    processingPredictions = performProcessingPredictions(peptidestretch = peptidestretch)
    
    # merge prediction info
    predictions[[paste0(allele,"affinity")]][match(x = unique(x = affinityPredictions,
                                                              by = "peptide")$peptide,
                                                   table = predictions$peptide)] = unique(x = affinityPredictions,
                                                                                          by = "peptide")[[paste0(allele,"affinity")]][match(x = predictions$peptide,
                                                                                                                                             table = unique(x = affinityPredictions,
                                                                                                                                                            by = "peptide")$peptide)]
    
    predictions = merge(x = predictions,
                        y = processingPredictions,
                        by = "c_term_pos")
    
    predictions[, hla_allele := gsub(pattern = "([A-C])([0-9]{2})([0-9]{2,3})",
                                     replacement = "HLA-\\1\\*\\2:\\3",
                                     x = allele)]
    predictions[, xmer := peptidelength]
  } else {
    predictions = emptyTableWithColumnNamesAndColumnClasses(colnames = c(names(peptides),
                                                                         "xmer", "hla_allele", paste0(allele, "affinity"), "c_term_aa", "processing_score"),
                                                            colclasses = c(unlist(lapply(peptides, class),
                                                                                  use.names=FALSE),
                                                                           "numeric", "character", "numeric", "character", "numeric"))
  }
  
  return(predictions)
}

# function for querying FASdb for peptide affinity
queryDatabaseWithPeptideForAffinityScore = function(index, peptides, allele) {
  res = NULL
  attempt = 1
  while(is.null(res) && attempt <= 10) {
    attempt = attempt + 1
    
    tryCatch(
      {
        dbConnection = dbConnect(MySQL(),
                                 host = sqlConfiguration$sqlhost,
                                 user = sqlConfiguration$sqluser,
                                 password = sqlConfiguration$sqlpass,
                                 dbname = sqlConfiguration$sqldbname)
        
        res = dbGetQuery(dbConnection, paste0('SELECT peptide, ',paste0(allele,"affinity"),' ',
                                              'FROM perPepAffinityScores_NEW',' ',
                                              'WHERE ', paste0('peptide="',peptides,'"',collapse = ' OR '), ';'))
        dbDisconnect(dbConnection)
        
        res = merge(x = data.table(peptide = peptides),
                    y = res,
                    by = "peptide",
                    all.x = TRUE)
        
        return(res)
      }
      ,
      error = function(err) {
        logQueryErrorToDisk(querytype = "affinityScore",
                            file = config$filename,
                            index = index,
                            error = err)
        
        if (!is.null(res)) {
          dbClearResult(res)  
        }
        
        res = NULL
        if(exists("dbConnection")) {
          dbDisconnect(dbConnection)  
        }
      }
    )
  }
}
