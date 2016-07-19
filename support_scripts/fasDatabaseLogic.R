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
        file = file.path(runParameters$filepath, 'predictions_logs', paste0(file, '_query_errors.log')),
        append = TRUE)
}

performFasDbPredictions = function(index, peptides, peptidestretch, allele, peptidelength, predictor) {
  if (nrow(peptides) > 0) {
    # do affinity lookups in FASdb
    affinityLookups = queryDatabaseWithPeptideForAffinityScore(index = index,
                                                               peptides = peptides$peptide,
                                                               allele = allele,
                                                               predictor = predictor)

    predictions = merge(x = peptides,
                        y = unique(x = affinityLookups,
                                   by = "peptide"),
                        by = "peptide",
                        all.x = TRUE)

    # perform affinity predictions on peptides missing affinity data
    if (any(is.na(predictions[[paste0(allele, "affinity")]]))) {
      affinityPredictions = performAffinityPredictions(peptides = subset(x = predictions,
                                                                         subset = is.na(predictions[[paste0(allele, "affinity")]]) == TRUE)$peptide,
                                                       allele = allele,
                                                       peptidelength = peptidelength)

      # merge missing affinity prediction info
      predictions = merge(x = predictions,
                          y = unique(x = affinityPredictions,
                                     by = "peptide"),
                          by = "peptide",
                          all.x = TRUE,
                          suffixes = c(".lookups", ".predictions"))

      predictions[, paste0(allele, "affinity") := ifelse(test = is.na(predictions[[paste0(allele, "affinity.lookups")]]),
                                                         yes = predictions[[paste0(allele, "affinity.predictions")]],
                                                         no = predictions[[paste0(allele, "affinity.lookups")]])]
      predictions[, c(paste0(allele, "affinity.lookups"), paste0(allele, "affinity.predictions")) := NULL]

      predictions[, paste0(allele, "percentile_rank") := ifelse(test = is.na(predictions[[paste0(allele, "percentile_rank.lookups")]]),
                                                                yes = predictions[[paste0(allele, "percentile_rank.predictions")]],
                                                                no = predictions[[paste0(allele, "percentile_rank.lookups")]])]
      predictions[, c(paste0(allele, "percentile_rank.lookups"), paste0(allele, "percentile_rank.predictions")) := NULL]
    }

    # perform processing predictions
    processingPredictions = performProcessingPredictions(peptidestretch = peptidestretch)

    # merge prediction info
    predictions = merge(x = predictions,
                        y = processingPredictions,
                        by = "c_term_pos")

    predictions[, hla_allele := gsub(pattern = "([A-C])([0-9]{2})([0-9]{2,3})",
                                     replacement = "HLA-\\1\\*\\2:\\3",
                                     x = allele)]
    predictions[, xmer := peptidelength]
  } else {
    predictions = emptyTableWithColumnNamesAndColumnClasses(colnames = c(names(peptides),
                                                                         "xmer", "hla_allele", paste0(allele, "affinity"), paste0(allele, "percentile_rank"), "c_term_aa", "processing_score"),
                                                            colclasses = c(unlist(lapply(peptides, class),
                                                                                  use.names = FALSE),
                                                                           "numeric", "character", "numeric", "numeric", "character", "numeric"))
  }

  return(predictions)
}

# function for querying FASdb for peptide affinity
queryDatabaseWithPeptideForAffinityScore = function(index, peptides, allele, predictor) {
  res = NULL
  attempt = 1

  tableName = paste('binding', paste('pan', gsub(pattern = '.',
                                                 replacement = '_',
                                                 x = predictor,
                                                 fixed = T),
                                     sep = '_'),
                    sep = '_')

  while (is.null(res) && attempt <= 10) {
    attempt = attempt + 1

    tryCatch(
      {
        dbConnection = dbConnect(MySQL(),
                                 host = sqlConfiguration$sqlhost,
                                 user = sqlConfiguration$sqluser,
                                 password = sqlConfiguration$sqlpass,
                                 dbname = sqlConfiguration$sqldbname)

        if (paste0(allele, 'affinity') %in% dbListFields(conn = dbConnection,
                                                         name = tableName)) {
          if (predictor != '2.4') {
            res = dbGetQuery(dbConnection, paste('SELECT peptide,', paste0(allele, 'affinity'), ',', paste0(allele, 'percentile_rank'),
                                                 'FROM', tableName,
                                                 'WHERE', paste0('peptide="', peptides, '"',
                                                                 collapse = ' OR '),
                                                 ';'))
          } else if (predictor == '2.4') {
            res = dbGetQuery(dbConnection, paste('SELECT peptide,', paste0(allele, 'affinity'),
                                                 'FROM', tableName,
                                                 'WHERE', paste0('peptide="', peptides, '"',
                                                                 collapse = ' OR '),
                                                 ';'))
            res[[paste0(allele, 'percentile_rank')]] = NA
          }
        } else if (predictor != '2.4') {
          res = emptyTableWithColumnNamesAndColumnClasses(colnames = c('peptide', paste0(allele, 'affinity'), paste0(allele, 'percentile_rank')),
                                                          colclasses = c('character', 'numeric', 'numeric'))
        } else {
          res = emptyTableWithColumnNamesAndColumnClasses(colnames = c('peptide', paste0(allele, 'affinity')),
                                                          colclasses = c('character', 'numeric'))
        }

        dbDisconnect(dbConnection)

        res = merge(x = data.table(peptide = peptides),
                    y = res,
                    by = "peptide",
                    all.x = TRUE)

        return(res)
      }
      ,
      error = function(err) {
        logQueryErrorToDisk(querytype = "fetchAffinityScore",
                            file = runParameters$filename,
                            index = index,
                            error = err)

        if (!is.null(res)) {
          dbClearResult(res)
        }

        res = NULL
        if (exists("dbConnection")) {
          dbDisconnect(dbConnection)
        }
      })
  }
}

writePeptideAffinityToDatabase = function(index, allele, predictions, predictor) {
  attempt = 1

  tableName = paste('binding', paste('pan', gsub(pattern = '.',
                                                 replacement = '_',
                                                 x = predictor,
                                                 fixed = T),
                                     sep = '_'),
                    sep = '_')

  while (attempt <= 10) {
    attempt = attempt + 1

    tryCatch(
      {
        dbConnection = dbConnect(MySQL(),
                                 host = sqlConfiguration$sqlhost,
                                 user = sqlConfiguration$sqluser,
                                 password = sqlConfiguration$sqlpass,
                                 dbname = sqlConfiguration$sqldbname)

        fieldTypes = list(peptide = 'VARCHAR(15)')
        fieldTypes[[paste0(allele, 'affinity')]] = 'DOUBLE(16,4)'
        fieldTypes[[paste0(allele, 'percentile_rank')]] = 'DOUBLE(6,2)'

        if (paste0(allele, 'affinity') %ni% dbListFields(conn = dbConnection,
                                                         name = tableName)) {
          dbSendQuery(conn = dbConnection,
                      statement = paste('ALTER TABLE', tableName, 'ADD COLUMN', names(fieldTypes)[2], fieldTypes[2], ';'))
          dbSendQuery(conn = dbConnection,
                      statement = paste('ALTER TABLE', tableName, 'ADD COLUMN', names(fieldTypes)[3], fieldTypes[3], ';'))
        }

        dbSendQuery(conn = dbConnection,
                    statement = paste('INSERT INTO', tableName, paste0('(peptide,', allele, 'affinity,', allele, 'percentile_rank)'),
                                      'VALUES', paste0('("',paste(predictions$peptide, collapse = '","'), '"),(', paste(predictions[[paste0(allele, "affinity")]], collapse = ","), '),(', paste(predictions[[paste0(allele, 'percentile_rank')]], collapse = ","), ')'),
                                      'ON DUPLICATE KEY UPDATE', paste0(allele, "affinity", '=VALUES(', paste0(allele, "affinity"), '), ', paste0(allele, "percentile_rank"), '=VALUES(', paste0(allele, "percentile_rank"), ')'), ';'))

        # dbWriteTable(conn = dbConnection,
        #              name = tableName,
        #              value = predictions,
        #              field.types = fieldTypes,
        #              row.names = FALSE,
        #              overwrite = FALSE,
        #              append = TRUE)

        dbDisconnect(dbConnection)
      }
      ,
      error = function(err) {
        logQueryErrorToDisk(querytype = "writeAffinityScore",
                            file = runParameters$filename,
                            index = index,
                            error = err)

        if (!is.null(res)) {
          dbClearResult(res)
        }

        res = NULL
        if (exists("dbConnection")) {
          dbDisconnect(dbConnection)
        }
      })
  }
}
