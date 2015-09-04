## this function allows writing to a table in a MySQL database and updating of existing records
dbUpdateRecord=function(dbtable, data=NULL, primary, vars) {
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
          
          if (!dbExistsTable(dbConnection, dbtable)) {
            stop("The target table \"", dbtable, "\" doesn't exist in the database \"", dbGetInfo(dbConnection)$dbname, "\"\n\n", call. = FALSE)
          }
          if (is.null(data)) {
            stop("The source dataframe is missing, with no default\n\n", call. = FALSE)
          }
          if (!(toupper(primary) %in% toupper(names(data)))) {
            stop("The primary key variable doesn't exist in the source dataframe\n\n", call. = FALSE)
          }
          if (!all(toupper(vars) %in% toupper(names(data)))) {
            stop("One or more variables don't exist in the source dataframe\n\n", call. = FALSE)
          }
          if (!(toupper(primary) %in% toupper(dbListFields(dbConnection, dbtable)))) {
            stop("The primary key variable doesn't exist in the target table\n\n", call. = FALSE)
          }
          if (!all(toupper(vars) %in% toupper(dbListFields(dbConnection, dbtable)))) {
            stop("One or more variables don't exist in the target table\n\n", call. = FALSE)
          }
          
          if (length(vars)>1){
            pastedvars=paste("'",apply(data[,vars,with=FALSE],1,paste,collapse="','"), "'",sep="")
          } else{
            pastedvars=paste("'",data[,vars,with=FALSE], "'",sep="")
          }
          
          varlist = paste(dbtable, "(", paste(c(primary, vars), collapse=", "), ")", sep="")
          datastring = paste("('", paste(paste(data[, primary, with=FALSE], pastedvars, sep="', "), collapse="), ("), ")", sep="")
          toupdate = paste(paste(vars, "=VALUES(", vars, ")", sep=""), collapse=", ")                          
          
          sqlstring = paste("INSERT INTO", varlist, "VALUES", datastring, "ON DUPLICATE KEY UPDATE", toupdate,";")
          res<-dbSendQuery(dbConnection, sqlstring)
          dbClearResult(res)
          
          dbDisconnect(dbConnection)
        }
        ,
        error=function(err)
        {
          write(x=paste0(format(Sys.Date(),"%Y%m%d")," - Problem with logging progress | file: ",fileName," | index: ",i," | error: ",err),
            file=paste(scriptPath,"/logs/",fileName,"_logProgress_errors.log",sep=""),
            append=TRUE)

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