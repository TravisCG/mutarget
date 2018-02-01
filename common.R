# Useful functions

# fetchDB: retreive data from the database using a query
fetchDB <- function(con, query){
        rs  <- dbSendQuery(con, query)
        raw <- fetch(rs, n=-1)
	return(raw)
}

# getExpMatrix: create an expression matrix for one cancer type and from one source
getExpMatrix <- function(con, cancerid, datasrc){
	query <- paste("select submitid, genename, value from expression inner join (individual, genetable) on (genetable_geneid = geneid and individual_patientid = patientid) where individual_cancer_cancerid = ",cancerid," and datasource_datasourceid = ",datasrc,";", sep = "");
	raw   <- fetchDB(con, query)
	count <- xtabs(value~genename+submitid, data = raw)
	count <- as.data.frame.matrix(count)
	count <- as.matrix(count)
	return(count)
}

# Common initialisations
con  <- dbConnect(MySQL(), user="XXXX", password="XXXX", dbname="mutarget", host="localhost")
