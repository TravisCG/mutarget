# Useful functions
fetchDB <- function(con, query){
        rs  <- dbSendQuery(con, query)
        raw <- fetch(rs, n=-1)
	return(raw)
}

getExpMatrix <- function(con, cancerid, datasrc){
	query <- paste("select individual.name, genename, value from expression left join (individual, genetable) on (genetable_geneid = geneid and individual_patientid = patientid) where individual_cancer_cancerid = ",cancerid," and datasource_datasourceid = ",datasrc,";", sep = "");
	raw   <- fetchDB(con, query)

}

# Common initialisations
con  <- dbConnect(MySQL(), user="XXXX", password="XXXX", dbname="mutarget", host="localhost")
