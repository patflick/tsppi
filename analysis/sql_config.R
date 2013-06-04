# Config for the SQL connection.
# 
# Author: flick
###############################################################################
DEFAULT_DATABASE_FILE = 'hpaDB.sqlite'




get_sql_conn <- function(database_file=DEFAULT_DATABASE_FILE)
{
	# backup working directory
	prior_wd <- getwd()
	
	# sets the path of the sqlite3 database file
	if (.Platform$OS.type == "unix")
	{
		#setwd("/home/patrick/dev/bio/data");
		setwd("/cygdrive/d/PPI")
	} else {
		setwd("D:\\PPI");
	}
	
	# adds the file to the path in a platform independent manner
	DATABASE_NAME = file.path(getwd(), database_file)
	
	# loads the sqlite3 database driver
	library('RSQLite')
	DATABASE_DRIVER <- dbDriver("SQLite")
	
	# connect
	con <- dbConnect(DATABASE_DRIVER, dbname = DATABASE_NAME)
	
	# restore working directory
	setwd(prior_wd)
	
	# return the connection object
	return(con)
}




