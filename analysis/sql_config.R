# Config for the SQL connection.
# 
# Author: flick
###############################################################################
DATABASE_FILE = 'hpaDB.sqlite'

# backup working directory
prior_wd <- getwd()

# sets the path of the sqlite3 database file
if (.Platform$OS.type == "unix")
{
	setwd("/home/patrick/dev/bio/data");
} else {
	setwd("D:\\PPI");
}

# adds the file to the path in a platform independent manner
DATABASE_NAME = file.path(getwd(), DATABASE_FILE)

# loads the sqlite3 database driver
library('RSQLite')
DATABASE_DRIVER <- dbDriver("SQLite")


get_sql_conn <- function()
{
	# connect
	con <- dbConnect(DATABASE_DRIVER, dbname = DATABASE_NAME)
	return(con)
}



# restore working directory
setwd(prior_wd)
