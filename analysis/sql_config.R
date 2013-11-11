# Config for the SQL connection.
# 
# Author: flick
###############################################################################
DEFAULT_DATABASE_FILE = 'hpaDB.sqlite'


get_sql_conn <- function(database_file=DEFAULT_DATABASE_FILE)
{

    # adds the file to the path in a platform independent manner
    #DATABASE_NAME = file.path(getwd(), database_file)
    DATABASE_NAME = database_file

    # loads the sqlite3 database driver
    library('RSQLite')
    DATABASE_DRIVER <- dbDriver("SQLite")

    # connect
    con <- dbConnect(DATABASE_DRIVER, dbname = DATABASE_NAME)

    # return the connection object
    return(con)
}




