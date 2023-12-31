#####################################################################
## This program is distributed in the hope that it will be useful, ##
## but WITHOUT ANY WARRANTY; without even the implied warranty of  ##
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the    ##
## GNU General Public License for more details.                    ##
#####################################################################

#-------------------------------------------------------------------------------
# gtoxSendQuery: Send query to the gtox databases
#-------------------------------------------------------------------------------

#' @rdname query_funcs
#'
#' @examples
#' 
#' ## Perform query
#' gtoxSendQuery(paste0("SELECT * FROM assay_source"))
#'
#' @import DBI
#' @importFrom RSQLite SQLite
#' @import data.table
#' @importFrom RMariaDB MariaDB
#' @importMethodsFrom RMariaDB dbSendQuery dbClearResult dbDisconnect dbConnect
#' @importFrom methods is
#' @export

gtoxSendQuery <- function(query, db=getOption("TCPL_DB"),
    drvr=getOption("TCPL_DRVR")) {

    ## Check for valid inputs
    if (length(query) != 1 | !is(query, "character")) {
        stop("The input 'query' must be a character of length one.")
    }
    if (length(db) != 1 | !is(db, "character")) {
        stop("The input 'db' must be a character of length one.")
    }

    db_pars <- NULL

    if (getOption("TCPL_DRVR") == "SQLite") {

        db_pars <- list(drv=SQLite(),
                        dbname=db)

    }

    if (getOption("TCPL_DRVR") == "MariaDB") {

        if (any(is.na(options()[c("TCPL_USER", "TCPL_HOST", "TCPL_PASS")]))) {
            stop(
                "Must configure TCPL_USER, TCPL_HOST, and TCPL_PASS options. ",
                "See ?gtoxConf for more details.")
        }

        db_pars <- list(drv=MariaDB(),
                        user=getOption("TCPL_USER"),
                        password=getOption("TCPL_PASS"),
                        host=getOption("TCPL_HOST"),
                        dbname=db)

    }

    if (is.null(db_pars)) {

        stop(
            getOption("TCPL_DRVR"), " is not a supported database system. See ",
            "?gtoxConf for more details.")

    }

    dbcon <- do.call(dbConnect, db_pars)
    temp <- try(dbSendQuery(dbcon, query), silent=TRUE)
    if (!is(temp, "try-error")) dbClearResult(temp)
    dbDisconnect(dbcon)

    if (!is(temp, "try-error")) return(TRUE)

    temp[1]

}

#-------------------------------------------------------------------------------
