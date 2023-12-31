#####################################################################
## This program is distributed in the hope that it will be useful, ##
## but WITHOUT ANY WARRANTY; without even the implied warranty of  ##
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the    ##
## GNU General Public License for more details.                    ##
#####################################################################

#-------------------------------------------------------------------------------
# gtoxLoadChem: Load sample/chemical information
#-------------------------------------------------------------------------------

#' @title Load sample/chemical information
#'
#' @description
#' \code{gtoxLoadChem} queries the gtox database and returns the chemcial
#' information for the given field and values.
#'
#' @param field Character of length 1, the field to query on
#' @param val Vector of values to subset on
#' @param exact Logical, should chemical names be considered exact?
#' @param include.spid Logical, should spid be included?
#'
#' @details
#' The 'field' parameter is named differently from the 'fld' parameter seen
#' in other functions because it only takes one input.
#'
#' The functionality of the 'exact' parameter cannot be demonstrated within
#' the SQLite environment. However, in the MariaDB environment the user should
#' be able to give parital chemcial name strings, to find chemicals with
#' similar names. For example, setting 'val' to "phenol" when 'field' is "chnm"
#' and 'exact' is \code{FALSE} might pull up the chemicals "mercury". More
#' technically, setting 'exact' to \code{FALSE} passes the string in 'val' to
#' an RLIKE statement within the MariaDB query.
#'
#' @examples
#' ## Store the current config settings, so they can be reloaded at the end
#' ## of the examples
#' conf_store <- gtoxConfList()
#' gtoxConfDefault()
#'
#' ## Passing no parameters gives all of the registered chemicals with their
#' ## sample IDs
#' gtoxLoadChem()
#'
#' ## Or the user can exclude spid and get a unique list of chemicals
#' gtoxLoadChem(include.spid = FALSE)
#'
#' ## Other examples:
#' gtoxLoadChem(field = "chnm", val = "chromium")
#'
#' ## Reset configuration
#' options(conf_store)
#'
#' @return A data.table with the chemical information for the given parameters
#'
#' @import data.table
#' @export

gtoxLoadChem <- function(field=NULL, val=NULL, exact=TRUE,
    include.spid=TRUE) {

    ## Variable-binding to pass R CMD Check
    code <- casn <- chid <- chnm <- NULL

    if (!is.null(field)) {
        vfield <- c("chid", "spid", "chnm", "casn", "code")
        if (!field %in% vfield) stop("Invalid 'field' value.")
    }

    if (getOption("TCPL_DRVR") == "SQLite") {
        if (!exact) {
            warning("The exact=FALSE option is not supported in SQLite.")
        }
        exact <- TRUE
    }

    qstring <- .ChemQ(field=field, val=val, exact=exact)

    dat <- gtoxQuery(query=qstring)

    if (nrow(dat) == 0) {
        warning("The given ", field,"(s) are not in the gtox database.")
        return(dat[])
    }

    dat[ , code := NA_character_]
    dat[!is.na(casn) , code := paste0("C", gsub("-|_", "", casn))]
    dat[ , chid := as.integer(chid)]
    dat[, chem_desc := as.character(chem_desc)]
    dat[, add_info := as.character(add_info)]
    dat <- unique(dat)
    
    if (include.spid) return (dat[])

    dat <- unique(dat[ , spid := NULL])

    dat[]

}
