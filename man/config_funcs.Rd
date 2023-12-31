% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/config_funcs.R, R/gtoxConf.R,
%   R/gtoxConfDefault.R, R/gtoxConfList.R, R/gtoxConfLoad.R, R/gtoxConfReset.R,
%   R/gtoxConfSave.R
\name{Configure functions}
\alias{Configure functions}
\alias{gtoxConf}
\alias{gtoxConfDefault}
\alias{gtoxConfList}
\alias{gtoxConfLoad}
\alias{gtoxConfReset}
\alias{gtoxConfSave}
\title{Functions for configuring the gtox package}
\usage{
gtoxConf(drvr = NULL, user = NULL, pass = NULL, host = NULL,
  db = NULL)

gtoxConfDefault()

gtoxConfList(show.pass = FALSE)

gtoxConfLoad(list.new = TRUE)

gtoxConfReset()

gtoxConfSave()
}
\arguments{
\item{drvr}{Character of length 1, which database driver to use}

\item{user}{Character of length 1, the database server username}

\item{pass}{Character of length 1, the database server password}

\item{host}{Character of length 1, the database server}

\item{db}{Character of length 1, the name of the gtox database}

\item{show.pass}{Logical, should the password be returned}

\item{list.new}{Logical of length 1, should the new settings be printed?}
}
\value{
None

\code{gtoxConf} changes \code{options} to set the gtox-specific options, 
most importantly to configure the connection to the gtox databases. 
\code{gtoxConf} will only change non-null values, and can be used to 
change a single value if needed. 

\code{gtoxConfSave} modifies the configuration file to reflect the current
gtox settings.

\code{gtoxConfList} lists the values assigned to the gtox global options.

\code{gtoxConfLoad} updates the gtox settings to reflect the current 
configuration file.

\code{gtoxConfDefault} changes the \code{options} to reflect the default
settings for the example SQLite database, but does not alter the 
configuration file.

\code{gtoxConfReset} is used to generate the initial configuration script,
and can be used to reset or regenerate the configuration script by the user.
}
\description{
These functions are used to configure the gtox settings.

Load the current configuration file
}
\details{
Currently, the gtox package only supports the "MariaDB" and "SQLite" database
drivers.

The settings can be stored in a configuration file to make the using the 
package more user-friendly. To create the configuration file, the user must
first create a system environment variable ('TCPL_CONF') that points to to 
the file. There is more information about system environment variables in
\code{\link{Startup}} and \code{\link{Sys.getenv}}. Briefly, the user 
needs to modify the '.Renviron' file in their home directory. If the file
does not exist, create it, and add the following line:

TCPL_CONF=path/to/confFile.conf

Here 'path/to/confFile.conf' can be any path to a file. One suggestion would
be to include gtoxConf in the home directory, eg. TCPL_CONF=~/gtoxConf. 
Note, '~' may not indicate the home directory on every operating system.
Once the environment variable is added, the user can change the settings
using \code{gtoxConf}, then save the settings to the file given by the
TCPL_CONF environment variable running \code{gtoxConfSave()}.
}
\examples{
gtoxConfList() # List configuration parameters


## Configure database
sqlite <- file.path(system.file(package="GladiaTOX"),
    "sql",
    "gladiatoxdb.sqlite")
gtoxConf(db=sqlite, user=NA, host=NA, drvr="SQLite")


## Configure database with default parameters
gtoxConfDefault()


## List configuration of database parameters
gtoxConfList()


## Set the environment variable pointing to the configuration file
Sys.setenv(TCPL_CONF=file.path(system.file(package="GladiaTOX"),"gtoxConf"))

## Configure database
gtoxConfLoad()


## Set the environment variable pointing to the configuration file
Sys.setenv(TCPL_CONF=file.path(system.file(package="GladiaTOX"),"gtoxConf"))

## Configure database
gtoxConfReset()


## Set the environment variable pointing to the configuration file
Sys.setenv(TCPL_CONF=file.path(system.file(package="GladiaTOX"),"gtoxConf"))

## Configure database
gtoxConfSave()

}
