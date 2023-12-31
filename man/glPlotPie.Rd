% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/glPlotPie.R
\name{glPlotPie}
\alias{glPlotPie}
\title{Pie chart for Minimal Effective Concentrations (MEC) and AC50 plot}
\usage{
glPlotPie(asid, chnms = NULL, acids = NULL, aeids = NULL,
  expos.time.ordr = NULL, stat = quote(modl_acc))
}
\arguments{
\item{asid}{Assay source id}

\item{chnms}{Character vector with list of chemical names}

\item{acids}{Numeric vector with list of acids}

\item{aeids}{Character vector with list of assay endpoints IDs}

\item{expos.time.ordr}{Character vector with sorted list of exposure times}

\item{stat}{Statistic to plot (e.g. MEC:modl_acc or modl_acb, AC50:modl_ga)}
}
\value{
None
}
\description{
This function plots MEC values
}
\details{
This funtion is useful to plot MEC or AC50 values
}
\examples{

## Create a pie plot of MEC values for all chemicals tested in the study
glPlotPie(asid=1L)

}
