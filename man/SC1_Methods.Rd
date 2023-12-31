% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sc1_mthds.R
\name{SC1_Methods}
\alias{SC1_Methods}
\alias{sc1_mthds}
\title{List of level 1 single-concentration normalization functions}
\usage{
sc1_mthds()
}
\value{
A list functions
}
\description{
\code{sc1_mthds} returns a list of functions to be used during level 1
single-concentration processing.
}
\details{
The functions contained in the list returned by \code{sc1_mthds} return
a list of expressions to be executed in the \code{sc2} (not exported)
function environment. The functions are described here for reference
purposes, The \code{sc1_mthds} function is not exported, nor is it
intended for use.

All available methods are described in the Available Methods section, listed
by the function/method name.
}
\note{
This function is not exported and is not intended to be used by the user.
}
\section{Available Methods}{


The methods are broken into three types, based on what fields they define.
Different methods are used to define "bval" (the baseline value), "pval"
(the positive control value), and "resp" (the final response value).

Although it does not say so specifically in each description, all methods
are applied by acid.

More information about the level 3 single-concentration processing is
available in the package vignette, "Pipeline_Overview."

\subsection{bval Methods}{
  \describe{
    \item{bval.apid.nwlls.med}{Calculate bval as the median of rval for
    wells with wllt equal to "n," by apid.}
    \item{bval.apid.twlls.med}{Calculate bval as the median of rval for
    wells with wllt equal to "t," by apid.}
    \item{bval.apid.tn.med}{Calculate bval as the median of rval for wells
    with wllt equal to "t" or "n," by apid.}
  }
}

\subsection{pval Methods}{
  \describe{
    \item{pval.apid.pwlls.med}{Calculate pval as the median of rval for
    wells with wllt equal to "p," by apid.}
    \item{pval.apid.mwlls.med}{Calculate pval as the median of rval for
    wells with wllt equal to "m," by apid.}
    \item{pval.apid.medpcbyconc.max}{First calculate the median of rval for
    wells with wllt equal to "p" or "c," by wllt, conc, and apid. Then
    calculate pval as the maximum of the calculated medians, by apid.}
    \item{pval.apid.medpcbyconc.min}{First calculate the median of rval for
    wells with wllt equal to "p" or "c," by wllt, conc, and apid. Then
    calculate pval as the minimum of the calculated medians, by apid.}
    \item{pval.apid.medncbyconc.min}{First calculate the median of rval for
    wells with wllt equal to "m" or "o," by wllt, conc, and apid. Then
    calculate pval as the minimum of the calculated medians, by apid.}
    \item{pval.zero}{Define pval as 0.}
  }
}

\subsection{resp Methods}{
  \describe{
    \item{resp.pc}{Calculate resp as \eqn{\frac{\mathit{rval} -
    \mathit{bval}}{\mathit{pval} - \mathit{bval}}100}{(rval - bval)/(pval
    - bval)*100}.}
    \item{resp.fc}{Calculate resp as \eqn{\mathit{rval}/\mathit{bval}}{
    rval/bval}.}
    \item{resp.logfc}{Calculate resp as \eqn{\mathit{rval} - \mathit{bval}}{
    rval - bval}.}
    \item{resp.log2}{Take the logarithm of resp with base 2.}
    \item{resp.multneg1}{Multiply resp by -1.}
    \item{none}{Do no normalization; make resp equal to rval.}
  }
}
}

\seealso{
\code{\link{sc1}}, \code{\link{Method functions}} to query what
methods get applied to each acid
}
\keyword{internal}
