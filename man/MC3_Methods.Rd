% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mc3_mthds.R
\name{MC3_Methods}
\alias{MC3_Methods}
\alias{mc3_mthds}
\title{List of level 3 multiple-concentration normalization methods}
\usage{
mc3_mthds()
}
\value{
A list of functions
}
\description{
\code{mc3_mthds} returns a list of normalization methods to be used
during level 3 multiple-concentration processing.
}
\details{
The functions contained in the list returned by \code{mc3_mthds} take
'aeids' (a numeric vector of aeid values) and returns a list of expressions
to be executed in the \code{mc3} (not exported) function environment. The
functions are described here for reference purposes, The
\code{mc3_mthds} function is not exported, nor is it intended for use.

All available methods are described in the Available Methods section, listed
by the type of function and the function/method name.
}
\note{
This function is not exported and is not intended to be used by the user.
}
\section{Available Methods}{


The methods are broken into three types, based on what fields they define.
Different methods are used to define "bval" (the baseline value), "pval"
(the positive control value), and "resp" (the final response value).

Although it does not say so specifically in each description, all methods
are applied by aeid.

More information about the level 3 multiple-concentration processing is
available in the package vignette, "Pipeline_Overview."

\subsection{bval Methods}{
  \describe{
    \item{bval.apid.nwlls.med}{Calculate bval as the median of cval for
    wells with wllt equal to "n," by apid.}
    \item{bval.apid.1owconc.med}{Calculate bval as the median of cval for
    wells with wllt equal to "t" and cndx equal to 1 or 2, by apid.}
    \item{bval.apid.twlls.med}{Calculate bval as the median of cval for
    wells with wllt equal to "t," by apid.}
    \item{bval.apid.tn.med}{Calculate bval as the median of cval for wells
    with wllt equal to "t" or "n," by apid.}
    \item{bval.apid.nwllslowconc.med}{Calculate bval as the median of cval
    for wells with wllt equal to "n" or wells with wllt equal to "t" and
    cndx equal to 1 or 2, by apid.}
    \item{bval.spid.lowconc.med}{Calculate bval as the median of cval for
    wells with wllt equal to "t" and cndx equal to 1, 2, or 3, by spid.}
  }
}

\subsection{pval Methods}{
  \describe{
    \item{pval.apid.pwlls.med}{Calculate pval as the median of cval for
    wells with wllt equal to "p," by apid.}
    \item{pval.apid.mwlls.med}{Calculate pval as the median of cval for
    wells with wllt equal to "m," by apid.}
    \item{pval.apid.medpcbyconc.max}{First calculate the median of cval for
    wells with wllt equal to "p" or "c," by wllt, conc, and apid. Then
    calculate pval as the maximum of the calculated medians, by apid.}
    \item{pval.apid.medpcbyconc.min}{First calculate the median of cval for
    wells with wllt equal to "p" or "c," by wllt, conc, and apid. Then
    calculate pval as the minimum of the calculated medians, by apid.}
    \item{pval.apid.medncbyconc.min}{First calculate the median of cval for
    wells with wllt equal to "m" or "o," by wllt, conc, and apid. Then
    calculate pval as the minimum of the calculated medians, by apid.}
    \item{pval.apid.pmv.min}{First calculate the median of cval for
    wells with wllt equal to "p," "m," or "v," by wllt, conc, and apid. Then
    calculate pval as the minimum of the calculated medians, by apid.}
    \item{pval.apid.pmv.max}{First calculate the median of cval for
    wells with wllt equal to "p," "m," or "v," by wllt, conc, and apid. Then
    calculate pval as the maximum of the calculated medians, by apid.}
    \item{pval.apid.f.max}{First calculate the median of cval for
    wells with wllt equal to "f," by wllt, conc, and apid. Then calculate
    pval as the maximum of the calculated medians, by apid.}
    \item{pval.apid.f.min}{First calculate the median of cval for
    wells with wllt equal to "f," by wllt, conc, and apid. Then calculate
    pval as the minimum of the calculated medians, by apid.}
    \item{pval.apid.p.max}{First calculate the median of cval for
    wells with wllt equal to "p," by wllt, conc, and apid. Then calculate
    pval as the maximum of the calculated medians, by apid.}
    \item{pval.apid.p.min}{First calculate the median of cval for
    wells with wllt equal to "p," by wllt, conc, and apid. Then calculate
    pval as the minimum of the calculated medians, by apid.}
    \item{pval.apid.v.min}{First calculate the median of cval for
    wells with wllt equal to "v," by wllt, conc, and apid. Then calculate
    pval as the minimum of the calculated medians, by apid.}
    \item{pval.zero}{Define pval as 0.}
  }
}

\subsection{resp Methods}{
  \describe{
    \item{resp.pc}{Calculate resp as \eqn{\frac{\mathit{cval} -
    \mathit{bval}}{\mathit{pval} - \mathit{bval}}100}{(cval - bval)/(pval
    - bval)*100}.}
    \item{resp.fc}{Calculate resp as \eqn{\mathit{cval}/\mathit{bval}}{
    cval/bval}.}
    \item{resp.logfc}{Calculate resp as \eqn{\mathit{cval} - \mathit{bval}}{
    cval - bval}.}
    \item{resp.log2}{Take the logarithm of resp with base 2.}
    \item{resp.mult25}{Multiply resp by 25.}
    \item{resp.scale.mad.log2fc}{Multiply resp by the scale factor
    \eqn{\frac{\log_2 \left( 1.2 \right)}{3\mathit{bmad}}}{
    log2(1.2)/(3*bmad)}.}
    \item{resp.scale.quant.log2fc}{Determine the maximum response
    \eqn{\mathit{md}}{md} where \eqn{\mathit{md}}{md} = abs(1st centile -
    50th centile) or abs(99th centile - 50th centile), whichever is greater.
    Scale the response such that 20 percent of md equals
    \eqn{\log_2 \left( 1.2 \right)}{log2(1.2)}.}
    \item{resp.multneg1}{Multiply resp by -1.}
    \item{resp.shiftneg.3bmad}{Shift all resp values less than 3*bmad to 0.}
    \item{resp.shiftneg.6bmad}{Shift all resp values less than 6*bmad to 0.}
    \item{resp.shiftneg.10bmad}{Shift all resp values less than 10*bmad to
    0.}
    \item{resp.blineshift.3bmad.repi}{Shift resp values with the blineShift
    function by repi, where the window (wndw) is 3*bmad.}
    \item{resp.blineshift.50.repi}{Shift resp values with the blineShift
    function by repi, where the window (wndw) is 50.}
    \item{resp.blineshift.3bmad.spid}{Shift resp values with the blineShift
    function by spid, where the window (wndw) is 3*bmad.}
    \item{resp.blineshift.50.spid}{Shift resp values with the blineShift
    function by spid, where the window (wndw) is 50.}
    \item{none}{Do no normalization; make resp equal to cval.}
  }
}
}

\seealso{
\code{\link{mc3}}, \code{\link{gtoxMthdLoad}} to query what methods
get applied to each aeid
}
\keyword{internal}
