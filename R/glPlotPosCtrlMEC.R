#####################################################################
## This program is distributed in the hope that it will be useful, ##
## but WITHOUT ANY WARRANTY; without even the implied warranty of  ##
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the    ##
## GNU General Public License for more details.                    ##
#####################################################################

#-------------------------------------------------------------------------------
# glPlotPosCtrlMEC: boxplot to check positive controls vs historical data
#-------------------------------------------------------------------------------

#' @title Box plot for positive control check
#' @description This function plots positive controls for study id asid as well
#' as boxplot historical positive control MECs
#'
#' @param asid Assay source id
#' @param masked Masking color
#'
#' @details
#' This funtion is useful to select plates to mask
#'
#' @examples
#' ## Store the current config settings, so they can be reloaded at the end
#' ## of the examples
#' conf_store <- gtoxConfList()
#' gtoxConfDefault()
#'
#' ## Create boxplot for all endpoints and chemicals tested. Useful to save
#' ## plots in a pdf file.
#' pp <- glPlotPosCtrlMEC(asid = 1L)
#' pp[[1]]
#'
#' ## Reset configuration
#' options(conf_store)
#'
#' @return A list of ggplot objects, one per assay X timepoint.
#'
#' @note PMI-specific
#' @import ggplot2
#' @export
#'

glPlotPosCtrlMEC <- function(asid, masked=NULL) {

    addFlds <- c("asid", "aid", "anm", "acnm")
    
    # get data for study asid
    acids = gtoxLoadAcid(fld="asid", val=asid)$acid
    lvl0 = unique(gtoxLoadData(lvl=0, fld=c("wllt", "acid"), val=list("c", acids))[, c("spid", "apid")])
    lvl0 = merge(lvl0, gtoxLoadApid(fld="apid", val=list(lvl0$apid)), by="apid")[, c("spid", "u_boxtrack")]
    aetbl <- gtoxLoadAeid(fld="asid", val=asid, add.fld=addFlds)
    aetbl <- aetbl[!grepl("Cell count", aetbl$aenm), ]
    aetbl <- aetbl[!grepl("Cytotoxicity (TIER1)_Cell membrane permeability", aetbl$aenm, fixed=TRUE), ]
    dat <- gtoxPrepOtpt(gtoxLoadData(lvl=5, fld=c("aeid", "spid"), val=list(aetbl$aeid, lvl0$spid)))
    dat <- merge(dat, aetbl, by=c("aeid", "aenm"))
    dat <- merge(dat, lvl0, by="spid")
    dat <- dat[ , .SD[which.min(modl_rmse)], by=c("spid", "acnm")]
    dat[ , aenm := sapply(strsplit(as.character(aenm), "_"), "[[", 2)]
    dat$aenm <- gsub(" changes", "", dat$aenm)
    dat$aenm_wrap <- str_wrap(string=dat$aenm, width=15)
    dat$expTm = getsplit(dat$anm, "_", last=TRUE)
    dat = dat[with(dat, order(aeid, u_boxtrack)), ]
    dat$replicate = unlist(lapply(table(dat$aeid), function(xx){1:xx}))
    
    # with no MEC
    dat2 = dat[is.na(dat$modl_acc), ]
    dat2$modl_acc = log10(470 - 20*dat2$replicate)
    dat2$masked = "yellow"
    
    # masked
    dat$masked = "lightgray"
    if(!is.null(masked)){
        dat$masked[dat$u_boxtrack%in%masked] = "tomato"
        dat2$masked[dat2$u_boxtrack%in%masked] = "tomato"
    }
    
    # historical data
    spids0 = unique(gtoxLoadData(lvl=0, fld="wllt", val="c")$spid)
    spids0 = spids0[!spids0%in%lvl0$spid]
    aetbl0 <- gtoxLoadAeid(add.fld=addFlds)
    aetbl0 <- aetbl0[!grepl("Cell count", aetbl0$aenm), ]
    aetbl0 <- aetbl0[!grepl("Cytotoxicity (TIER1)_Cell membrane permeability", aetbl0$aenm, fixed=TRUE), ]
    dat0 <- gtoxPrepOtpt(gtoxLoadData(lvl=5, fld=c("aeid", "spid"), 
                                      val=list(aetbl0$aeid, spids0)))
    dat0 <- merge(dat0, aetbl0, by=c("aeid", "aenm"))
    dat0 <- dat0[ , .SD[which.min(modl_rmse)], by=c("spid", "acnm")]
    dat0$aenm <- gsub(" changes", "", dat0$aenm)
    dat0[ , aenm := sapply(strsplit(as.character(aenm), "_"), "[[", 2)]
    dat0$aenm_wrap <- str_wrap(string=dat0$aenm, width=15)
    dat0$expTm = getsplit(dat0$anm, "_", last=TRUE)
    
    # create plot
    pp = ggplot(dat, aes(x=aenm_wrap, y=10^modl_acc)) +
            geom_boxplot(data = dat0, aes(x=aenm_wrap, y=10^modl_acc), 
                         alpha=0.5, lwd=0.3, col="darkgray",
                         outlier.shape = NA) +
            geom_point(shape=dat$replicate, size=2) +
            geom_label_repel(aes(label=u_boxtrack), size=2, data=dat2, fill=dat2$masked,
                             fontface="plain", seed=4778) +
            geom_label_repel(aes(label=u_boxtrack), size=2, data=dat, fill=dat$masked,
                             fontface="plain", seed=4778) +
            scale_y_continuous(name="Minimal Effective Concentration") +
            xlab("") +
            facet_grid(expTm~.)
    return(pp)
}

