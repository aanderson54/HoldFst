##' Top Windows for Heterozygosity and Fixation  (Description)
##'
##' Input results from running HoldFst to get the top most significant windows (Details)
##' @title Find top windows of the genome associated with a difference between the normal and variant population
##' @param object The object containing the results from HoldFst
##' @param sig Which scores should be used to find the most significant windows ('Fst'=Fixation Index, 'Phet'=heterozygosity of P population,'Qhet'=heterozygosyity of Q population,'Diff'=difference in heterozygosity,'All'=all scores(default))
##' @author Ashlyn Anderson
##' @export
##' @examples
##' #topwindows(results,sig='Qhet')
##'



topwindows <- function(object, sig = "all") {
    ifelse(is.null(sig), "all", sig)
    if (sig != "Fst" & sig != "Qhet" & sig != "Phet" & sig != "Diff" & sig != "all")
        stop("Available sig options are ('Fst'=Fixation Index, 'Phet'=heterozygosity of P population,'Qhet'=heterozygosyity of Q population,'Diff'=difference in heterozygosity,'All'=all scores(default))")
    if (is.list(results)) {
        data = results[[1]]
        data$QhetP <- 2 * pnorm(-abs(data$Qz))
        data$PhetP <- 2 * pnorm(-abs(data$Pz))
        data$hetDiff <- data$Qhet - data$Phet
        data$zDiff <- (data$hetDiff - mean(data$hetDiff)/sd(data$hetDiff))
        data$pDiff <- 2 * pnorm(-abs(data$zDiff))
        data <- data[, c(1, 2, 3, 4, 5, 6, 7, 8, 11, 9, 10, 12, 13, 14, 15)]
        sigFst <- head(data[order(data$FstP), ], 20)
        sigPhet <- head(data[order(data$PhetP), ], 20)
        sigQhet <- head(data[order(data$QhetP), ], 20)
        sigDiff <- head(data[order(data$pDiff), ], 20)
        if (sig == "Fst")
            return(sigFst) else if (sig == "Qhet")
            return(sigQhet) else if (sig == "Phet")
            return(sigPhet) else if (sig == "Diff")
            return(sigDiff) else if (sig == "all") {
            a <- merge(sigFst, sigPhet, by = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
                13, 14, 15))
            b <- merge(a, sigQhet, by = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13,
                14, 15))
            c <- merge(b, sigQhet, by = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13,
                14, 15))
            return(c)
        }
    } else if (is.data.frame(object)) {
        data <- object
        data$QhetP <- 2 * pnorm(-abs(data$Qz))
        data$PhetP <- 2 * pnorm(-abs(data$Pz))
        data$hetDiff <- data$Qhet - data$Phet
        data$zDiff <- (data$hetDiff - mean(data$hetDiff)/sd(data$hetDiff))
        data$pDiff <- 2 * pnorm(-abs(data$zDiff))
        data <- data[, c(1, 2, 3, 4, 5, 6, 7, 8, 11, 9, 10, 12, 13, 14, 15)]
        sigFst <- head(data[order(data$FstP), ], 20)
        sigPhet <- head(data[order(data$PhetP), ], 20)
        sigQhet <- head(data[order(data$QhetP), ], 20)
        sigDiff <- head(data[order(data$pDiff), ], 20)
        if (sig == "Fst")
            return(sigFst) else if (sig == "Qhet")
            return(sigQhet) else if (sig == "Phet")
            return(sigPhet) else if (sig == "Diff")
            return(sigDiff) else if (sig == "all") {
            a <- merge(sigFst, sigPhet, by = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
                13, 14, 15))
            b <- merge(a, sigQhet, by = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13,
                14, 15))
            c <- merge(b, sigQhet, by = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13,
                14, 15))
            return(c)
        }
    }
}
