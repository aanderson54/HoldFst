##' Calculate Heterozygosity scores and Fixation Index for low-coverage whole genome sequencing (Description)
##'
##' Calculate and plot scores from low-coverage sequencing from a vcf input file (Details)
##' @title Heterozygosity and Fixation Index Scores
##' @param file file path to vcf file containing both the normal and variant population
##' @param P The name of the normal population (ex. 'Normal')
##' @param Q The name of the population with the phenotypic difference (ex. 'Variant')
##' @param bin Size of windows to calculate Heterozygosity and Fst scores over (default=500kb)
##' @param plot Logical. Should the function return manhattan plots (T/F)
##' @author Ashlyn Anderson
##' @export
##' @examples
##' #HoldFst('input.vcf','Normal','Variant',500000,plot=T)
##' @import ggplot2
##' @import vcfR
##' @import gridExtra

HoldFst <- function(file, P, Q, bin = 5e+05, plot = TRUE) {
    if (bin < 1e+05 | bin > 1e+06)
        stop("Bin size must be between 100kb and 1Mb")
    if (is.character(P) == F)
        stop("The P population must be passed as a charaster argument")
    if (is.character(P) == F)
        stop("The Q population must be passed as a charaster argument")

    list.of.packages <- c("vcfR", "gridExtra", "ggplot2")
    new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,
        "Package"])]
    if (length(new.packages))
        stop("vcfR is required")
    lapply(list.of.packages, require, character.only = T)


    if(is.character(file)) vcf <- read.vcfR(file, verbose=F)
    else if(is.character(file)==FALSE) vcf<-file
    a <- as.data.frame(vcf@gt)

    ### Get ref and alt allele counts for P population( i.e normal)
    if (is.null(a[[P]]))
        stop("The P population must be present in the provided vcf")
    striped <- a[[P]]
    Pref <- as.numeric(sapply(strsplit(striped, ":"), "[[", 4))
    Palt <- as.numeric(sapply(strsplit(striped, ":"), "[[", 6))

    ### Get ref and alt allele counts for Q population (i.e. exhibit feature)
    if (is.null(a[[Q]]))
        stop("The Q population must be present in the provided vcf")
    spotted <- a[[Q]]
    Qref <- as.numeric(sapply(strsplit(spotted, ":"), "[[", 4))
    Qalt <- as.numeric(sapply(strsplit(spotted, ":"), "[[", 6))



    #### MINOR ALLELE FREQUENCY FOR P POPULATION ######
    chrom <- as.data.frame(vcf@fix)$CHROM
    chrom <- as.numeric(ifelse(chrom == "X", 32, chrom))
    pos <- as.numeric(as.data.frame(vcf@fix)$POS)
    pos1 <- pos - 1
    P.minorAlleleFreq <- data.frame(chrom, pos1, pos, pmax(Pref, Palt), pmin(Pref,
        Palt))
    colnames(P.minorAlleleFreq)[c(4, 5)] <- c("PMaf", "Pmaf")
    P.minorAlleleFreq <- na.omit(P.minorAlleleFreq[with(P.minorAlleleFreq, order(chrom,
        pos1)), ])


    #### MINOR ALLELE FREQUENCY FOR Q POPULATION ######
    Q.minorAlleleFreq <- data.frame(chrom, pos1, pos, pmax(Qref, Qalt), pmin(Qref,
        Qalt))
    colnames(Q.minorAlleleFreq)[c(4, 5)] <- c("QMaf", "Qmaf")
    Q.minorAlleleFreq <- na.omit(Q.minorAlleleFreq[with(Q.minorAlleleFreq, order(chrom,
        pos1)), ])


    ################# GET HETEROZYGOSITY BINS for P POPULATION#################
    start = P.minorAlleleFreq$pos1[1]
    totalMajor = P.minorAlleleFreq$PMaf[1]
    totalMinor = P.minorAlleleFreq$Pmaf[1]
    j = 1
    k = 1
    df = matrix(ncol = 5, nrow = length(chrom))

    message(paste("Creating windows for ", P, "Heterozygosity", sep = ""), "\r", appendLF = T)
    flush.console()
    for (i in 2:nrow(P.minorAlleleFreq)) {
        if (P.minorAlleleFreq$pos[i] < (start + bin) & P.minorAlleleFreq$chrom[i] ==
            P.minorAlleleFreq$chrom[i - 1]) {
            totalMajor = totalMajor + as.numeric(P.minorAlleleFreq$PMaf[i])
            totalMinor = totalMinor + as.numeric(P.minorAlleleFreq$Pmaf[i])
        } else {
            row <- c(P.minorAlleleFreq$chrom[i - 1], P.minorAlleleFreq$pos[j], P.minorAlleleFreq$pos[j] +
                bin, totalMajor, totalMinor)
            j = i
            k = k + 1
            totalMajor = 0
            totalMinor = 0
            start = P.minorAlleleFreq$pos[i]
            df[k, ] <- row
        }
    }

    hetP <- as.data.frame(na.omit(df))
    names(hetP) <- c("chr", "start", "end", "major", "minor")
    hetP$Phet <- (hetP$major * hetP$minor * 2)/(hetP$major + hetP$minor)^2
    hetP$Pz <- (hetP$Phet - mean(hetP$Phet))/sd(hetP$Phet)


    ################## GET HETEROZYGOSITY BINS for Q POPULATION#################
    start = Q.minorAlleleFreq$pos1[1]
    totalMajor = Q.minorAlleleFreq$QMaf[1]
    totalMinor = Q.minorAlleleFreq$Qmaf[1]
    j = 1
    k = 1
    df = matrix(ncol = 5, nrow = length(chrom))

    message(paste("Creating windows for ", Q, "Heterozygosity", sep = ""), "\r", appendLF = T)
    flush.console()
    for (i in 2:nrow(Q.minorAlleleFreq)) {
        if (Q.minorAlleleFreq$pos[i] < (start + bin) & Q.minorAlleleFreq$chrom[i] ==
            Q.minorAlleleFreq$chrom[i - 1]) {
            totalMajor = totalMajor + as.numeric(Q.minorAlleleFreq$QMaf[i])
            totalMinor = totalMinor + as.numeric(Q.minorAlleleFreq$Qmaf[i])
        } else {
            row <- c(Q.minorAlleleFreq$chrom[i - 1], Q.minorAlleleFreq$pos[j], Q.minorAlleleFreq$pos[j] +
                bin, totalMajor, totalMinor)
            j = i
            k = k + 1
            totalMajor = 0
            totalMinor = 0
            start = Q.minorAlleleFreq$pos[i]
            df[k, ] <- row
        }
    }

    hetQ <- as.data.frame(na.omit(df))
    names(hetQ) <- c("chr", "start", "end", "major", "minor")
    hetQ$Qhet <- (hetQ$major * hetQ$minor * 2)/(hetQ$major + hetQ$minor)^2
    hetQ$Qz <- (hetQ$Qhet - mean(hetQ$Qhet))/sd(hetQ$Qhet)


    ################## Fixation Index###################
    Fst <- Q.minorAlleleFreq[, c(1, 2, 3)]
    Fst$QAC <- Qalt
    Fst$QRC <- Qref
    Fst$QT <- Qref + Qalt
    Fst$PRC <- Pref
    Fst$PAC <- Palt
    Fst$PT <- Pref + Palt

    Fst$QAF <- Fst$QAC/Fst$QT
    Fst$QRF <- Fst$QRC/Fst$QT
    Fst$PAF <- Fst$PAC/Fst$PT
    Fst$PRF <- Fst$PRC/Fst$PT

    n <- (Fst$QAF - Fst$PRF)^2 - ((Fst$QAF - Fst$QRF)/(Fst$QT - 1)) - ((Fst$PRF - Fst$PAF)/(Fst$PT -
        1))
    d <- Fst$QAF * Fst$PAF + Fst$QRF * Fst$PRF

    temp <- Fst[, c(1, 2, 3)]
    temp$n <- n
    temp$d <- d


    ##### BIN n & d########
    start = temp$pos1[1]
    N = temp$n[1]
    D = temp$d[1]
    j = 1  #index for start of bin
    k = 1  #index for adding row to new data frame
    df = matrix(ncol = 5, nrow = nrow(temp))

    message("Creating windows for Fst", "\r", appendLF = T)
    flush.console()
    for (i in 2:nrow(temp)) {
        if (temp$pos[i] < (start + bin) & temp$chrom[i] == temp$chrom[i - 1]) {
            N = N + as.numeric(temp$n[i])
            D = D + as.numeric(temp$d[i])
        } else {
            row <- c(temp$chrom[i - 1], temp$pos[j], temp$pos[j] + bin, N, D)
            j = i
            k = k + 1
            N = 0
            D = 0
            start = temp$pos[i]  #reset start of next bin to first row with value greater than previous start +bin
            df[k, ] <- row
        }
    }
    Fst <- as.data.frame(na.omit(df))
    colnames(Fst) <- c("chr", "start", "end", "N", "D")
    Fst <- Fst[which(Fst$N > -Inf & Fst$N < Inf), ]  #precautionary filter in case variants were called that were fixed in both pops
    Fst$Fst <- as.numeric(Fst$N)/as.numeric(Fst$D)
    Fst$FstZ <- (Fst$Fst - mean(Fst$Fst))/sd(Fst$Fst)
    Fst$FstP <- 2 * pnorm(-abs(Fst$FstZ))

    Fst <- Fst[, -c(4, 5)]
    hetQ <- hetQ[, c(1, 2, 3, 6, 7)]
    hetP <- hetP[, c(1, 2, 3, 6, 7)]

    data <- merge(Fst, hetQ, by = c("chr", "start", "end"))
    data <- merge(data, hetP, by = c("chr", "start", "end"))
    data$start <- as.numeric(data$start)
    data <- data[with(data, order(chr, start)), ]

    if (plot == TRUE) {

        ############# Distribution of heterozygosity and Fst #########
        h1 = hist(data$Phet)
        h2 = hist(data$Qhet)
        dist <- {
            par(mfrow = c(1, 2))
            hist(data$Phet, xlim = c(min(data$Phet, data$Qhet), max(data$Phet, data$Qhet)),
                ylim = c(0, max(h1$counts, h2$counts)), col = rgb(0, 0, 1, 0.25), main = NULL,
                xlab = "Heterozygosity", cex.lab = 1.5, cex.axis = 1)
            par(new = T)
            hist(data$Qhet, xlim = c(min(data$Phet, data$Qhet), max(data$Phet, data$Qhet)),
                ylim = c(0, max(h1$counts, h2$counts)), col = rgb(1, 0, 0, 0.25), main = NULL,
                xlab = "", xaxt = "n", ylab = "", yaxt = "n")
            legend("topright", legend = c(P, Q), col = c(rgb(0, 0, 1, 0.25), rgb(1,
                0, 0, 0.25)), pch = 15, 15)
            hist(data$Fst, xlim = c(min(data$Fst), max(data$Fst)), main = NULL, xlab = "Fst",
                cex.lab = 1.5, cex.axis = 1)
        }


        ############### GWAS ##################
        scaler <- as.data.frame(aggregate(end ~ chr, data, max))
        colnames(scaler)[2] <- "max"
        scaledData <- merge(data, scaler, by = "chr")


        scaledData$max <- as.numeric(scaledData$max)
        scaledData$end <- as.numeric(scaledData$end)
        scaledData$ScalePos <- scaledData$chr + (scaledData$end/scaledData$max)
        scaledData$hetDiff <- scaledData$Qhet - scaledData$Phet
        scaledData$zDiff <- (scaledData$hetDiff - mean(scaledData$hetDiff)/sd(scaledData$hetDiff))
        scaledData$pDiff <- 2 * pnorm(-abs(scaledData$zDiff))


        p1 <- ggplot(scaledData, aes(ScalePos, -log10(FstP), col = as.factor(chr))) +
            geom_point() + theme_bw() + theme(legend.position = "none") + xlab("") +
            geom_hline(yintercept = qnorm(0.95), color = "red", )
        p2 <- ggplot(scaledData, aes(ScalePos, -log10(pDiff), col = as.factor(chr))) +
            geom_point() + theme_bw() + theme(legend.position = "none") + xlab("Chomosome") +
            ylab("-log10(p) Difference in Het")
        gwas <- grid.arrange(p1, p2, ncol = 1)

        return(list(data = data, distribution = dist, gwas = gwas))
    } else return(data)

}
