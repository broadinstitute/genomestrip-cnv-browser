AUTO_MERGE_TARGET = 100

cmdDefaults = list(PLOT_DEPTH_LIMIT = 10, XSCALETEXT = "(Kb)", MERGEBINS = "auto", MAXGAPLENGTH = 1000)

cmdPlotGroups = NULL

cmdMinimumPadding = 1000
cmdMinimumPaddingFraction = 0.5

cmdMergeBins = NULL
cmdMergeOffset = NULL
cmdPlotGroups = NULL
cmdBackgroundSampleList = NULL

endsWith <- function(string, suffix) {
    n = nchar(string)
    m = nchar(suffix)
    return (n >= m && substr(string, n-m+1, n) == suffix)
}

getCopyNumber <- function(profileData, regionStart, regionEnd, variantSamples) {
    # Use only completely contained bins
    binMask = (profileData$BINSTARTS >= regionStart) & (profileData$BINENDS <= regionEnd)
    numerator = sum(profileData$COVNUMERATORS[binMask,variantSamples], na.rm=T)
    denominator = sum(profileData$COVDENOMINATORS[binMask,variantSamples], na.rm=T)
    return(ifelse(denominator == 0, NA, numerator/denominator))
}

computeMergeIndex <- function(profileData, binsToMerge, maxGapLength, mergeOffset) {
    binStarts = profileData$BINSTARTS
    binEnds = profileData$BINENDS
    binGaps = profileData$BINGAPS
    nbins = length(binStarts)
    binIndex = c()
    if (mergeOffset > 0) {
        binIndex = rep(NA, mergeOffset)
    }
    idxValue = 1
    idxCount = 0
    firstBin = 1 + mergeOffset
    if (firstBin <= nbins) {
        for (i in firstBin:nbins) {
            binIndex = c(binIndex, idxValue)
            idxCount = idxCount + 1
            if (idxCount == binsToMerge || binGaps[i] > maxGapLength) {
                idxValue = idxValue + 1
                idxCount = 0
            }
        }
    }
    return(binIndex)
}

mergeProfileBins <- function(profileData, binsToMerge, mergeOffset) {
    # for now, we drop the GCRATIOs (these should really be parsed in the original reading code)
    chrom = profileData$CHROM
    binStarts = profileData$BINSTARTS
    binEnds = profileData$BINENDS
    covnums = profileData$COVNUMERATORS
    covdenoms = profileData$COVDENOMINATORS
    maxGapLength = cmdDefaults$MAXGAPLENGTH
    binIndex = computeMergeIndex(profileData, binsToMerge, maxGapLength, mergeOffset)
    outBins = max(binIndex)
    binStarts = tapply(binStarts, binIndex, min)
    binEnds = tapply(binEnds, binIndex, max)
    binGaps = c(tail(binStarts, -1) - head(binEnds, -1) - 1, 0)
    binGaps = ifelse(binGaps < 0, 0, binGaps)
    covnums = apply(covnums, 2, function(col) { tapply(col, binIndex, sum) })
    if (is.null(dim(covnums))) {
        covnums = matrix(covnums, nrow=1)
    }
    covdenoms = apply(covdenoms, 2, function(col) { tapply(col, binIndex, sum) })
    if (is.null(dim(covdenoms))) {
        covdenoms = matrix(covdenoms, nrow=1)
    }

    profileData = list(CHROM=chrom,
                       MERGEDBINS=binsToMerge,
                       BINSTARTS=binStarts,
                       BINENDS=binEnds,
                       BINGAPS=binGaps,
                       COVNUMERATORS=covnums,
                       COVDENOMINATORS=covdenoms)
    return(profileData)
}

runCommand <- function(profileFile, window, field) {
    cmd = ifelse(endsWith(profileFile, ".gz"), "zcat", "cat")
    cmd = sprintf("%s %s | head -1 && tabix  %s %s:%.0f-%0.f | awk -v field=%d -f extract_profile_data.awk", cmd, profileFile, profileFile, window$CHROM, window$START, window$END, field)
    read.table(pipe(cmd), header=T, stringsAsFactors=F, row.names=1, check.names=F, sep="\t")
}

readProfileDataForRegion <- function(profileFile, window) {
    counts = runCommand(profileFile, window, 1)
    expected = runCommand(profileFile, window, 2)

    binStarts = counts$START
    binEnds = counts$END
    binGaps = c(tail(binStarts, -1) - head(binEnds, -1) - 1, 0)
    binGaps = ifelse(binGaps < 0, 0, binGaps)
    counts = counts[, 6:ncol(counts)]
    expected = expected[, 6:ncol(expected)]
    profileData = list(CHROM=window$CHROM,
                       MERGEDBINS=1,
                       BINSTARTS=binStarts,
                       BINENDS=binEnds,
                       BINGAPS=binGaps,
                       COVNUMERATORS=counts,
                       COVDENOMINATORS=expected)

    return(profileData)
}

formatLengthText <- function(length) {
    if (length > 1000000) {
        return(sprintf("%1.1f Mb", length/1000000))
    } else if (length > 1000) {
        return(sprintf("%1.1f Kb", length/1000))
    } else {
        return(sprintf("%d bp", round(length)))
    }
}

getBinsToMerge <- function(numBins, binsToMerge) {
    binsToMerge = ifelse(is.null(binsToMerge), 1, binsToMerge)
    if (binsToMerge == "auto") {
        binsToMerge = max(1, floor(numBins / AUTO_MERGE_TARGET))
    }
    binsToMerge
}

plotAuxDelData <- function(delData, window, plotSection, dataColor) {
    if (is.null(delData)) {
        return()
    }

    df = delData[delData$RIGHTEND >= window[1] & delData$LEFTSTART <= window[2], ]
    if (nrow(df) == 0) {
        return()
    }

    offsets = df$OFFSET - .5
    segments(
        x0 = df$LEFTSTART / PROFILE_PLOT_XSCALE,
        x1 = df$RIGHTEND / PROFILE_PLOT_XSCALE,
        y0 = plotSection$Y0 + offsets / plotSection$SCALE,
        col = dataColor
    )
    rect(
        xleft = c(df$LEFTSTART, df$RIGHTSTART) / PROFILE_PLOT_XSCALE,
        xright = c(df$LEFTEND, df$RIGHTEND) / PROFILE_PLOT_XSCALE,
        ytop = plotSection$Y0 + (offsets + .1) / plotSection$SCALE,
        ybottom = plotSection$Y0 + (offsets - .1) / plotSection$SCALE,
        col = dataColor
    )
}

plotProfileData <- function(profileData,gtData,  maxCopyNumber, plotMedian, variantSamples, selectedSamples, selectedSampleColor, plotSection) {
    writeLog("plotProfileData")
    y0 = plotSection$Y0
    yScale = plotSection$SCALE
    coverage = profileData$COVERAGE

    selectedSampleColor = ifelse(is.null(selectedSampleColor), "red", selectedSampleColor)

    # lines to help with visualization
    abline(h=(y0 + 1:floor(maxCopyNumber) / yScale), col="black", lwd=1, lty="dashed")

    # the actual depth data
    xvals = as.vector(mapply(c,profileData$BINSTARTS,profileData$BINENDS))
    xvals = xvals / PROFILE_PLOT_XSCALE
    coverage[coverage > maxCopyNumber] = maxCopyNumber
    apply(coverage, 2, function(y) {
        lines(xvals, rep(y0 + y / yScale, each=2), col="lightgrey")
    })

    if (plotMedian) {
        yvals = y0 + rep(apply(coverage,1,median), each=2) / yScale
        lines(xvals, yvals, col="black")
    }

    # plot the target and selected sample
    for (sample in c(variantSamples, selectedSamples)) {
        if (sample %in% colnames(coverage)) {
            sampleColor = ifelse(sample %in% selectedSamples, selectedSampleColor, colors[paste0("CN", gtData[gtData$SAMPLE == sample, "CN"])])
            if (is.na(sampleColor)) {
                sampleColor = "blue"
            }
            yvals = y0 + rep(coverage[, sample], each=2) / yScale
            lines(xvals, yvals, col=sampleColor, lwd=1.5)
        }
    }
}

makeProfilePlot <- function(searchResults, plotDims, xRange, maxCopyNumber, plotMedian, variantSamples, selectedSamples, selectedSampleColor, gtData, plotLayout, plotHeight) {
    writeLog("makeProfilePlot")

    par(mar=c(3, 2, .2, 0),xpd=F)

    xlim = xRange / PROFILE_PLOT_XSCALE
    plot(x=NULL, y=NULL,
         xlim=xlim,
         ylim=c(0, 1),
         axes=F,
         xaxs="i",
         yaxs="i",
         xlab="",
         ylab=""
    )

    axis(1, at=axTicks(1), format(axTicks(1), nsmall=0, big.mark=",", scientific=FALSE))

    profilePlotLayout = plotLayout["Profile", ]
    axis(2, labels=seq(1, maxCopyNumber, 2), at=(profilePlotLayout$Y0 + seq(1, maxCopyNumber, 2)/profilePlotLayout$SCALE))

    mtext(sprintf("Chromosome %s %s", searchResults$profileData$CHROM, cmdDefaults$XSCALETEXT), side=1, line=2, at=mean(xlim))

    abline(h=plotLayout[, "Y0"], col="blue")

    gridXcoords = seq(xlim[1], xlim[2], diff(xlim)/80)
    plotGeneData(searchResults, xRange, gridXcoords, PROFILE_PLOT_XSCALE, plotLayout["Gene", ], plotHeight)

    plotProfileData(searchResults$profileData, gtData, maxCopyNumber, plotMedian, variantSamples, selectedSamples, selectedSampleColor, plotLayout["Profile", ])

    plotAuxDelData(searchResults$discoveryPairData, xRange, plotLayout["Discovery", ], "navyblue")
    plotAuxDelData(searchResults$genotypesPairData, xRange, plotLayout["Genotypes", ], "navyblue")

    for (sectionIdx in 1:nrow(plotLayout)) {
        layoutSection = plotLayout[sectionIdx, ]
        mtext(layoutSection$SECTION_TITLE, side=2, at=layoutSection$Y0 + .5 * layoutSection$SCALED_HEIGHT, padj=0, xpd=T, cex=.8)
    }

    for (sectionName in c("Discovery", "Genotypes")) {
        plotSectionControls(plotDims, xRange, plotLayout[sectionName, ])
    }

    for (call in searchResults$selectedCall[, c("START", "END")]) {
        abline(v=call/PROFILE_PLOT_XSCALE, lty="dashed", col="red")
    }

    box()

    writeLog("Done makeProfilePlot")
}

getProfilePlotSamples <- function(x, y, searchResults) {
    yBandHeight = 1 / 40

    profileData = searchResults$profileData
    coverage = profileData$COVERAGE[profileData$BINSTARTS <= x & profileData$BINENDS >= x, , drop=FALSE]

    samples = unlist(sapply(colnames(coverage), function(sample) {
        sampleCoverage = coverage[, sample]
        if (any(sampleCoverage - yBandHeight <= y & sampleCoverage + yBandHeight >= y)) {
            return(sample)
        }
    }))
    samples = samples[!is.null(samples)]
    if (length(samples) > 0) {
        return(samples)
    }

    xBandWidth = (searchResults$window$END - searchResults$window$START + 1) / 500
    binIndices = which(profileData$BINSTARTS - xBandWidth <= x & profileData$BINSTARTS + xBandWidth >= x | profileData$BINENDS - xBandWidth <= x & profileData$BINENDS + xBandWidth >= x)

    if (length(binIndices) != 2) {
         return(c())
    }

    cov1 = unlist(profileData$COVERAGE[binIndices[1], ])
    cov2 = unlist(profileData$COVERAGE[binIndices[2], ])
    isOpposite = (cov1 < y & cov2 > y | cov1 > y & cov2 < y)
    names(isOpposite[isOpposite])
}

getProfileTooltip <- function(x, y, plotSection, searchResults) {
    y = (y - plotSection$Y0) * plotSection$SCALE
    samples = getProfilePlotSamples(x, y, searchResults)
    paste0(samples, collapse="; ")
}

getAuxDataItem <- function(x, offset, df) {
    if (is.null(df)) {
        return(NULL)
    }
    dataItem = df[df$OFFSET == offset & df$LEFTSTART <= x & df$RIGHTEND >= x, ]
    if (nrow(dataItem) == 0) {
        return(NULL)
    }
    dataItem
}

getAuxItemSpan <- function(auxItem) {
    c(
        formatCallInterval(auxItem$LEFTCHR, auxItem$LEFTSTART, auxItem$LEFTEND),
        formatCallInterval(auxItem$RIGHTCHR, auxItem$RIGHTSTART, auxItem$RIGHTEND))
}

getDiscoveryPairTooltip <- function(x, y, plotSection, searchResults) {
    offset = getPlotOffsetByCoordinate(y, plotSection)
    auxItem = getAuxDataItem(x, offset, searchResults$discoveryPairData)
    if (is.null(auxItem)) {
        return(NULL)
    }

    span = auxItem$RIGHTSTART - auxItem$LEFTEND - 1
    tooltipItems = c(
        auxItem$READNAME,
        formatCallInterval(auxItem$CHR, auxItem$LEFTSTART, auxItem$LEFTEND),
        sprintf("%s, span=%d", formatCallInterval(auxItem$CHR, auxItem$RIGHTSTART, auxItem$RIGHTEND), span),
        sprintf("%s %s", auxItem$SAMPLE, auxItem$TYPE)
    )
}

getGenotypesPairTooltip <- function(x, y, plotSection, searchResults) {
    offset = getPlotOffsetByCoordinate(y, plotSection)
    auxItem = getAuxDataItem(x, offset, searchResults$genotypesPairData)
    if (is.null(auxItem)) {
        return(NULL)
    }
    span = auxItem$RIGHTSTART - auxItem$LEFTEND - 1
    tooltipItems = c(
        auxItem$READNAME,
        formatCallInterval(auxItem$LEFTCHR, auxItem$LEFTSTART, auxItem$LEFTEND),
        sprintf("%s, span=%d", formatCallInterval(auxItem$RIGHTCHR, auxItem$RIGHTSTART, auxItem$RIGHTEND), span),
        sprintf("%s %s", auxItem$SAMPLE, auxItem$PAIRTYPE)
    )

    if (!is.na(searchResults$selectedCall$GSBKPT)) {
        tooltipItems = c(tooltipItems, searchResults$selectedCall$GSBKPT)
    }
    tooltipItems
}

getSectionControlDims <- function(plotDims, xRange, layoutSection) {
    controlXCoord = xRange[1] / PROFILE_PLOT_XSCALE
    controlYCoord = layoutSection$Y0 + layoutSection$SCALED_HEIGHT - 0.2 / plotDims[2]
    controlWidth = .15 / plotDims[1] * (xRange[2] - xRange[1]) / PROFILE_PLOT_XSCALE
    controlHeight = .15 / plotDims[2]

    list(x=controlXCoord, y=controlYCoord, width=controlWidth, height=controlHeight, gap=.2)
}

plotSectionControls <- function(plotDims, xRange, layoutSection) {
    controlDims = getSectionControlDims(plotDims, xRange, layoutSection)
    rasterImage(minusImage, controlDims$x, controlDims$y, controlDims$x + controlDims$width, controlDims$y + controlDims$height)
    rasterImage(plusImage, controlDims$x + (1 + controlDims$gap) * controlDims$width, controlDims$y, controlDims$x + (2 + controlDims$gap) * controlDims$width, controlDims$y + controlDims$height)
}

processAuxDataClick <- function(x, y, plotDims, xRange, layoutSection) {
    controlDims = getSectionControlDims(plotDims, xRange, layoutSection)

    if (y < controlDims$y || y > controlDims$y + controlDims$height) {
        return(NULL)
    }

    button = NULL
    if (x >= controlDims$x && x <= controlDims$x + controlDims$width) {
        button = "minus"
    } else if (x >= controlDims$x + (1 + controlDims$gap) * controlDims$width && x <= controlDims$x + (2 + controlDims$gap) * controlDims$width) {
        button = "plus"
    }
    return(button)
}

