XSCALE <- 1000

createLocusPlotLayout <- function(searchResult) {
    # Calls and genes have the height of 1, segdups are .25
    geneHeight = getPlotSectionHeight(searchResult$geneTxDF)
    segdupHeight = getPlotSectionHeight(searchResult$segdupDF)
    callHeight = getPlotSectionHeight(searchResult$callsDF)

    # First allocate space between segdups and calls+genes, then between calls and genes
    callHeightRange = c(.3, .7)
    geneHeightRange = c(.3, .7)
    segdupHeightRange = c(.1, .2)

    segdupFraction = .2
    if (segdupHeight != 0) {
        segdupFraction = .25*segdupHeight / (callHeight + geneHeight + .25*segdupHeight)
        segdupFraction = min(max(segdupFraction, segdupHeightRange[1]), segdupHeightRange[2])
    }

    callFraction = .5
    geneFraction = .5
    if (callHeight != 0 || geneHeight != 0) {
        callFraction = callHeight / (callHeight + geneHeight)
        geneFraction = geneHeight / (callHeight + geneHeight)
    }
    callFraction = (1 - segdupFraction) * min(max(callFraction, callHeightRange[1]), callHeightRange[2])
    geneFraction = (1 - segdupFraction) * min(max(geneFraction, geneHeightRange[1]), geneHeightRange[2])

    plotLayout = data.frame(
        SECTION_NAME = c("GENE", "SEGDUP", "CALL"),
        SCALED_HEIGHT = c(geneFraction, segdupFraction, callFraction),
        HEIGHT = c(geneHeight, segdupHeight, callHeight),
        DIRECTION = c(-1, 1, 1),
        TOOLTIP_FUNCTION = c("getGeneTooltip", "getSegdupTooltip", "getCallTooltip"),
        stringsAsFactors=F
    )
    rownames(plotLayout) = plotLayout$SECTION_NAME

    plotLayout$Y0 <- cumsum(c(0, head(plotLayout[, "SCALED_HEIGHT"], -1)))
    plotLayout$SCALE <- plotLayout$HEIGHT / plotLayout$SCALED_HEIGHT

    plotLayout
}

getCallColors <- function(categories) {
    ifelse(categories == "DEL", DEL_COLOR,
        ifelse(categories == "DUP", DUP_COLOR, MIXED_COLOR))
}

getCallHeights <- function(calls, useSmoothFreqHeight=TRUE) {
    if (useSmoothFreqHeight) {
        (log10(calls$NDEL+calls$NDUP) + 1) / (1+log10(2*calls$NCALL)) / 3.5
    } else {
        ifelse(calls$NDEL + calls$NDUP == 1,
            0.075,
            ifelse(calls$FREQ <= 0.1,
                0.125,
                0.225)
        )
    }
}

plotCallData <- function(callsDF, selectedCallID, window, plotSection) {
    if (is.null(callsDF)) {
        return()
    }

    callColors = getCallColors(callsDF$CATEGORY)
    callHeights = getCallHeights(callsDF)

    y0 = plotSection$Y0
    yScale = plotSection$SCALE
    rect(
        xleft = pmax(window$START, callsDF$START) / XSCALE,
        xright = (callsDF$END+1) / XSCALE,
        ybottom = y0 + (callsDF$OFFSET - .5 - callHeights) / yScale,
        ytop = y0 + (callsDF$OFFSET - .5 + callHeights) / yScale,
        col = callColors,
        border = ifelse(callsDF$ID == selectedCallID, SELECTED_COLOR, NA)
    )

    if (!is.null(blacklistedRegions)) {
        bRegions = blacklistedRegions[blacklistedRegions$CHROM == window$CHROM & blacklistedRegions$END >= window$START & blacklistedRegions$START <= window$END, ]
        if (nrow(bRegions) > 0) {
            rect(
                xleft = pmax(window$START, bRegions$START) / XSCALE,
                xright = pmin(window$END, bRegions$END) / XSCALE,
                ybottom = y0,
                ytop = y0 + plotSection$SCALED_HEIGHT,
                col = "lightgrey",
                border = NA
            )
        }
    }
}

plotSegdupData <- function(segdupDF, window, gridXcoords, plotSection) {
    if (is.null(segdupDF)) {
        return()
    }

    y0 = plotSection$Y0
    yScale = plotSection$SCALE
    rect(
        xleft = pmax(window[1], segdupDF$START) / XSCALE,
        xright = (segdupDF$END+1) / XSCALE,
        ybottom = y0 + (segdupDF$OFFSET - 0.6) / yScale,
        ytop = y0 + (segdupDF$OFFSET - 0.4) / yScale,
        col = "darkorange3",
        border = NA
    )

    if (gridXcoords[2] - gridXcoords[1] <= 10) {
        plotOrientationArrows(segdupDF, window, gridXcoords, XSCALE, y0 + (segdupDF$OFFSET - .5) / yScale, "white")
    }
}

makeLocusPlot <- function(searchResults, plotLayout, plotHeight, selectedCallID=NA) {
    writeLog(sprintf("makeLocusPlot, interval: %s:%s-%s", searchResults$window$CHROM, searchResults$window$START, searchResults$window$END))

    window = c(searchResults$window$START, searchResults$window$END)

    par(mar=c(3,1,0,0),xpd=T)

    plot(x=NULL, y=NULL,
         xlim=window/XSCALE,
         ylim=c(0, 1),
         axes=F,
         xaxs="i",
         yaxs="i",
         xlab="",
         ylab=""
    )

    # Make the gridlines
    gridXcoords = seq(window[1], window[2], diff(window)/80) / XSCALE
    segments(x0=gridXcoords, y0=0, y1=1, col = "#00609f", lwd=0.1)

    # Plot track names
    textStart = (81*window[1] - window[2]) / 80 / XSCALE
    text(textStart, (1+plotLayout["CALL", "Y0"])/2-.1, "CNV calls", pos=4, srt=90)
    text(textStart, (plotLayout["CALL", "Y0"] + plotLayout["SEGDUP", "Y0"])/2-.05, "Seg dups", pos=4, srt=90)
    text(textStart, plotLayout["SEGDUP", "Y0"]/2-.02, "Genes", pos=4, srt=90)

    # Plot track separation lines
    segments(x0=window[1]/XSCALE, y0=c(plotLayout["CALL", "Y0"], plotLayout["SEGDUP", "Y0"]), x1=window[2]/XSCALE, y1=c(plotLayout["CALL", "Y0"], plotLayout["SEGDUP", "Y0"]), col="#7e725b")

    plotGeneData(searchResults, window, gridXcoords, XSCALE, plotLayout["GENE", ], plotHeight)
    plotCallData(searchResults$callsDF, selectedCallID, searchResults$window, plotLayout["CALL", ])
    plotSegdupData(searchResults$segdupDF, window, gridXcoords, plotLayout["SEGDUP", ])

    axis(1, at=axTicks(1), format(axTicks(1), nsmall=0, big.mark=",", scientific=FALSE))
    title(xlab=paste("HG38 Chromosome", searchResults$window$CHROM,"(kb)"),line=2)
    title(line=0.2)
    box()
}

makeLosusLegend <- function() {
    par(mar=c(0,0,0,0))

    plot(0, type='n', axes=FALSE, ann=FALSE, xlim=c(0, 10), ylim=c(0, 10))

    xLeft = 0
    yCoords = c(5, 7, 3)
    rect(xLeft, yCoords, xLeft+1, yCoords+1, col=c(DUP_COLOR, DEL_COLOR, "white"),border = c(NA,NA,SELECTED_COLOR))
    text(x=2, y=9, "Type")
    text(x=xLeft+2, y=yCoords+0.5, labels=c("Duplication", "Deletion", "Selected"), adj=0)
}

getLocusPlotMouseItem <- function(searchResults, plotLayout, mouseEvent) {
    x = mouseEvent$x * 1000
    y = mouseEvent$y

    plotSection = plotLayout[plotLayout$Y0 == max(plotLayout[y >= plotLayout$Y0, "Y0"]), ]
    if (nrow(plotSection) == 0) {
        return(NULL)
    }
    offset = getPlotOffsetByCoordinate(y, plotSection)

    type = plotSection$SECTION_NAME
    name = NULL
    exonDescr = NULL
    if (type == "CALL") {
        call = getOverlappingItem(x, offset, searchResults$callsDF)
        if (!is.null(call) && nrow(call) == 1) {
            name = call$ID
        }
    } else if (type == "GENE") {
        tx = getOverlappingItem(x, offset, searchResults$geneTxDF)
        if (!is.null(tx) && nrow(tx) == 1) {
            name = tx$GENE
            geneExonDF = searchResults$geneExonDF
            exon = geneExonDF[geneExonDF$TX_ID == rownames(tx) & geneExonDF$START <= x & geneExonDF$END >= x, ]
            if (!is.null(exon) && nrow(exon) == 1) {
                exonDescr = sprintf("%d/%d", exon$ORDINAL, tx$EXON_COUNT)
            } 
        }
    }

    mouseItem = NULL
    if (!is.null(name)) {
        mouseItem = list(type=type, name=name, exonDescr=exonDescr)
    }
    return(mouseItem)
}

getCallTooltip <- function(x, y, plotSection, searchResults) {
    offset = getPlotOffsetByCoordinate(y, plotSection)
    call = getOverlappingItem(x, offset, searchResults$callsDF)
    if (!is.null(call) && nrow(call) == 1) {
        sprintf("Freq: %.3f%%, %s (%s)", 100*call$FREQ, formatInterval(call), call$ID)
    }
}

getSegdupTooltip <- function(x, y, plotSection, searchResults) {
    offset = getPlotOffsetByCoordinate(y, plotSection)
    segdup = getOverlappingItem(x, offset, searchResults$segdupDF)
    if (!is.null(segdup) && nrow(segdup) == 1) {
        c(sprintf("%s", formatInterval(segdup)),
          sprintf("%s", formatCallInterval(segdup$CHROM.2, segdup$START.2, segdup$END.2)),
          sprintf("Match: %.2f%%", 100*segdup$FRAC_MATCH))
    }
}
