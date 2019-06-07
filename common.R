require(png)

plusImage <- readPNG("www/plus.png")
minusImage <- readPNG("www/minus.png")

# write either warning or message to console with a time stamp
writeLog <- function(..., warn=FALSE) {
    level.func <- ifelse(warn, warning, message)
    level.func(Sys.time(),": ", ...)
}

formatInterval <- function(interval) {
    sprintf("%s:%s-%s", interval$CHROM, format(interval$START, nsmall=0, big.mark=",", scientific=FALSE), format(interval$END, nsmall=0, big.mark=",", scientific=FALSE))
}

formatLocusInterval <- function(interval) {
    intervalLength = interval$END - interval$START + 1
    sprintf("%s -- %s bp", formatInterval(interval), format(intervalLength, nsmall=0, big.mark=",", scientific=FALSE))
}

formatCallInterval <- function(chrom, start, end) {
    sprintf("%s:%s-%s", chrom, format(start, nsmall=0, big.mark=","), format(end, nsmall=0, big.mark=","))
}

readTabixedData <- function(dataFile, window) {
    # Note that the header start with #, so it needs to be stripped
    cmd = sprintf("zcat -f %s | head -1 | sed -e '1s/^\\#//' && tabix %s %s:%.0f-%0.f", dataFile, dataFile, window$CHROM, window$START, window$END)
    read.table(pipe(cmd), header=T, sep="\t", stringsAsFactors=F)
}

getPlotSection <- function(y, plotLayout) {
    if (y < 0) {
        return(NULL)
    }
    plotLayout[plotLayout$Y0 == max(plotLayout[y >= plotLayout$Y0, "Y0"]), ]
}

getPlotSectionHeight <- function(df) {
    sectionHeight = 0
    if (!is.null(df) && nrow(df) > 0) {
        sectionHeight = max(df$OFFSET)
    }
    return(sectionHeight)
}

getPlotOffsetByCoordinate <- function(y, plotSection) {
    offset = round(plotSection$SCALE * (y - plotSection$Y0) + .5)
    if (plotSection$DIRECTION == -1) {
        offset = plotSection$HEIGHT - offset + 1
    }
    offset
}

getOverlappingItem <- function(x, offset, df) {
    if (is.null(df)) {
        return(NULL)
    }
    df[df$OFFSET == offset & df$START <= x & df$END >= x, ]
}

computeSegmentOffsets <- function(segmentStarts, segmentEnds, segmentStackingGap) {
    if (is.null(segmentStarts) || is.null(segmentEnds)) {
        return(NULL)
    }
    if (length(segmentEnds) != length(segmentStarts)) {
        writeLog("ERROR: length(segmentEnds) != length(segmentStarts)") 
        return(NULL)
    }

    # The n-th eleement contains the end position for the last call with OFFSET=n
    segmentEndList = integer()
    offsets = integer(length(segmentStarts))
    for (idx in 1:length(segmentStarts)) {
        for (offset in 1:(length(segmentEndList)+1)) {
            endPos = segmentEndList[offset]
            if (is.na(endPos) || (segmentStarts[idx] - endPos) >= segmentStackingGap) {
                offsets[idx] = offset
                segmentEndList[offset] = segmentEnds[idx]
                break
            }
        }
    }
    return(offsets)
}

plotOrientationArrows <- function(df, window, gridXcoords, xScale, yCoords, arrowColor) {
    arrowLength = (window[2] - window[1]) / xScale / 100
    for (idx in 1:nrow(df)) {
        arrowX = gridXcoords[gridXcoords - arrowLength > max(window[1], df[idx, "START"])/xScale & gridXcoords + arrowLength < df[idx, "END"]/xScale]
        if (length(arrowX) > 0) {
            direction = ifelse(df[idx, "STRAND"] == "+", 1, -1)
            arrows(x0=arrowX, x1=arrowX + direction*arrowLength, y0=yCoords[idx], length=.05, col=arrowColor)
        }
    }
}

