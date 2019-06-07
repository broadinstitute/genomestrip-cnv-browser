plotGeneData <- function(geneData, window, gridXcoords, xScale, plotSection, plotHeight) {
    sectionHeight = plotSection$SCALED_HEIGHT * plotHeight

    geneTxDF = geneData$geneTxDF
    geneExonDF = geneData$geneExonDF

    if (is.null(geneTxDF)) {
        return()
    }

    geneExonDF = geneExonDF[geneExonDF$END >= window[1] & geneExonDF$START <= window[2], ]

    y0 = plotSection$Y0
    yScale = plotSection$SCALE
    numRows = max(geneTxDF$OFFSET)
    txYCoords = y0 + (numRows - geneTxDF$OFFSET + .5) / yScale
 
    segments(
        x0 = pmax(window[1], geneTxDF$START) / xScale,
        x1 = geneTxDF$END / xScale,
        y0 = txYCoords,
        col = "#00609f"
    )

    if (gridXcoords[2] - gridXcoords[1] <= 10) {
        plotOrientationArrows(geneTxDF, window, gridXcoords, xScale, txYCoords, "#00609f")
    }

    if (!is.null(geneExonDF)) {
        # .4 is the fraction of the gene height from the full row height
        # Gowever, don't let it be more than .2 in in real height
        heightFactor = 1 / max(1, 2 * sectionHeight / numRows)
        halfHeight = heightFactor * ifelse(geneExonDF$IS_UTR, 0.1, 0.2) / yScale
        exonOffsets = geneTxDF[geneExonDF$TX_ID, "OFFSET"]
        exonYCoords = y0 + (numRows - exonOffsets + .5) / yScale
        rect(
            xleft = pmax(window[1], geneExonDF$START) / xScale,
            xright = (geneExonDF$END+1) / xScale,
            ytop = exonYCoords + halfHeight,
            ybottom = exonYCoords - halfHeight,
            col = "#00609f",
            border="#00609f"
        )
    }

    textCex = 2.5 * sectionHeight / numRows
    textCex = min(textCex, .9)
    if (textCex >= .4) {
        text(
            x = (pmax(window[1], geneTxDF$START) + pmin(window[2], geneTxDF$END+1)) / 2 / xScale,
            y = txYCoords + 0.3 / yScale,
            labels = geneTxDF$GENE,
            cex=textCex
        )
    }
}

getGeneTooltip <- function(x, y, plotSection, searchResults) {
    tooltip = NULL

    geneTxDF = searchResults$geneTxDF
    geneExonDF = searchResults$geneExonDF

    offset = getPlotOffsetByCoordinate(y, plotSection)
    tx = geneTxDF[geneTxDF$OFFSET == offset & geneTxDF$START <= x & geneTxDF$END >= x, ]
    if (!is.null(tx) && nrow(tx) == 1) {
        gene = tx$GENE
        exon = tail(geneExonDF[geneExonDF$TX_ID == rownames(tx) & geneExonDF$START <= x, ], 1)
        if (nrow(exon) == 1) {
            if (exon$END >= x) {
                tooltip = sprintf("%s; exon %d/%d, %s", gene, exon$ORDINAL, tx$EXON_COUNT, formatInterval(exon))
            } else {
                intronOrdinal = ifelse(tx$STRAND == "+", exon$ORDINAL, exon$ORDINAL - 1)
                tooltip = sprintf("%s; intron %d/%d", gene, intronOrdinal, tx$EXON_COUNT-1)
            }
        }
    }
    return(tooltip)
}

