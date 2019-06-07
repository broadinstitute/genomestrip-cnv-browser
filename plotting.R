plotGenotypesHistogram <- function(data, baseHist, greyNonConfidentCallss, yAxisLimit) {
    if (is.null(data)) {
        return(NULL)
    }

    copyNumbers = sort(unique(data$CN))
    if (greyNonConfidentCallss) {
        data[data$CNQ < 13, "CN"] <- NA
        copyNumbers = c(copyNumbers, NA)
    }

    copyNumberColors = colors[ifelse(is.na(copyNumbers), "NC", paste0("CN", copyNumbers))]

    yAxisLimit = suppressWarnings(as.numeric(yAxisLimit))
    if (length(yAxisLimit) == 0 || is.na(yAxisLimit)) {
        yLim = NULL
    } else {
        yLim = c(0, yAxisLimit)
    }

    #Makes the histogram with enough room for legend
    par(mar=c(2.5,3.5,0.4,4.6), xpd=T)
    
    plot(baseHist, ylim=yLim, xlab="", ylab="", main="", axes=F, family="sans", font.main=1, border=NA)

    title(xlab="Copy number",line=1.5)
    title(ylab="Count",line=2.5)

    axis(1,las=1, at=seq(0,round(max(data[,"CNF"]))+1,1), pos=c(0))
    axis(2,las=1)

    legend(x=(max(baseHist$breaks)-.5), y=par("yaxp")[2],
        legend=ifelse(names(copyNumberColors) == "CN-1", "* allele", names(copyNumberColors)),
        col=copyNumberColors, pch=15, bty="n",
        title="Copy number\ncalled as:", title.adj=0.5,
        y.intersp=1.25, pt.cex=1.25)

    NUM_BINS = length(baseHist$breaks) - 1
    baseHeights = rep(0, NUM_BINS)
    binLeft = baseHist$breaks[1:NUM_BINS]
    binRight = baseHist$breaks[2:(NUM_BINS+1)]
    binWidth = binRight[1] - binLeft[1]
    for (idx in 1:length(copyNumbers)) {
        cn = copyNumbers[idx]
        if (is.na(cn)) {
            selected = is.na(data$CN)
        } else {
            selected = !is.na(data$CN) & data$CN == cn
        }
        cnData = data[selected, "CNF"]
        cnCounts = tapply(cnData, pmax(1, ceiling(cnData/binWidth)), length)
        binData = rep(NA, NUM_BINS)
        names(binData) = seq(1, NUM_BINS)
        binData[names(cnCounts)] = cnCounts

        rect(
            xleft=binLeft,
            xright=binRight,
            ybottom=baseHeights,
            ytop=baseHeights+binData,
            col=copyNumberColors[idx],
            border=NA
        )
        baseHeights = baseHeights + ifelse(is.na(binData), 0, binData)
    }
}

plotPopulationHistogram <- function(PopCounts, plotCounts) {
    if (plotCounts) {
        df = PopCounts
        yLabel = "Count"
    } else {
        df = PopCounts / rowSums(PopCounts)
        yLabel = "Population Frequency"
    }

    #Makes the histogram with enough room for legend
    par(mar=c(2.5,3.5,0.4,8), xpd=T)

    #Designing the colors of the barplot
    BarColor=c("#dca419","#cf22cb","#19dca4", "#08f2f3", "#ff7b1e", "#D3D3D3")
    names(BarColor)=c("AfAm", "Eur","Lat", "Dutch", "Fin", "Other")
    BorderColor = rep(BarColor[rownames(df)],ncol(df))
    BorderColor[which(as.vector(df)==0)] = NA

    #Generating barplots
    mp=barplot(as.matrix(df), beside=T, las=1, ylab="", xlab="", main="", font.main=1, col=BarColor[rownames(df)], axisnames=FALSE)
    title(xlab="Copy number", line=1.5)
    title(ylab=yLabel, line=2.5)
    mtext(text = colnames(df), side = 1, at = colMeans(mp), line = 0.5)
    
    legend(x=ncol(df)*nrow(df)+ncol(df), y=par("yaxp")[2], legend=rownames(df),col=BarColor[rownames(df)], pch=15, bty="n", title="Population ancestry:", y.intersp=1.25, pt.cex=1.25, title.adj=0)
}

makeStackedBarplot <- function(PopCounts, plotCounts) {
    if (plotCounts) {
        df = PopCounts
        title = "Count"
    } else {
        df = PopCounts / rowSums(PopCounts)
        title = "Frequency"
    }

    copyNumbers = colnames(df)
    copyNumberColors = colors[ifelse(is.na(copyNumbers), "NC", paste0("CN", copyNumbers))]

    par(mar=c(3, 5, 1, 5), xpd=T)

    df = t(df)
    barplot(df[, ncol(df):1], horiz=T, las=1, col=copyNumberColors, xlab="", ylab="",
        space=0,
        axes=F
    )
    axis(1, line=0)
    title(1, title, line=2)

    cnRange = 1:min(length(copyNumbers), 7)
    legend(
        x=par("xaxp")[2] + 0.04,
        y=par("yaxp")[2], 
        legend=ifelse(names(copyNumberColors) == "CN-1", "* allele", names(copyNumberColors)),
        col=copyNumberColors[cnRange],
        pch=15, 
        bty="n",
        x.intersp=0.6,
        y.intersp=0.8,
        pt.cex=1.25,
        title.adj=0.5
    )
}

