getCategoryDescr <- function(categories) {
    ifelse(categories == "DEL", "Deletion",
        ifelse(categories == "DUP", "Duplication", "Mixed"))
}

getLocusPlotTooltip <- function(searchResults, plotLayout, mouseEvent) {
    x = mouseEvent$x * 1000
    y = mouseEvent$y

    plotSection = plotLayout[plotLayout$Y0 == max(plotLayout[y >= plotLayout$Y0, "Y0"]), ]
    tooltipLines = do.call(plotSection$TOOLTIP_FUNCTION, list(x, y, plotSection, searchResults))
    HTML(paste(tooltipLines, collapse='<br/>'))
}

getUCSCHref <- function(window) {
    if (is.null(window)) {
        return(NULL)
    }

    url = sprintf("%s&position=%s:%d-%d", UCSC_URL, window$CHROM, window$START, window$END)
    tags$a("View region in UCSC", title="View in UCSC Genome Browser", target="_blank", href=url)
}

createGenotypesHistogram <- function(data) {
    NUM_BINS = 100
    maxX = ceiling(max(0, data[!is.na(data$CNF), "CNF"]) + 0.5)
    binWidth = maxX / NUM_BINS
    binBreaks = seq(0, maxX, maxX / NUM_BINS)

    hist(data$CNF, breaks=binBreaks, plot=FALSE)
}

output$searchWarnings = renderText({
    MainPlotValues$window$ERROR_MSG
})

output$callName = renderText({
    getSelectedCall()$ID
})

output$callCategory = renderText({
    call = getSelectedCall()
    segdupText = ""
    if (!is.na(call$CHROM.2)) {
       segdupText = " (segdup)"
    }
    paste0(getCategoryDescr(getSelectedCall()$CATEGORY), segdupText)
})

output$callInterval = renderText({
    call = getSelectedCall()
    callInterval = formatCallInterval(call$CHROM, call$START, call$END)
    ifelse(is.na(call$CHROM.2),
        callInterval,
        paste0(callInterval, "; ", formatCallInterval(call$CHROM.2, call$START.2, call$END.2))
    )
})

output$callLength = renderText({
    sprintf("%s (%s alignable)", format(getSelectedCall()$LENGTH, nsmall=0, big.mark=","), format(getSelectedCall()$GSELENGTH, nsmall=0, big.mark=","))
})

output$callFreq = renderText({
    call = getSelectedCall()
    sprintf("%.3f%% (%d dels, %d dups)",
        100*call$FREQ, call$NDEL, call$NDUP)
})

output$callrate = renderText({
    paste0(round(100*getSelectedCall()$CALLRATE, 2), "%")
})

output$breakpoint = renderText({
    getSelectedCall()$GSBKPT
})

output$callLink = renderUI({
    clientData = session$clientData
    url = sprintf("%s//%s:%s%s?id=%s", clientData$url_protocol, clientData$url_hostname, clientData$url_port, clientData$url_pathname, getSelectedCall()$ID)
    tags$a(getSelectedCall()$ID, href=url, target="_blank")
})

output$geneList = renderUI({
    genes = getOverlappingGenes()
    if (length(genes) == 0) {
        return(NULL)
    }

    geneItems = lapply(1:length(genes), function(idx) {
        geneSeparator = ifelse(idx < length(genes), ",", "")
        geneUrl = sprintf("http://www.genecards.org/cgi-bin/carddisp.pl?gene=%s", genes[idx])
        geneLink = tags$a(href=geneUrl, paste0(genes[idx], geneSeparator), target="_blank")
        geneLink
    })

    geneItems
})

output$MainPlotInterval <- renderText({
    formatLocusInterval(MainPlotValues$window)
})

output$MainPlot.UCSCLink <- renderUI({
    getUCSCHref(MainPlotValues$window)
})
outputOptions(output, "MainPlot.UCSCLink", suspendWhenHidden = FALSE)

output$MainPlotHeight <- renderUI({
    sliderInput("mainPlotSlider", "Plot height:", min=100, max=1000, step=100, value=MainPlotValues$plotHeight, ticks=FALSE)
})

output$locusPlotMain <- renderPlot({
    writeLog("Render locusPlotMain")
    MainPlotValues$redisplay
    if (!is.null(searchResultsMain())) {
        plotLayout = createLocusPlotLayout(searchResultsMain())
        MainPlotValues$plotLayout = plotLayout
        makeLocusPlot(searchResultsMain(), plotLayout, par("fin")[2])
    }
})

output$locusPlotMainUI <- renderUI({
    plotMainHover = DEFAULT_HOVER
    plotMainHover$id = "locusPlotMain_hover"
    plotMainBrush = DEFAULT_BRUSH
    plotMainBrush$id = "locusPlotMain_brush"

    plotOutput("locusPlotMain", click="locusPlotMain_click", hover=plotMainHover, brush=plotMainBrush, height=MainPlotValues$plotHeight)
})

output$MainPlotLegend <- renderPlot({
    makeLosusLegend()
})

output$searchResultsTable = DT::renderDataTable({
    writeLog("Render searchResults table")
    if (is.null(searchResultsMain()$callsDF)) {
        return(NULL)
    }

    colNames = c("ID", "CHROM", "START", "END", "LENGTH", "CALLRATE", "NDEL", "NDUP", "FREQ", "CATEGORY")

    df = searchResultsMain()$callsDF[, colNames]
    df$CATEGORY = getCategoryDescr(df$CATEGORY)

    datatable(
        data = df,
        rownames = FALSE,
        selection = 'single',
        colnames = c("CNV name","Chr.","Start", "End", "Length", "Callrate (%)", "# Dels", "# Dups", "Freq. (%)", "Type"),
        extensions = 'Buttons',
        options = list(
            dom = 'l<"#callDownloadButton"B>"frtip',
            lengthMenu = c(25, 50, 100), pageLength = 25,
            buttons = list(
                list(extend="csv", text="Download", filename="gs_search_results.csv")
            )
        )
    ) %>%
        formatStyle('ID', color="#00609F", cursor="pointer") %>% 
        formatPercentage('CALLRATE', 2) %>% 
        formatPercentage('FREQ', 3) %>% 
        formatCurrency('LENGTH', '', digits=0)
}, server=TRUE)

output$genotypingPlotsHeading = renderText({
    paste0("Genotype data for ", getSelectedCall()$ID)
})

output$searchResultsText = renderText({
    writeLog("searchResultsText")
    if (!is.null(MainPlotValues$window) & is.null(searchResultsMain()$callsDF)) {
        "No CNVs Found"
    } else {
        NULL
    }
})

output$locusPlotOneInterval <- renderText({
    formatLocusInterval(PlotOneValues$window)
})

output$locusPlotOne.UCSCLink <- renderUI({
    getUCSCHref(PlotOneValues$window)
})
 
output$locusPlotOneHeight <- renderUI({
    sliderInput("locusPlotOneSlider", "Plot height:", min=100, max=1000, step=100, value=PlotOneValues$plotHeight, ticks=FALSE)
})

output$locusPlotOne <- renderPlot({
    PlotOneValues$redisplay
    if (!is.null(PlotOneValues$window)) {
        plotLayout = createLocusPlotLayout(searchResultsOne())
        PlotOneValues$plotLayout = plotLayout
        makeLocusPlot(searchResultsOne(), plotLayout, par("fin")[2], selectedCallID=SelectedCall$ID)
    }
})

outputLocusPlot <- function(plotName, plotClick, plotHeight) {
    plotHover = DEFAULT_HOVER
    plotHover$id = paste0(plotName, "_hover")
    plotBrush = DEFAULT_BRUSH
    plotBrush$id = paste0(plotName, "_brush")

    plotOutput(plotName, click=plotClick, hover=plotHover, brush=plotBrush, height=plotHeight)
}

output$locusPlotOneUI <- renderUI({
    plotClick = clickOpts(id="locusPlotOne_click")
    outputLocusPlot("locusPlotOne", plotClick, PlotOneValues$plotHeight)
})

output$locusPlotOneLegend <- renderPlot({
    makeLosusLegend()
})

output$locusPlotTwoHeight <- renderUI({
    sliderInput("locusPlotTwoSlider", "Plot height:", min=100, max=1000, step=100, value=PlotTwoValues$plotHeight, ticks=FALSE)
})

output$locusPlotTwo <- renderPlot({
    PlotTwoValues$redisplay
    if (!is.null(PlotTwoValues$window)) {
        plotLayout = createLocusPlotLayout(searchResultsTwo())
        PlotTwoValues$plotLayout = plotLayout
        makeLocusPlot(searchResultsTwo(), plotLayout, par("fin")[2], selectedCallID=SelectedCall$ID)
    }
})

output$locusPlotTwoUI <- renderUI({
    outputLocusPlot("locusPlotTwo", NULL, PlotTwoValues$plotHeight)
})

output$locusPlotTwoLegend <- renderPlot({
    makeLosusLegend()
})

output$locusPlotTwoInterval <- renderText({
    formatLocusInterval(PlotTwoValues$window)
})

output$locusPlotTwo.UCSCLink <- renderUI({
    getUCSCHref(PlotTwoValues$window)
})
 
# Histograms
output$genotypesHistogram = renderPlot({
    req(getGtData())

    gtData = getGtData()
    if (input$hideStarAlleleSamples) {
        gtData = gtData[gtData$ACN1 != -1 & gtData$ACN2 != -1, ]
    }

    baseHist = createGenotypesHistogram(gtData)

    plotGenotypesHistogram(gtData, baseHist, input$greyNonConfidentCallss, input$histHeightController)
}, height=350)

output$histHeightControllerUI <- renderUI({
    req(getGtData())

    baseHist = createGenotypesHistogram(getGtData())
    histHeight = max(baseHist$counts)
    sliderInput("histHeightController", "Y axis coordinates", min=10, max=histHeight, value=histHeight, step=10, ticks=FALSE, width=150)
})

output$populationAlleleUI <- renderUI({
    req(getGtData())

    if (input$showHaploidData && length(setdiff(unique(c(getGtData()$ACN1, getGtData()$ACN2)), -1)) > 2) {
        HTML("Allelic counts are unavailable, pending statistical phasing")
    } else if (input$showPopulationAlleleTable) {
        tableOutput(outputId="populationAlleleTable")
    } else {
        plotOutput(outputId="populationAlleleHistogram", height=350)
    }
})

output$populationAlleleTable= renderTable({
    writeLog("populationAlleleTable")
    req(getCopyNumbersByAncestry())

    tableData = getCopyNumbersByAncestry()
    columnPrefix = ifelse(input$showHaploidData, "ACN", "DCN")
    if (input$showCountsData) {
        colnames(tableData) = ifelse(colnames(tableData) == "-1", "Star", paste0(columnPrefix, colnames(tableData)))
    } else {
        tableData = round(100 * tableData / rowSums(tableData), 2)
        colnames(tableData) = ifelse(colnames(tableData) == "-1", "Star (%)", paste0(columnPrefix, colnames(tableData), " (%)"))
    }

    cbind("Ancestry" = rownames(tableData), tableData)
}, align="l")

output$populationAlleleHistogram = renderPlot({
    writeLog("populationAlleleHistogram")
    req(getCopyNumbersByAncestry())

    tableData = getCopyNumbersByAncestry()
    columnPrefix = ifelse(input$showHaploidData, "ACN", "DCN")
    colnames(tableData) = ifelse(colnames(tableData) == "-1", "Star", paste0(columnPrefix, colnames(tableData)))
    plotPopulationHistogram(tableData, input$showCountsData)
}, height=350)

output$locusPlotMain_tooltip <- renderUI({
    hover <- input$locusPlotMain_hover
    if (is.null(hover)) {
        return(NULL)
    }
    getLocusPlotTooltip(searchResultsMain(), MainPlotValues$plotLayout, hover)
})

output$locusPlotOne_tooltip <- renderUI({
    hover <- input$locusPlotOne_hover
    if (is.null(hover)) {
        return(NULL)
    }
    getLocusPlotTooltip(searchResultsOne(), PlotOneValues$plotLayout, hover)
})

output$locusPlotTwo_tooltip <- renderUI({
    hover <- input$locusPlotTwo_hover
    if (is.null(hover)) {
        return(NULL)
    }
    getLocusPlotTooltip(searchResultsTwo(), PlotTwoValues$plotLayout, hover)
})

# dev plots
output$perSampleGenotypesTable = DT::renderDataTable({
    req(getGtData())

    tableData = getPerSampleGenotypesTable()
    tableData[!is.na(tableData$ACN1) & tableData$ACN1 == -1, "ACN1"] = "*"
    tableData[!is.na(tableData$ACN2) & tableData$ACN2 == -1, "ACN2"] = "*"
    tableData$ACN <- paste(tableData$ACN1, tableData$ACN2, sep="+")

    selectedColumns = c("BATCH", "SAMPLE", "ACN", "CNF")
    columnNames = c("Batch", "Sample", "ACN", "CNF")

    if ("GSPC" %in% colnames(tableData)) {
        tableData[is.na(tableData$GSPC) | tableData$GSPC == ".", "GSPC"] = 0
        selectedColumns = c(selectedColumns, "GSPC", "GSPC2")
        columnNames = c(columnNames, "PC", "PC2")
    }

    tableData = tableData[, selectedColumns]
    datatable(
        data = tableData,
        rownames = FALSE,
        selection = 'single',
        colnames = columnNames,
        extensions = 'Buttons',
        options = list(
            dom = 'l<"#sampleGtDownloadButton"B>"frtip',
            lengthMenu = c(5, 10, 25, 50, 100), pageLength = 5,
            buttons = list(
                list(extend="csv", text="Download", filename="gs_sample_gts.csv")
            )
        )
    ) %>%
        formatStyle('BATCH', color="#00609F", cursor="pointer")
}, server=TRUE)

output$batchSelector <- renderUI({
    selectInput("genotypingBatch", "Genotyping Batch", batches, selected="batch1")
})

output$batchGenotypesHistogram = renderPlot({
    req(getGtData())
    req(input$genotypingBatch)
    writeLog("renderPlot batchGenotypesHistogram")

    batch = as.numeric(sub("batch", "", input$genotypingBatch))
    data = getGtData()
    data = data[data$BATCH == batch, ]

    baseHist = createGenotypesHistogram(data)
    plotGenotypesHistogram(data, baseHist, input$greyNonConfidentCallss, max(baseHist$counts))
}, height=350)

output$profilePlotUI = renderUI({
    writeLog("RenderUI profilePlotUI")

    plotClick = clickOpts(id="profilePlot_click")
    plotHover = DEFAULT_HOVER
    plotHover$id = "profilePlot_hover"
    plotBrush = DEFAULT_BRUSH
    plotBrush$id = "profilePlot_brush"

    plotOutput("profilePlot", click=plotClick, hover=plotHover, brush=plotBrush, height=ProfilePlotValues$plotHeight)
})
outputOptions(output, "profilePlotUI", priority=2)

output$profilePlot = renderPlot({
    writeLog("Render profilePlot")

    req(getProfilePlotData())
    req(getGtData())
    req(ProfilePlotValues$plotLayout)

    plotDims = isolate(ProfilePlotValues$plotDims)
    if (is.null(plotDims) || any(plotDims != par("fin"))) {
        writeLog("Resetting ProfilePlot plotDims")
        plotDims = par("fin")
        ProfilePlotValues$plotDims = plotDims
    }
    xRange = c(ProfilePlotValues$window$START, ProfilePlotValues$window$END)

    makeProfilePlot(getProfilePlotData(), plotDims, xRange, input$profileMaxDepth, input$plotProfileMedian, getBatchVariantSamples(), getSelectedProfileSamples(), ProfilePlotValues$sampleColor, getGtData(), isolate(ProfilePlotValues$plotLayout), par("fin")[2])
}, execOnResize=TRUE)
outputOptions(output, "profilePlot", priority=1)

output$profilePlot_tooltip <- renderUI({
    hover <- input$profilePlot_hover
    if (is.null(hover)) {
        return(NULL)
    }
    x = PROFILE_PLOT_XSCALE * hover$x
    y = hover$y

    plotLayout = isolate(ProfilePlotValues$plotLayout)
    plotSection = getPlotSection(y, plotLayout)
    if (is.null(plotSection)) {
        return()
    }

    tooltipLines = do.call(plotSection$TOOLTIP_FUNCTION, list(x, y, plotSection, getProfilePlotData()))
    HTML(paste(tooltipLines, collapse='<br/>'))
})

output$profileBinSize = renderText({
    paste("Using bin size:", getProfileBinSize())
})

output$profileHeightControllerUI <- renderUI({
    sliderInput("profileHeightController", "Profile plot height", min=100, max=1000, step=100, value=ProfilePlotValues$plotHeight, ticks=FALSE)
})


