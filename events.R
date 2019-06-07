
setInputSearch <- function(search) {
    InputSearch$search = NULL
    InputSearch$search = search
    if (!is.null(search)) {
        writeLog(sprintf("Set input search: %s", search))
        # Force search panel to display
        toggleVisibility("searchInfoOverlay", TRUE)
        updateTabsetPanel(session, "mainTabset", selected = "panelSearch")
    } else {
        writeLog("Set input search: NULL")
    }
}

observeEvent(input$homeLink,{
    writeLog("Home link!")
    updateTabsetPanel(session, "mainTabset", selected = "panelSearch")
    toggleVisibility("searchInfoOverlay", FALSE);
    toggleVisibility("siteInfoOverlay", TRUE);
})

observe ({
    toggleVisibility("navBarSearchWidget", input$mainTabset != "panelSearch")
})

observeEvent(input$geneClick,{
    setInputSearch(defaultGene)
    updateTypeaheadInput(session, inputId='search', value=defaultGene, choices=unname(unlist(choices))) 
    writeLog(sprintf("Gene %s clicked", defaultGene))
})

observeEvent(input$CNVClick,{
    setInputSearch(defaultCall)
    updateTypeaheadInput(session,inputId='search',value=defaultCall,choices=unname(unlist(choices))) 
    writeLog(sprintf("CNV %s clicked", defaultCall))
})

observeEvent(input$PositionClick,{
    setInputSearch(defaultInterval)
    updateTypeaheadInput(session,inputId='search',value=defaultInterval,choices=unname(unlist(choices)))
    writeLog(sprintf("Interval %s clicked", defaultInterval))
})

observeEvent(input$searchButton, {
    writeLog("observeEvent input$searchButton, updating Input$Search")
    if (input$search != "") {
        setInputSearch(input$search)
    }
})

observeEvent(input$navBarSearchButton, {
    writeLog("observeEvent input$navBarSearchButton, updating Input$Search")
    if (input$navBarSearch != "") {
        updateTypeaheadInput(session, inputId='search', value=input$navBarSearch, choices=unname(unlist(choices)))
        setInputSearch(input$navBarSearch)
    }
})

observeEvent(input$filterSearch, {
    toggleVisibility("searchFiltersArea", input$filterSearch)
    setInputSearch(InputSearch$search)
})

observeEvent(input$filterType, {
    setInputSearch(InputSearch$search)
})

observeEvent(InputSearch$search,{
    writeLog("observeEvent InputSearch$search")
    MainPlotValues$window = getSearchWindow(InputSearch$search)
    SelectedCall$ID = NULL

    toggleVisibility("summaryResultsArea", is.null(MainPlotValues$window$ERROR_MSG))
})

observe({
    searchResultsMain()
    writeLog("observe searchResultsMain()")

    if (is.null(searchResultsMain()$callsDF)) {
        writeLog("No calls, setting PlotOneValues to NULL")
        PlotOneValues$window = NULL
        ProfilePlotValues$window = NULL
    }
})

# Main locusPlot events
observeEvent(input$locusPlotMain_click,{
    writeLog("x,y coordinates =", input$locusPlotMain_click$x, ",", input$locusPlotMain_click$y)
    mouseItem = getLocusPlotMouseItem(searchResultsMain(), MainPlotValues$plotLayout, input$locusPlotMain_click)
    writeLog(sprintf("Clicked %s %s", mouseItem$type, mouseItem$name))
    if (!is.null(mouseItem) && mouseItem$type == "CALL" ) {
        SelectedCall$ID = mouseItem$name
        SelectedCall$adjustWindow = TRUE

        updateTabsetPanel(session, "mainTabset", selected = "panelDetails")
    }
})
 
processBrushEvent <- function(brushEvent, plotValues, plotScale) {
    start = round(brushEvent$xmin*plotScale)
    end = round(brushEvent$xmax*plotScale)
    if ((end - start)/ (plotValues$window$END - plotValues$window$START) >= BRUSH_FRACTION) {
        plotValues$window$START = start
        plotValues$window$END = end
    } else {
        plotValues$redisplay = NULL
        plotValues$redisplay = TRUE
    }
}

observeEvent(input$locusPlotMain_brush,{
    processBrushEvent(input$locusPlotMain_brush, MainPlotValues, LOCUS_PLOT_SCALE)
})

#  locusPlotOne events
observeEvent(input$locusPlotOne_click,{
    mouseItem = getLocusPlotMouseItem(searchResultsOne(), PlotOneValues$plotLayout, input$locusPlotOne_click)
    writeLog(sprintf("Clicked %s %s", mouseItem$type, mouseItem$name))
    if (!is.null(mouseItem) && mouseItem$type == "CALL" ) {
        SelectedCall$ID = mouseItem$name
        SelectedCall$adjustWindow = FALSE
    }
})
 
observeEvent(input$locusPlotOne_brush,{
    processBrushEvent(input$locusPlotOne_brush, PlotOneValues, LOCUS_PLOT_SCALE)
})
 
#  locusPlotTwo events
observeEvent(input$locusPlotTwo_click,{
    mouseItem = getLocusPlotMouseItem(searchResultsTwo(), PlotTwoValues$plotLayout, input$locusPlotTwo_click)
    writeLog(sprintf("Clicked %s %s", mouseItem$type, mouseItem$name))
    if (!is.null(mouseItem) && mouseItem$type == "CALL" ) {
        selectedCall = callData[callData$ID == mouseItem$name, ][1, ]
        if (nrow(selectedCall) > 0) {
            window2 = padWindow(selectedCall[, c("CHROM", "START", "END")])
            PlotTwoValues$window = window2
        }
    }
})
 
observeEvent(input$locusPlotTwo_brush,{
    processBrushEvent(input$locusPlotTwo_brush, PlotTwoValues, LOCUS_PLOT_SCALE)
})
 
observeEvent(input$searchResultsTable_rows_selected, {
    writeLog("observeEvent searchResultsTable_rows_selected: ", input$searchResultsTable_rows_selected)

    SelectedCall$ID = isolate(searchResultsMain()$callsDF[input$searchResultsTable_rows_selected, "ID"])
    SelectedCall$adjustWindow = TRUE

    updateTabsetPanel(session, "mainTabset",
        selected = "panelDetails")
})

observe({
    PlotOneValues$window
    writeLog("observe PlotOneValues$window changed to ", PlotOneValues$window)

    showCallDetails = !is.null(PlotOneValues$window)
    toggleVisibility("locusPlotOneArea", showCallDetails)
})

observe({
    PlotTwoValues$window
    writeLog("observe PlotTwoValues$window changed to ", PlotTwoValues$window)

    toggleVisibility("locusPlotTwoArea", !is.null(PlotTwoValues$window))
})

#### Button functions ####

shiftLeft <- function(window, fraction) {
    windowLength = window$END - window$START + 1
    window$START = max(window$START - round(fraction * windowLength), 1)
    window$END = window$START + windowLength - 1
    window
}

shiftRight <- function(window, fraction) {
    windowLength = window$END - window$START + 1
    window$END = min(window$END + round(fraction * windowLength), chromData[window$CHROM, "END"])
    window$START = window$END - windowLength + 1
    window
}

zoomIn <- function(window, fraction) {
    windowMiddle = (window$START + window$END) / 2
    windowLength = (window$END - window$START + 1) / fraction
    window$START = round(windowMiddle - windowLength / 2)
    window$END = round(windowMiddle + windowLength / 2)
    return(window)
}

zoomOut <- function(window, fraction) {
    windowMiddle = (window$START + window$END) / 2
    windowLength = fraction * (window$END - window$START + 1)

    window$START = round(max(1, windowMiddle - windowLength / 2))
    window$END = round(min(window$START + windowLength - 1, chromData[window$CHROM, "END"]))
    if (window$END == chromData[window$CHROM, "END"]) {
        window$START = max(1, window$END - windowLength + 1)
    }
    return(window)
}

updateHistoryButtons <- function() {
    if (historyIdx > 1) {
        enable("backButton")
    } else {
        disable("backButton")
    }
    if (historyIdx < length(callHistory)) {
        enable("forwardButton")
    } else {
        disable("forwardButton")
    }
}

addCallToHistory <- function(callId) {
    writeLog("addCallToHistory")
    if (is.null(callId)) {
        return()
    }
    if (historyIdx > 0 && callHistory[historyIdx] == callId) {
        return()
    }

    if (historyIdx > 0 && historyIdx < length(callHistory)) {
        callHistory <<- callHistory[1:historyIdx]
    }
    callHistory <<- c(callHistory, callId)
    historyIdx <<- length(callHistory)

    updateHistoryButtons()
}

observe({
    SelectedCall$ID

    addCallToHistory(SelectedCall$ID)

    if (is.null(SelectedCall$ID)) {
         hideTab(inputId="mainTabset", target="panelDetails")

        PlotOneValues$window = NULL
        PlotTwoValues$window = NULL
        ProfilePlotValues$window = NULL

        return(NULL)
    }

    showTab(inputId="mainTabset", target="panelDetails")

    selectedCall = callData[callData$ID == SelectedCall$ID, ][1, ]

    if (SelectedCall$adjustWindow) {
        window1 = padWindow(selectedCall[, c("CHROM", "START", "END")])
        PlotOneValues$window = window1
    }
    ProfilePlotValues$window = zoomOut(selectedCall[, c("CHROM", "START", "END")], 10)

    if (!is.na(selectedCall$CHROM.2)) {
        window2 = list(CHROM=selectedCall$CHROM.2, START=selectedCall$START.2, END=selectedCall$END.2)
        window2 = padWindow(window2)
        PlotTwoValues$window = window2

        toggleVisibility("locusPlotOneHeading", TRUE)
    } else {
        PlotTwoValues$window = NULL

        toggleVisibility("locusPlotOneHeading", FALSE)
    }
})

toggleVisibility <- function(id, isShow) {
    writeLog(sprintf("toggleVisibility for %s to %s", id, isShow))
    if (isShow) {
        addClass(id, "show")
        removeClass(id, "Hide")
    } else {
        removeClass(id, "show")
        addClass(id, "Hide")
    }
}

observeEvent(input$mainPlotSlider, {
    writeLog("observeEvent input$mainPlotSlider")
    MainPlotValues$plotHeight = input$mainPlotSlider
})

observeEvent(input$locusPlotOneSlider, {
    PlotOneValues$plotHeight = input$locusPlotOneSlider
})

observeEvent(input$locusPlotTwoSlider, {
    PlotTwoValues$plotHeight = input$locusPlotTwoSlider
})

observeEvent(input$backButton, {
    writeLog("observeEvent input$backButton")
    if (length(callHistory) > 0 & historyIdx > 1) {
        historyIdx <<- historyIdx - 1
        SelectedCall$ID = callHistory[historyIdx]
    }

    updateHistoryButtons()
})

observeEvent(input$forwardButton, {
    writeLog("observeEvent input$forwardButton")
    if (length(callHistory) > 0 & historyIdx < length(callHistory)) {
        historyIdx <<- historyIdx + 1
        SelectedCall$ID = callHistory[historyIdx]
    }

    updateHistoryButtons()
})

observeEvent(input$genotypingBatch, {
    req(getGtData())
    gtData = getGtData()$Data

    batch = as.numeric(sub("batch", "", input$genotypingBatch))
    updateTextAreaInput(session, "profileSamples", value=ProfilePlotValues$selectedSamples[batch])
})

observeEvent(input$perSampleGenotypesTable_rows_selected, {
    sample = getPerSampleGenotypesTable()[input$perSampleGenotypesTable_rows_selected, "SAMPLE"]
    batch = getPerSampleGenotypesTable()[input$perSampleGenotypesTable_rows_selected, "BATCH"]

    updateSelectInput(session, "genotypingBatch", selected=paste0("batch", batch))
    updateTextAreaInput(session, "profileSamples", value=sample)
    ProfilePlotValues$selectedSamples[batch] = sample
})

createProfilePlotLayout <- function() {
    plotLayout = data.frame(
        SECTION_NAME = c("Gene", "Profile", "Discovery", "Genotypes"),
        SECTION_TITLE = c("Gene", "Read Profile", "Discovery", "Pairs"),
        SCALED_HEIGHT = c(.3, .5, .1, .1),
        DIRECTION = c(-1, 1, 1, 1),
        TOOLTIP_FUNCTION = c("getGeneTooltip", "getProfileTooltip", "getDiscoveryPairTooltip", "getGenotypesPairTooltip"),
        stringsAsFactors=F
    )
    rownames(plotLayout) = plotLayout$SECTION_NAME

    plotLayout
}

observe({
    writeLog("observe getProfilePlotData()")

    plotData = getProfilePlotData()
    geneSectionHeight = getPlotSectionHeight(plotData$geneTxDF)
    profileSectionHeight = input$profileMaxDepth
    discoveryPairSectionHeight = getPlotSectionHeight(plotData$discoveryPairData)
    genotypesPairSectionHeight =  getPlotSectionHeight(plotData$genotypesPairData)

    plotLayout = ProfilePlotValues$plotLayout
    if (is.null(plotLayout)) {
        plotLayout = createProfilePlotLayout()
    }

    plotLayout$HEIGHT <- c(geneSectionHeight, profileSectionHeight, discoveryPairSectionHeight, genotypesPairSectionHeight)
    plotLayout$Y0 <- cumsum(c(0, head(plotLayout[, "SCALED_HEIGHT"], -1)))
    plotLayout$SCALE <- plotLayout$HEIGHT / plotLayout$SCALED_HEIGHT

    ProfilePlotValues$plotLayout = plotLayout
}, priority=10)

rescaleProfilePlotSection <- function(plotSection, scale) {
    plotLayout = ProfilePlotValues$plotLayout
    plotLayout[plotSection$SECTION_NAME, "SCALED_HEIGHT"] = max(.02, scale * plotLayout[plotSection$SECTION_NAME, "SCALED_HEIGHT"])
    plotHeight = ProfilePlotValues$plotHeight * sum(plotLayout$SCALED_HEIGHT)
    plotLayout$SCALED_HEIGHT = plotLayout$SCALED_HEIGHT / sum(plotLayout$SCALED_HEIGHT)

    ProfilePlotValues$plotHeight = plotHeight
    ProfilePlotValues$plotLayout = plotLayout
}

observeEvent(input$profilePlot_click,{
    x = input$profilePlot_click$x
    y = input$profilePlot_click$y

    plotSection = getPlotSection(y, ProfilePlotValues$plotLayout)
    if (is.null(plotSection) || plotSection$SECTION_NAME == "Gene") {
        return()
    }

    if (plotSection$SECTION_NAME %in% c("Discovery", "Genotypes")) {
        xRange = c(ProfilePlotValues$window$START, ProfilePlotValues$window$END)
        button = processAuxDataClick(x, y, ProfilePlotValues$plotDims, xRange, plotSection)

        if (is.null(button)) {
            return()
        } else if (button == "plus") {
            rescaleProfilePlotSection(plotSection, 2)
        } else if (button == "minus") {
            rescaleProfilePlotSection(plotSection, .5)
        }
  
        return()
    }

    x = x * PROFILE_PLOT_XSCALE
    y = (y - plotSection$Y0) * plotSection$SCALE
    selectedSamples = getProfilePlotSamples(x, y, getProfilePlotData())
    if (length(selectedSamples) > 0) {
        selectedSamples = paste0(selectedSamples, collapse=";")
        updateTextAreaInput(session, "profileSamples", value=selectedSamples)

        batch = as.numeric(sub("batch", "", isolate(input$genotypingBatch)))
        ProfilePlotValues$selectedSamples[batch] = selectedSamples
    }
})
 
observeEvent(input$profilePlot_brush,{
    processBrushEvent(input$profilePlot_brush, ProfilePlotValues, PROFILE_PLOT_XSCALE)
})
 
observeEvent(input$profileSampleColor,{
    ProfilePlotValues$sampleColor = input$profileSampleColor
})

observeEvent(input$refreshProfilePlot,{
    batch = as.numeric(sub("batch", "", isolate(input$genotypingBatch)))
    ProfilePlotValues$selectedSamples[batch] = input$profileSamples

    if (!input$useAutoBinning) {
        ProfilePlotValues$binsToMerge = input$profileBinning
    }
})

observeEvent(input$useAutoBinning,{
    if (!input$useAutoBinning) {
        ProfilePlotValues$binsToMerge = input$profileBinning
    }
})

observeEvent(input$profileHeightController, {
    if (ProfilePlotValues$plotHeight != input$profileHeightController) {
        ProfilePlotValues$plotHeight = input$profileHeightController
    }
})

helpTargetId <- NULL

observeEvent(input$mainTabset, {
    # Watch for switches to the help page and move to the requested topic (if any).
    # This code navigates to help in the same window. We should experiment with using a popup window instead.
    if (input$mainTabset == "panelHelp") {
        if (!is.null(helpTargetId)) {
            #writeLog(sprintf("Moving to help panel, targetID = %s", helpTargetId))
            runjs(sprintf("window.location.replace('#%s');", helpTargetId));
            helpTargetId <<- NULL
        }
    }
})

# I could not find a simple way to make these methods data driven (e.g. by a property on the element).

observeEvent(input$helpAlleleFrequencyTable, {
    helpTargetId <<- "help-topic-allele-frequency-table"
    updateTabsetPanel(session, "mainTabset", selected = "panelHelp")
})
observeEvent(input$helpSampleGenotypes, {
    helpTargetId <<- "help-topic-sample-genotypes"
    updateTabsetPanel(session, "mainTabset", selected = "panelHelp")
})
observeEvent(input$helpCopyNumberPlot, {
    helpTargetId <<- "help-topic-copy-number-plot"
    updateTabsetPanel(session, "mainTabset", selected = "panelHelp")
})
observeEvent(input$helpCopyNumberBatchPlot, {
    helpTargetId <<- "help-topic-copy-number-batch-plot"
    updateTabsetPanel(session, "mainTabset", selected = "panelHelp")
})
observeEvent(input$helpBatchProfilePlot, {
    helpTargetId <<- "help-topic-batch-profile-plot"
    updateTabsetPanel(session, "mainTabset", selected = "panelHelp")
})
