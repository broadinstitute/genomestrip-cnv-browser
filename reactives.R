getFilterConditions <- function() {
    if (!input$filterSearch) {
        return(NULL)
    }

    conditions = c()
    if (input$filterType != "any") {
        conditions = c(conditions, sprintf("CATEGORY == '%s'", input$filterType))
    }

    if (input$filterFreq[1] > 0 || input$filterFreq[2] < 100) {
        freqCondition = sprintf("FREQ >= %f & FREQ <= %f", input$filterFreq[1]/100, input$filterFreq[2]/100)
        conditions = c(conditions, freqCondition)
    }
    if (input$filterCallrate[1] > 0 || input$filterCallrate[2] < 100) {
        callrateCondition = sprintf("CALLRATE >= %f & CALLRATE <= %f", input$filterCallrate[1]/100, input$filterCallrate[2]/100)
        conditions = c(conditions, callrateCondition)
    }

    conditions
}

getSelectedCall <- reactive({
    writeLog("reactive getSelectedCall()")

    if (is.null(SelectedCall$ID)) {
        return(NULL)
    }

    callData[callData$ID == SelectedCall$ID, ][1, ]
})

searchResultsMain <- reactive({
    #MainPlotValues$window

    writeLog("reactive MainPlotValues$window, generating searchResults")

    if (is.null(MainPlotValues$window)) {
        return(NULL)
    }

    searchCalls(MainPlotValues$window, getFilterConditions())
})

searchResultsOne <- reactive({
    #PlotOneValues$window

    writeLog("reactive PlotOneValues$searchResults")
    if (is.null(PlotOneValues$window)) {
        NULL
    } else {
        searchCalls(PlotOneValues$window)
    }
})
 
searchResultsTwo <- reactive({
    #PlotTwoValues$window

    writeLog("reactive PlotTwoValues$searchResults")
    if (is.null(PlotTwoValues$window)) {
        NULL
    } else {
        searchCalls(PlotTwoValues$window)
    }
})
 
getProfileGeneData <- reactive({
    #ProfilePlotValues$window

    writeLog("reactive getProfileGeneData")
    if (!is.null(ProfilePlotValues$window)) {
        searchGenes(ProfilePlotValues$window)
    }
})

getGtData <- reactive({
    writeLog("reactive getGtData()")
    req(getSelectedCall())

    call = getSelectedCall()
    if (is.null(call)) {
        writeLog("call = NULL")
        return(NULL)
    }

    fileName = paste0("data/seq_", call$CHROM, "/", call$ID, ".Genotypes.txt")
    data = read.table(fileName, header=T, sep="\t", stringsAsFactors=F)
    data = merge(data, sampleBatch, by="SAMPLE")
    data
})
 
getCopyNumbersByAncestry <- reactive({
    writeLog("reactive getCopyNumbersByAncestry()")

    data = getGtData()
    if (input$showHaploidData) {
        tableData = rbind(
            data.frame(POP=data$POP, ACN=data$ACN1, stringsAsFactors=F),
            data.frame(POP=data$POP, ACN=data$ACN2, stringsAsFactors=F)
        )
    } else {
        tableData = data.frame(POP=data$POP, DCN=ifelse(data$ACN1 == -1 | data$ACN2 == -1, -1, data$ACN1 + data$ACN2), stringsAsFactors=F)
    }
    cnSummary = as.data.frame.matrix(table(tableData))
    colNames = colnames(cnSummary)
    colOrder = as.numeric(colNames)
    colOrder = ifelse(colOrder == -1, NA, colOrder)
    cnSummary = cnSummary[ , colNames[order(colOrder)]]
})
 
getPerSampleGenotypesTable <- reactive({
    tableData = getGtData()

    contentType = input$showSampleGenotypes
    if (contentType == "carrier") {
        tableData = tableData[tableData$ACN1 != 1 & tableData$ACN1 != -1 | tableData$ACN2 != 1 & tableData$ACN2 != -1, ]
    } else if (contentType == "nonref") {
        tableData = tableData[tableData$ACN1 != 1 | tableData$ACN2 != 1, ]
    }
    tableData[with(tableData, order(BATCH, CN)), ]
})

getOverlappingGenes <- reactive({
    df = callData[callData$ID == SelectedCall$ID, ]
    if (nrow(df) == 0) {
        return(NULL)
    }

    searchCondition = do.call(paste, c(lapply(1:nrow(df), function(idx) {
        sprintf("CHROM == '%s' & END >= %d & START <= %d", df[idx, "CHROM"], df[idx, "START"], df[idx, "END"])
    }), sep=" | "))

    genes = subset(
        geneTxData,
        eval(parse(text=searchCondition))
    )

    unique(genes$GENE)
})

getProfilePlotData <- reactive({
    req(getSelectedCall())

    list(
        profileData = getProfileData(),
        geneTxDF = getProfileGeneData()$geneTxDF,
        geneExonDF = getProfileGeneData()$geneExonDF,
        discoveryPairData = getDiscoveryPairData(),
        genotypesPairData = getGenotypesPairData(),
        window = ProfilePlotValues$window,
        selectedCall = getSelectedCall())
})

getProfileBinSize <- reactive({
    windowSize = ProfilePlotValues$window$END - ProfilePlotValues$window$START + 1
    avalableBinSizes = PROFILE_DIRS$BIN_SIZE

    whichBins = which(windowSize / avalableBinSizes <= 500)
    if (length(whichBins) == 0) {
        binSize = tail(avalableBinSizes, 1)
    } else {
        binSize = avalableBinSizes[min(whichBins)]
    }
    binSize
})

loadProfileData <- reactive({
    writeLog("reactive loadProfileData()")

    req(input$genotypingBatch, ProfilePlotValues$window)

    binSize = getProfileBinSize()
    profilesPath = PROFILE_DIRS[PROFILE_DIRS$BIN_SIZE == binSize, "PATH"]

    profileFile = sprintf("%s/md_%s/profile_seq_%s_%d.dat.gz", profilesPath, input$genotypingBatch, ProfilePlotValues$window$CHROM, binSize)

    readProfileDataForRegion(profileFile, ProfilePlotValues$window)
})

getProfileData <- reactive({
    writeLog("reactive getProfileData()")
    
    profileData = loadProfileData()
    binsToMerge = ProfilePlotValues$binsToMerge
    if (input$useAutoBinning || is.null(binsToMerge)) {
        numBins = length(profileData$BINSTARTS)
        binsToMerge = getBinsToMerge(numBins, "auto")
        updateTextInput(session, inputId="profileBinning", value=binsToMerge)
    }

    if (binsToMerge > 1) {
        profileData = mergeProfileBins(profileData, binsToMerge, 0)
    }

    coverage = profileData$COVNUMERATORS / profileData$COVDENOMINATORS
    profileData$COVERAGE = coverage

    profileData
})

getBatchVariantSamples <- reactive({
    batch = as.numeric(sub("batch", "", input$genotypingBatch))

    gtData = getGtData()
    gtData[gtData$BATCH == batch & gtData$CN != -1 & gtData$CN != 2, "SAMPLE"]
})

getSelectedProfileSamples <- reactive({
    batch = as.numeric(sub("batch", "", input$genotypingBatch))
    if (is.null(ProfilePlotValues$selectedSamples[batch])) {
        selectedSamples = c()
    } else {
        selectedSamples = trimws(strsplit(ProfilePlotValues$selectedSamples[batch], ";")[[1]])
    }
})

getPairSupportData <- function(dataFile, samples, siteId=NULL) {
    if (!file.exists(dataFile) || is.null(samples)) {
        return(NULL)
    }

    call = getSelectedCall()
    pairSupportData = readTabixedData(dataFile, call)
    if (nrow(pairSupportData) == 0) {
        return(NULL)
    }

    pairSupportData = pairSupportData[pairSupportData$SAMPLE %in% samples, ]
    if (!is.null(siteId)) {
        pairSupportData = pairSupportData[pairSupportData$ID == siteId, ] 
    }
    if (nrow(pairSupportData) == 0) {
        return(NULL)
    }

    segmentStackingGap = round(call$END - call$START + 1) / PROFILE_PLOT_XSCALE / 100
    pairSupportData$OFFSET <- computeSegmentOffsets(pairSupportData$LEFTSTART, pairSupportData$RIGHTEND, segmentStackingGap)
    pairSupportData
}

getDiscoveryPairData <- reactive({
    writeLog("reactive getDiscoveryPairData")
    samples = getSelectedProfileSamples()
    if (length(samples) == 0) {
        return(NULL)
    }
    dataFile = sprintf("data/del_aux_data/%s/md_%s/gs_dels.discovery.pairs.gz", getSelectedCall()$CHROM, input$genotypingBatch)

    getPairSupportData(dataFile, samples)
})

getGenotypesPairData <- reactive({
    writeLog("reactive getGenotypesPairData")
    samples = getSelectedProfileSamples()
    if (length(samples) == 0) {
        return(NULL)
    }
    dataFile = sprintf("data/del_aux_data/%s/md_%s/gs_dels.genotypes.pairs.gz", getSelectedCall()$CHROM, input$genotypingBatch)

    getPairSupportData(dataFile, samples, SelectedCall$ID)
})

