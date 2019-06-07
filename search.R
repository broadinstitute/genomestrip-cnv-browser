getSearchInterval <- function(chrom, startPos, endPos) {
    sprintf("%s:%s-%s", chrom, startPos, endPos)
}

padWindow <- function(window) {
    windowLength = window$END - window$START + 1
    windowMiddle = (window$START + window$END) / 2
    window$START = round(max(windowMiddle - 0.54 * windowLength, 1))
    window$END = round(min(windowMiddle + 0.54 * windowLength, chromData[window$CHROM, "END"]))
    return(window)
}

getIntervalWindow <- function(input) {
    s = strsplit(input, ":")
    chrom = s[[1]][1]
    chromEnd = chromData[chrom, "END"]
    if (is.na(chromEnd)) {
        return(list(ERROR_MSG = sprintf("There is no chromosome named %s", chrom)))
    } else if (length(s[[1]]) == 1) {
        coordinates = c(1, chromEnd)
    } else {
        coordinates = strsplit(s[[1]][2], "-")[[1]]
        if (length(coordinates) == 1) {
            coordinates = c(coordinates[1], chromEnd)
        }
        coordinates = suppressWarnings(as.numeric(coordinates))
        if (is.na(coordinates[1])) {
            return(list(ERROR_MSG = sprintf("Error: interval start is not numeric")))
        } else if (is.na(coordinates[2])) {
            return(list(ERROR_MSG = sprintf("Error: interval end is not numeric")))
        } else {
            if (coordinates[2] < coordinates[1]) {
                return(list(ERROR_MSG = "Error: interval start is larger than interval end"))
            }
            if (coordinates[2] > chromEnd) {
                coordinates = c(coordinates[1], chromEnd)
            }
        }
    }

    return(list(
        CHROM = chrom,
        START = coordinates[1],
        END = coordinates[2]))
}

getCallWindow <- function(input) {
    call = callData[callData$ID == input, ]
    if (nrow(call) == 0) {
        return(list(ERROR_MSG = sprintf("There is no call named %s", input)))
    }

    window = call[1, c("CHROM", "START", "END")]
    padWindow(window)
}

getGeneWindow <- function(input) {
    geneInfo = geneTxData[geneTxData$GENE == input, ]
    if (nrow(geneInfo) == 0) {
        return(list(ERROR_MSG = sprintf("There is no gene named %s", input)))
    }

    window = list(CHROM = geneInfo[1,"CHROM"], START = min(geneInfo$START), END = max(geneInfo$END))
    padWindow(window)
}

getSearchWindow <- function(input) {
    input = gsub("^\\s+|\\s+$", "", input)
    if (grepl(INTERVAL_REGEX, tolower(input))) {
        getIntervalWindow(input)
    } else if (grepl(CALL_REGEX, input)) {
        getCallWindow(input)
    } else {
        getGeneWindow(input)
    }
}

getOverlappingData <- function(df, searchWindow, filterConditions=NULL) {
    segmentStackingGap = round((searchWindow$END - searchWindow$START) / 100)

    searchCondition = "CHROM == searchWindow$CHROM & END >= searchWindow$START & START <= searchWindow$END"
    if (!is.null(filterConditions)) {
        searchCondition = paste(searchCondition, paste(filterConditions, collapse="&"), sep="&")
    }

    searchResult = subset(
        df,
        eval(parse(text=searchCondition))
    )

    if (nrow(searchResult) == 0) {
        searchResult = NULL
    } else {
        searchResult$OFFSET = computeSegmentOffsets(searchResult$START, searchResult$END, segmentStackingGap)
    }
    searchResult
}

searchCalls <- function(searchWindow, filterConditions=NULL) {
    writeLog(sprintf("searchCalls: %s:%d-%d", searchWindow$CHROM, searchWindow$START, searchWindow$END))

    callsDF = getOverlappingData(callData, searchWindow, filterConditions)
    segdupDF = getOverlappingData(segdupData, searchWindow)
    geneTxDF = getOverlappingData(geneTxData, searchWindow)
    geneExonDF = geneExonData[geneExonData$CHROM  == searchWindow$CHROM & geneExonData$END >= searchWindow$START & geneExonData$START <= searchWindow$END, ]

    writeLog("# callsDF: ", ifelse(is.null(callsDF), 0, nrow(callsDF)))

    if (!is.null(geneTxDF)) {
        writeLog("# genes: ", length(unique(geneTxDF$GENE)))
    }

    searchResults = list(
        callsDF = callsDF,
        segdupDF = segdupDF,
        geneTxDF = geneTxDF,
        geneExonDF = geneExonDF,
        window = searchWindow)

    return(searchResults)
}

searchGenes <- function(searchWindow) {
    writeLog(sprintf("searchGenes: %s:%d-%d", searchWindow$CHROM, searchWindow$START, searchWindow$END))

    geneTxDF = getOverlappingData(geneTxData, searchWindow)
    geneExonDF = geneExonData[geneExonData$TX_ID %in% rownames(geneTxDF), ]

    if (!is.null(geneTxDF)) {
        writeLog("# genes: ", length(unique(geneTxDF$GENE)))
    }

    searchResults = list(
        geneTxDF = geneTxDF,
        geneExonDF = geneExonDF,
        window = searchWindow)

    return(searchResults)
}

