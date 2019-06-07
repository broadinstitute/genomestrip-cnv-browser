# The shiny server log files don't seem to be reliably writing to 
# /home/unix/dkulp/bobh/shiny/var/log/shiny-server
zz <- file("app.log")
sink(zz, append=T)
sink(zz, append=T, type = "message")

#Server function
#### Loading Stuff In ####################################################

require(shiny)
require(shinyjs)
require(shinyFiles)
require(shinysky)
require(shinyTypeahead)
require(DT)

load("data/CallData.RData", envir=.GlobalEnv)
load("data/GeneData.RData", envir=.GlobalEnv)

sampleBatch <- read.table("data/sample_batch.dat", header=TRUE, sep="\t", stringsAsFactors=F)
batches <- paste0("batch", unique(sampleBatch$BATCH))

colors <- read.table("www/CNVcolors.txt", sep="\t", col.names=c("CN", "COLOR"), comment.char="", stringsAsFactors=F)
colors <- structure(colors$COLOR, names=colors$CN)

segdupData <- read.table("data/ucsc_segdups_hg38.dat", header=TRUE, sep="\t", stringsAsFactors=F)
segdupData <- segdupData[with(segdupData, order(CHROM, START, END)), ]

blacklistedRegions <- NULL
# Read the blacklisted regions file, if it exists
if (file.exists("data/blacklisted_regions.dat")) {
    blacklistedRegions <- read.table("data/blacklisted_regions.dat", header=TRUE, sep="\t", stringsAsFactors=F)
}

source("default_props.R", local=TRUE)
source("common.R", local=TRUE)
source("gene_plotting.R", local=TRUE)
source("locus_plotting.R", local=TRUE)
source("plotting.R", local=TRUE)
source("profile_plotting.R", local=TRUE)
source("search.R", local=TRUE)

shinyServer(function(input, output, session) {
    writeLog("ShinyServer")

    historyIdx <- 0
    callHistory <- c()

    observe ({
        query <- parseQueryString(session$clientData$url_search)

        if (!is.null(query$search)) {
            InputSearch$search = query$search
        } else if (!is.null(query$id)) {
            SelectedCall$ID = query$id
            SelectedCall$adjustWindow = TRUE

            updateTabsetPanel(session, "mainTabset", selected = "panelDetails")
        }
    })

    useShinyjs()

    source("nav_buttons.R", local=TRUE)
    source("reactiveValues.R", local=TRUE)
    source("reactives.R", local=TRUE)
    source("events.R", local=TRUE)
    source("outputs.R", local=TRUE)

    updateTypeaheadInput(session, inputId='search', choices=unname(unlist(choices)))
    updateTypeaheadInput(session, inputId='navBarSearch', choices=unname(unlist(choices)))
})

