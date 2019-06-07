require(shinysky)
require(shiny)
require(shinyFiles)
require(shinyjs)
require(shinyTypeahead)
require(colourpicker)

source("default_props.R")
source("custom/site_props.R")

getHeadTag <- function() {
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "CNVBrowser.css"),
    tags$link(rel="shortcut icon", href="gsfavicon-16x16.png", type="image/png"),
    tags$script(type="text/javascript", src="http://software.broadinstitute.org/software/genomestrip/misc/jquery.once.js?v=1.2"),
    #tags$script(type="text/javascript", src="http://software.broadinstitute.org/software/genomestrip/sites/all/modules/google_analytics/googleanalytics.js?nirlpt"),
    tags$head(tags$script(src="CNVBrowser.js"))
  )
}

createSearchArea <- function() {
  div(class="SearchArea",
      h1(siteName),
      typeaheadInput(inputId="search",label=NULL,choices=NULL,items=3,minLength=2),
      actionButton(inputId="searchButton",label="Search"),
      div(class="tinyText",
          p(class="examples"," Search by gene name ",
            actionButton(inputId="geneClick", label=paste0("(", defaultGene, ")"), class="hyperlink"),
            ", CNV identifier ",
            actionButton(inputId="CNVClick", label=paste0("(", defaultCall, ")"), class="hyperlink"),
            ", or hg38 genomic position ",
            actionButton(inputId="PositionClick",label=paste0("(", defaultInterval, ")"), class="hyperlink")
          )
      ),
      textOutput(outputId="searchWarnings")
  )
}

createFilterArea <- function() {
  list(
    div(id="searchFilterCheckBox", class="searchCheckBox",
        checkboxInput(inputId="filterSearch", label="Apply search filters")),
    fluidRow(id="searchFiltersArea", class="Hide",
             column(4,
                    radioButtons(inputId="filterType", label="CNV Type",
                                 choices=c("any" = "any",
                                           "Deletion" = "DEL",
                                           "Duplication" = "DUP",
                                           "Mixed" = "MIXED"),
                                 selected = "any"
                    )
             ),
             column(8,
                    sliderInput(inputId="filterFreq", label="Frequency (%)", min=0, max=100, value=c(0, 100), step=10)
             ),
             column(8,
                    sliderInput(inputId="filterCallrate", label="Callrate (%)", min=0, max=100, value=c(0, 100), step=10)
             )
    )
  )
}

createNavButtons <- function(plotName) {
  list(
    actionButton(inputId=paste0(plotName, "LeftButton3"),label="<<<",class="NoRight"),
    actionButton(inputId=paste0(plotName, "LeftButton2"),label="<<",class="NoRight"),
    actionButton(inputId=paste0(plotName, "LeftButton1"),label="<",class="NoRight"),
    actionButton(inputId=paste0(plotName, "RightButton1"), label=">", class="NoRight"),
    actionButton(inputId=paste0(plotName, "RightButton2"),label=">>",class="NoRight"),
    actionButton(inputId=paste0(plotName, "RightButton3"),label=">>>",class="NoRight"),
    "Zoom in",
    actionButton(inputId=paste0(plotName, "ZoomInButton1"), label="1.5x",class="NoRight"),
    actionButton(inputId=paste0(plotName, "ZoomInButton2"), label="3x",class="NoRight"),
    actionButton(inputId=paste0(plotName, "ZoomInButton3"), label="10x",class="NoRight"),
    "Zoom out",
    actionButton(inputId=paste0(plotName, "ZoomOutButton1"),label="1.5x",class="NoRight"),
    actionButton(inputId=paste0(plotName, "ZoomOutButton2"), label="3x",class="NoRight"),
    actionButton(inputId=paste0(plotName, "ZoomOutButton3"), label="10x",class="NoRight")
  )
}

createPlotControls <- function(plotName) {
  fluidRow(
    column(12, id="locusPlotNavDiv",
           #div(id="locusPlotNavDiv",
           div(class="locusWindowInterval",
               textOutput(outputId=paste0(plotName, "Interval"))
           ),
           div(class="ucscLink",
               uiOutput(outputId=paste0(plotName, ".UCSCLink"))
           ),
           div(class="navbuttons",
               createNavButtons(plotName)
               #),
           )
    )
  )
}

createLocusPlotArea <- function(plotName,  locusNumber) {
  div(id=paste0(plotName, "Area"), class="secondaryLocusArea",
      h1(id=paste0(plotName, "Heading"), sprintf("Locus %d", locusNumber), class="locusHeading"),
      createPlotControls(plotName),
      uiOutput(paste0(plotName, "_tooltip")),
      fluidRow(
        column(11,
               uiOutput(paste0(plotName, "UI"))
        ),
        column(1,
               plotOutput(outputId=paste0(plotName, "Legend"), height=150, width=80),
               tags$br(),
               uiOutput(outputId=paste0(plotName, "Height"))
        )
      )
  )
}

createSummaryResultsArea <- function() {
  plotMainHover = DEFAULT_HOVER
  plotMainHover$id = "locusPlotMain_hover"
  plotMainBrush = DEFAULT_BRUSH
  plotMainBrush$id = "locusPlotMain_brush"
  
  div(id="summaryResultsArea", class="Hide",
      createPlotControls("MainPlot"),
      uiOutput("locusPlotMain_tooltip"),
      fluidRow(id="mainPlotArea",
               column(11,
                      uiOutput(outputId="locusPlotMainUI")
               ),
               column(1,
                      plotOutput(outputId="MainPlotLegend", height=150,width=80),
                      tags$br(),
                      uiOutput(outputId="MainPlotHeight")
               )
      ),
      DT::dataTableOutput(outputId="searchResultsTable")
  )
}

createCallDetailsArea <- function() {
  list(
    h2("Basic Information"),
    fluidRow(
      column(6,
             tags$dl(
               tags$dt("Type:"),
               tags$dd(textOutput(outputId="callCategory")),
               tags$dt("Est. Length:"),
               tags$dd(textOutput(outputId="callLength")),
               tags$dt("Est. Interval:"),
               tags$dd(textOutput(outputId="callInterval")),
               tags$dt("Breakpoint:"),
               tags$dd(uiOutput(outputId="breakpoint"))
             )
      ),
      column(5,
             tags$dl(
               tags$dt("Callrate:"),
               tags$dd(textOutput(outputId="callrate")),
               tags$dt("Frequency:"),
               tags$dd(textOutput(outputId="callFreq")),
               tags$dt("CNV URL:"),
               tags$dd(uiOutput(outputId="callLink")),
               tags$dt("Gene(s):"),
               tags$dd(uiOutput(outputId="geneList"))
             )
      )
    ))
}

helpIcon <- function() {
    icon("question-sign", lib="glyphicon", class="help-icon")
}

createGenotypingPlots <- function() {
  div(id="genotypingPlotsArea",
      fluidRow(
        column(6,
               h3("Allele frequency by population"),
               div(class="xtooltip",
               actionLink(inputId="helpAlleleFrequencyTable", label="", icon=helpIcon()),
               p(class="tooltiptext", "Allele frequency help")),
               div(class="controls",
                   column(4, class="rightalign",
                          checkboxInput(inputId="showHaploidData", label="Haploid / Diploid", value=TRUE)),
                   column(4, class="centeralign",
                          checkboxInput(inputId="showPopulationAlleleTable", label="Table / plot", value=TRUE)),
                   column(4,
                          checkboxInput(inputId="showCountsData", label="Count / frequency", value=TRUE))
               )),
        column(5,
               h3("Sample genotypes"),
               div(class="xtooltip",
                   actionLink(inputId="helpSampleGenotypes", label="", icon=helpIcon()),
                   p(class="tooltiptext", "Sample Genotypes help")),
               div(class="controls",
                   radioButtons(inputId="showSampleGenotypes", inline=TRUE, label="Show samples:",
                                choices=c(
                                  "Carrier" = "carrier",
                                  "Non-ref" = "nonref",
                                  "All" = "all"),
                                selected = "carrier"
                   )))),
      fluidRow(
        column(6,
               uiOutput(outputId="populationAlleleUI", class="dataTable")
        ),
        column(5, class="tableholder",
               DT::dataTableOutput(outputId="perSampleGenotypesTable")
        ))
  )
}

createMainTab <- function() {
  list(
    createSearchArea(),
    div(id="searchInfoOverlay",
        class="Hide",
        createFilterArea(),
        createSummaryResultsArea()
    ),
    div(id="siteInfoOverlay",
        class="info",
        includeHTML("custom/www/site_description.html")
    )
  )
}

createDetailsTab <- function() {
  div(id="callDetailsTab",
      div(class="centered",
          div(class="history",
              "History:",
              actionButton(inputId="backButton", label="Back"),
              actionButton(inputId="forwardButton", label="Forward")
          ),
          textOutput(outputId="callName")),
      h2("Locus Map"),
      createLocusPlotArea("locusPlotOne", 1),
      createLocusPlotArea("locusPlotTwo", 2),
      div(id="callDetailsArea",
          createCallDetailsArea()
      ),
      h2("Genotypes"),
      createGenotypingPlots(),
      fluidRow(
        column(6,
               h3("Copy number across samples"),
               div(class="xtooltip",
                   actionLink(inputId="helpCopyNumberPlot", label="", icon=helpIcon()),
                   p(class="tooltiptext", "Copy number help")),
               div(class="controls",
                   column(7,
                          checkboxInput(inputId="greyNonConfidentCallss", label="Show non-confident calls in grey", value=TRUE),
                          checkboxInput(inputId="hideStarAlleleSamples", label="Hide star-allele samples", value=FALSE)
                   ),
                   column(4,
                          uiOutput(outputId="histHeightControllerUI"))
               )),
        column(5,
               h3("Copy number by genotyping batch"),
               div(class="xtooltip",
                   actionLink(inputId="helpCopyNumberBatchPlot", label="", icon=helpIcon()),
                   p(class="tooltiptext", "Copy number batch help")),
               div(class="controls batch",
                   uiOutput(outputId="batchSelector")))),
      fluidRow(
        column(6,
               plotOutput(outputId="genotypesHistogram", height=350)
        ),
        column(5,
               plotOutput(outputId="batchGenotypesHistogram", height=350)
        )
      ),
      
      h3(id="ProfilePlot","Batch profile plot"),
      div(class="xtooltip",
          actionLink(inputId="helpBatchProfilePlot", label="", icon=helpIcon()),
          p(class="tooltiptext", "Profile plot help")),
      uiOutput("profilePlot_tooltip"),
      div(id="profilePlotControls",
                   column(12,
                 div(class="controls",
                   column(2, class="smaller",
                          div(id="plotProfileMedianDiv",
                              checkboxInput(inputId="plotProfileMedian", value=FALSE, label="Plot median")
                          )),
                   column(4, class="smaller",
                          div(id="profileSamplesDiv",
                              textAreaInput(inputId="profileSamples", label="Highlight samples", value="", height="100%", resize="none")
                          ),
                          div(id="colorInputDiv",
                              colourInput(inputId="profileSampleColor", label="Color", value="blue", showColour="background")
                          ),
                          div(id="refreshProfilePlotDiv",
                              actionButton(inputId="refreshProfilePlot", label="Refresh profile plot")
                          )
                   ),
                   column(2,
                          div(id="profileBinningDiv",
                              textInput(inputId="profileBinning", label="Binning", value="", width="150px"),
                              checkboxInput(inputId="useAutoBinning", value=TRUE, label="Use auto-binning", width="150px"),
                              textOutput(outputId="profileBinSize")
                          )
                   ),
                   column(2,
                          sliderInput(inputId="profileMaxDepth", label="Max read depth", min=0, max=10, value=5, step=.2)
                   ),
                   column(2,
                              uiOutput(outputId="profileHeightControllerUI")
                          )
                   )
               )),
      fluidRow(
        column(12, id="profilePlotNavDiv",
               createNavButtons("ProfilePlot")
        )
      ),
      fluidRow(
        column(12, id="profilePlotDiv",
               uiOutput(outputId="profilePlotUI")
        )
      )
  )
}

createTabset <- function() {
  div(class="topnavbar",
      div(id="logoArea",
          class="logo",
          img(src="gs_logo_tag.png"),
          # this needs to be styled differently
          actionLink("homeLink","Genome STRiP",class="logo")
      ),
      div(id="navBarSearchWidget", class="Hide",
          div(class="SearchAreaNavBar",
              typeaheadInput(inputId="navBarSearch",label=NULL,choices=NULL,items=3,minLength=2),
              actionButton(inputId="navBarSearchButton",label="Search")
          )
      ),
      tabsetPanel(id="mainTabset", type= "tabs",
                  tabPanel(title="Search & Results", value="panelSearch",
                           createMainTab()
                  ),
                  tabPanel(title="Details", value="panelDetails",
                           createDetailsTab()
                  ),
                  tabPanel(title="Help", value="panelHelp",
                           div(id="Help", class="doc",
                               includeHTML("www/help.html")
                           )
                  ),
                  tabPanel(title="About", value="panelAbout",
                           div(id="About", class="doc",
                               includeHTML("www/about.html")
                           )
                  )
      )
  )
}


function(request) {fluidPage(title="CNV Browser | GenomeSTRiP", useShinyjs(),
                             
                             getHeadTag(),
                             tags$div(id="container",  
                                      createTabset()
                             ),
                             busyIndicator(wait = 3000)
)}

