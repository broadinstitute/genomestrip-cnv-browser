INTERVAL_REGEX = "^(chr){0,1}(\\d{1,2}|X|Y):(\\d+)-(\\d+)$"
CALL_REGEX = "^(SEG_|DEL_|GS_SD_)"

DEL_COLOR = "#b00b35"
DUP_COLOR = "#0ee912"
MIXED_COLOR = "black"
SELECTED_COLOR="black"

UCSC_URL = "http://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOtherUser=submit&hgS_otherUserName=usher&hgS_otherUserSessionName=CNV_Browser&db=hg38"

LOCUS_PLOT_HEIGHT <<- 400
LOCUS_PLOT_SCALE <<- 1000

PROFILE_PLOT_HEIGHT <<- 500

DEFAULT_HOVER <<- hoverOpts(id="", delay=100, delayType="throttle")

DEFAULT_BRUSH <<- brushOpts(id="", fill = "steelblue", stroke = "#036", opacity = 0.4, delay = 1000, delayType = "debounce", clip = TRUE, direction = "x", resetOnNew = TRUE)

# If the brush window's fraction of the whole locus window is less than this value, Ignore the brush event
BRUSH_FRACTION = 0.0125

PROFILE_PLOT_XSCALE = 1000

APP_DATA_DIR = "/stanley/genome_strip/wgspd"

PROFILE_DIRS = data.frame(
    BIN_SIZE = 100 * c(1, 10, 100, 1000),
    PATH = paste0(APP_DATA_DIR, "/", "profiles_", c("100", "1000", "10000", "100Kb")),
    stringAsFactors=FALSE
)

