$(document).keyup(function(event) {
    if ($("#search").is(":focus") && (event.keyCode == 13)) {
        $("#searchButton").click();
    }
});

$(document).keyup(function(event) {
    if ($("#navBarSearch").is(":focus") && (event.keyCode == 13)) {
        $("#navBarSearchButton").click();
    }
});

$(document).keyup(function(event) {
    if ($("#profileBinning").is(":focus") && (event.keyCode == 13)) {
        $("#refreshProfilePlot").click();
    }
});

$(document).on("shiny:value",function(event){
    if(event.name === "CNVName"){
        console.log("CNVName happened!");
        $(".Population").css("float","none");
        $(".Population").css("margin","auto");
        if($("#CNVTwoPlot").css("display")=="none"){
            console.log("dup hidden");
            $("div.Centerable.SecondaryMap-Area").addClass("Centered")
            $("button.Mini.Resize").removeClass("Mini")
        }
        if($("#CNVTwoPlot").css("display")=="block"){
            console.log("dup showing");
            $("div.Centerable.SecondaryMap-Area.Centered").removeClass("Centered")
            $("button.Resize").addClass("Mini")
        }
    }
});

$(document).on("shiny:value",function(event){
    if(event.name === "MainPlotInterval"){
    	 console.log("Something happened!");
    $(".info").css("display","none");
    }});
    	
    	
$(document).on("shiny:value",function(event){
    if(event.name === "spiderplotGraph"){
        console.log("spiderplotchange");
        $(".Population").css("float","left");
        $(".Population").css("margin-top","6em");
    }
});

$(document).ready(function() {
    $("#search").focus();

    $("#locusPlotMainUI").mousemove(function(e) {
        $("#locusPlotMain_tooltip").show();         
        $("#locusPlotMain_tooltip").css({
            top:  (e.pageY + 25) + "px",
            left: (e.pageX - 100) + "px"
        });
    });

    $("#locusPlotOneUI").mousemove(function(e) {
        $("#locusPlotOne_tooltip").show();         
        $("#locusPlotOne_tooltip").css({
            top:  (e.pageY + 25) + "px",
            left: (e.pageX - 100) + "px"
        });
    });

    $("#locusPlotTwoUI").mousemove(function(e) {
        $("#locusPlotTwo_tooltip").show();         
        $("#locusPlotTwo_tooltip").css({
            top:  (e.pageY + 25) + "px",
            left: (e.pageX - 100) + "px"
        });
    });

    $("#profilePlotUI").mousemove(function(e) {
        $("#profilePlot_tooltip").show();         
        $("#profilePlot_tooltip").css({
            top:  (e.pageY + 25) + "px",
            left: (e.pageX - 100) + "px"
        });
    });
})

