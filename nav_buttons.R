SHIFT_1_FRACTION = 0.1
SHIFT_2_FRACTION = 0.475
SHIFT_3_FRACTION = 0.95

ZOOM_1_FRACTION = 1.5
ZOOM_2_FRACTION = 3
ZOOM_3_FRACTION = 10

observeEvent(input$MainPlotLeftButton1 ,{
    MainPlotValues$window = shiftLeft(MainPlotValues$window, SHIFT_1_FRACTION)
})

observeEvent(input$MainPlotLeftButton2 ,{
    MainPlotValues$window = shiftLeft(MainPlotValues$window, SHIFT_2_FRACTION)
})

observeEvent(input$MainPlotLeftButton3 ,{
    MainPlotValues$window = shiftLeft(MainPlotValues$window, SHIFT_3_FRACTION)
})

observeEvent(input$MainPlotRightButton1,{
    writeLog('MainPlot.RightButton1')
    MainPlotValues$window = shiftRight(MainPlotValues$window, SHIFT_1_FRACTION)
})

observeEvent(input$MainPlotRightButton2,{
    MainPlotValues$window = shiftRight(MainPlotValues$window, SHIFT_2_FRACTION)
})

observeEvent(input$MainPlotRightButton3,{
    MainPlotValues$window = shiftRight(MainPlotValues$window, SHIFT_3_FRACTION)
})

observeEvent(input$MainPlotZoomInButton1,{
    MainPlotValues$window = zoomIn(MainPlotValues$window, ZOOM_1_FRACTION)
})

observeEvent(input$MainPlotZoomInButton2,{
    MainPlotValues$window = zoomIn(MainPlotValues$window, ZOOM_2_FRACTION)
})

observeEvent(input$MainPlotZoomInButton3,{
    MainPlotValues$window = zoomIn(MainPlotValues$window, ZOOM_3_FRACTION)
})

observeEvent(input$MainPlotZoomOutButton1,{
    MainPlotValues$window = zoomOut(MainPlotValues$window, ZOOM_1_FRACTION)
})

observeEvent(input$MainPlotZoomOutButton2,{
    MainPlotValues$window = zoomOut(MainPlotValues$window, ZOOM_2_FRACTION)
})

observeEvent(input$MainPlotZoomOutButton3,{
    MainPlotValues$window = zoomOut(MainPlotValues$window, ZOOM_3_FRACTION)
})

observeEvent(input$locusPlotOneLeftButton1 ,{
    PlotOneValues$window = shiftLeft(PlotOneValues$window, SHIFT_1_FRACTION)
})

observeEvent(input$locusPlotOneLeftButton2 ,{
    PlotOneValues$window = shiftLeft(PlotOneValues$window, SHIFT_2_FRACTION)
})

observeEvent(input$locusPlotOneLeftButton3 ,{
    PlotOneValues$window = shiftLeft(PlotOneValues$window, SHIFT_3_FRACTION)
})

observeEvent(input$locusPlotOneRightButton1,{
    writeLog('locusPlotOne.RightButton1')
    PlotOneValues$window = shiftRight(PlotOneValues$window, SHIFT_1_FRACTION)
})

observeEvent(input$locusPlotOneRightButton2,{
    PlotOneValues$window = shiftRight(PlotOneValues$window, SHIFT_2_FRACTION)
})

observeEvent(input$locusPlotOneRightButton3,{
    PlotOneValues$window = shiftRight(PlotOneValues$window, SHIFT_3_FRACTION)
})

observeEvent(input$locusPlotOneZoomInButton1,{
    PlotOneValues$window = zoomIn(PlotOneValues$window, ZOOM_1_FRACTION)
})

observeEvent(input$locusPlotOneZoomInButton2,{
    PlotOneValues$window = zoomIn(PlotOneValues$window, ZOOM_2_FRACTION)
})

observeEvent(input$locusPlotOneZoomInButton3,{
    PlotOneValues$window = zoomIn(PlotOneValues$window, ZOOM_3_FRACTION)
})

observeEvent(input$locusPlotOneZoomOutButton1,{
    PlotOneValues$window = zoomOut(PlotOneValues$window, ZOOM_1_FRACTION)
})

observeEvent(input$locusPlotOneZoomOutButton2,{
    PlotOneValues$window = zoomOut(PlotOneValues$window, ZOOM_2_FRACTION)
})

observeEvent(input$locusPlotOneZoomOutButton3,{
    PlotOneValues$window = zoomOut(PlotOneValues$window, ZOOM_3_FRACTION)
})

observeEvent(input$locusPlotTwoLeftButton1 ,{
    PlotTwoValues$window = shiftLeft(PlotTwoValues$window, SHIFT_1_FRACTION)
})

observeEvent(input$locusPlotTwoLeftButton2 ,{
    PlotTwoValues$window = shiftLeft(PlotTwoValues$window, SHIFT_2_FRACTION)
})

observeEvent(input$locusPlotTwoLeftButton3 ,{
    PlotTwoValues$window = shiftLeft(PlotTwoValues$window, SHIFT_3_FRACTION)
})

observeEvent(input$locusPlotTwoRightButton1,{
    writeLog('locusPlotTwo.RightButton1')
    PlotTwoValues$window = shiftRight(PlotTwoValues$window, SHIFT_1_FRACTION)
})

observeEvent(input$locusPlotTwoRightButton2,{
    PlotTwoValues$window = shiftRight(PlotTwoValues$window, SHIFT_2_FRACTION)
})

observeEvent(input$locusPlotTwoRightButton3,{
    PlotTwoValues$window = shiftRight(PlotTwoValues$window, SHIFT_3_FRACTION)
})

observeEvent(input$locusPlotTwoZoomInButton1,{
    PlotTwoValues$window = zoomIn(PlotTwoValues$window, ZOOM_1_FRACTION)
})

observeEvent(input$locusPlotTwoZoomInButton2,{
    PlotTwoValues$window = zoomIn(PlotTwoValues$window, ZOOM_2_FRACTION)
})

observeEvent(input$locusPlotTwoZoomInButton3,{
    PlotTwoValues$window = zoomIn(PlotTwoValues$window, ZOOM_3_FRACTION)
})

observeEvent(input$locusPlotTwoZoomOutButton1,{
    PlotTwoValues$window = zoomOut(PlotTwoValues$window, ZOOM_1_FRACTION)
})

observeEvent(input$locusPlotTwoZoomOutButton2,{
    PlotTwoValues$window = zoomOut(PlotTwoValues$window, ZOOM_2_FRACTION)
})

observeEvent(input$locusPlotTwoZoomOutButton3,{
    PlotTwoValues$window = zoomOut(PlotTwoValues$window, ZOOM_3_FRACTION)
})

observeEvent(input$ProfilePlotLeftButton1 ,{
    ProfilePlotValues$window = shiftLeft(ProfilePlotValues$window, SHIFT_1_FRACTION)
})

observeEvent(input$ProfilePlotLeftButton2 ,{
    ProfilePlotValues$window = shiftLeft(ProfilePlotValues$window, SHIFT_2_FRACTION)
})

observeEvent(input$ProfilePlotLeftButton3 ,{
    ProfilePlotValues$window = shiftLeft(ProfilePlotValues$window, SHIFT_3_FRACTION)
})

observeEvent(input$ProfilePlotRightButton1,{
    writeLog('ProfilePlot.RightButton1')
    ProfilePlotValues$window = shiftRight(ProfilePlotValues$window, SHIFT_1_FRACTION)
})

observeEvent(input$ProfilePlotRightButton2,{
    ProfilePlotValues$window = shiftRight(ProfilePlotValues$window, SHIFT_2_FRACTION)
})

observeEvent(input$ProfilePlotRightButton3,{
    ProfilePlotValues$window = shiftRight(ProfilePlotValues$window, SHIFT_3_FRACTION)
})

observeEvent(input$ProfilePlotZoomInButton1,{
    ProfilePlotValues$window = zoomIn(ProfilePlotValues$window, ZOOM_1_FRACTION)
})

observeEvent(input$ProfilePlotZoomInButton2,{
    ProfilePlotValues$window = zoomIn(ProfilePlotValues$window, ZOOM_2_FRACTION)
})

observeEvent(input$ProfilePlotZoomInButton3,{
    ProfilePlotValues$window = zoomIn(ProfilePlotValues$window, ZOOM_3_FRACTION)
})

observeEvent(input$ProfilePlotZoomOutButton1,{
    ProfilePlotValues$window = zoomOut(ProfilePlotValues$window, ZOOM_1_FRACTION)
})

observeEvent(input$ProfilePlotZoomOutButton2,{
    ProfilePlotValues$window = zoomOut(ProfilePlotValues$window, ZOOM_2_FRACTION)
})

observeEvent(input$ProfilePlotZoomOutButton3,{
    ProfilePlotValues$window = zoomOut(ProfilePlotValues$window, ZOOM_3_FRACTION)
})


