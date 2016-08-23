# server.R
library(rgl)
library(shinyRGL)
library(shiny)
library(rglwidget)
source("helpers.R")

jet.colors <-   # function from grDevices package
colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                   "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
colorzjet <- jet.colors(100)  # 100 separate color 

shinyServer(function(input, output) {
    
    # plot P-wave radiation pattern
    vals_p <- eventReactive(input$plot, {
         rp(input$u, input$v, input$sigma,input$phi, input$delta,input$theta,input$eta)
    })
    output$plotRP <- renderWebGL({
        rad <- vals_p()
        minmin <- 1.1*min(c(rad$x, rad$y, rad$z))
        maxmax <- 1.1*max(c(rad$x, rad$y, rad$z))
        axes3d()
        surface3d(rad$x,rad$y,rad$z,
                  color = colorzjet[findInterval(rad$r, seq(min(rad$r), max(rad$r), length=100))])
        title3d('P-wave', input$date, 'North', 'East', 'Down')
        rgl.lines(c(minmin, maxmax), c(0, 0), c(0, 0), color = "black")
        rgl.lines(c(0, 0), c(minmin,maxmax), c(0, 0), color = "red")
        rgl.lines(c(0, 0), c(0, 0), c(minmin,maxmax), color = "green")
        
    })
    
    # plot S-wave radiation pattern
    vals_s <- eventReactive(input$plot, {
        rs(input$u, input$v, input$sigma,input$phi, input$delta,input$theta,input$eta)
    })
    output$plotRS <- renderWebGL({
        rad <- vals_s()
        minmin <- 1.1*min(c(rad$x, rad$y, rad$z))
        maxmax <- 1.1*max(c(rad$x, rad$y, rad$z))

        axes3d()
        surface3d(rad$x,rad$y,rad$z,
                  color = colorzjet[findInterval(rad$r, seq(min(rad$r), max(rad$r), length=100))])
        title3d('S-wave', input$date, 'North', 'East', 'Down')
        rgl.lines(c(minmin, maxmax), c(0, 0), c(0, 0), color = "black")
        rgl.lines(c(0, 0), c(minmin,maxmax), c(0, 0), color = "red")
        rgl.lines(c(0, 0), c(0, 0), c(minmin,maxmax), color = "green")
    })

    # plot SV-wave radiation pattern
    vals_sv <- eventReactive(input$plot, {
        rsv(input$u, input$v, input$sigma,input$phi, input$delta,input$theta,input$eta)
    })
    output$plotRSV <- renderWebGL({
        rad <- vals_sv()
        minmin <- 1.1*min(c(rad$x, rad$y, rad$z))
        maxmax <- 1.1*max(c(rad$x, rad$y, rad$z))
        bbox3d()
        surface3d(rad$x,rad$y,rad$z,
                  color = colorzjet[findInterval(rad$r, seq(min(rad$r), max(rad$r), length=100))])
        title3d('SV-wave', input$date, 'North', 'East', 'Down')
        rgl.lines(c(minmin, maxmax), c(0, 0), c(0, 0), color = "black")
        rgl.lines(c(0, 0), c(minmin,maxmax), c(0, 0), color = "red")
        rgl.lines(c(0, 0), c(0, 0), c(minmin,maxmax), color = "green")
    })

    # plot SH-wave radiation pattern
    vals_sh <- eventReactive(input$plot, {
        rsh(input$u, input$v, input$sigma,input$phi, input$delta,input$theta,input$eta)
    })
    output$plotRSH <- renderWebGL({
        rad <- vals_sh()
        minmin <- 1.1*min(c(rad$x, rad$y, rad$z))
        maxmax <- 1.1*max(c(rad$x, rad$y, rad$z))
        axes3d()
        surface3d(rad$x,rad$y,rad$z,
                  color = colorzjet[findInterval(rad$r, seq(min(rad$r), max(rad$r), length=100))])
        title3d('SH-wave', input$date, 'North', 'East', 'Down')
        rgl.lines(c(minmin, maxmax), c(0, 0), c(0, 0), color = "black")
        rgl.lines(c(0, 0), c(minmin,maxmax), c(0, 0), color = "red")
        rgl.lines(c(0, 0), c(0, 0), c(minmin,maxmax), color = "green")
    })

})