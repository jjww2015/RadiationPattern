library(rgl)
library(shinyRGL)
library(shiny)
library(rglwidget)

shinyUI(fluidPage(
  titlePanel("Earthquake Energy Radiation Pattern Visualization"),
  fluidRow(
      column(width =12,
             column(3,
                    h4("Define earthquake source to visualize its radiation pattern."),
                    h5("Click Plot to update"),
                    helpText("Hover over the text to see more help info. See the cartoon on the right for the geometry"),
                    wellPanel(
                            tags$div(title="controls the energy released by shear source:",
                                     numericInput("u", label = h5(withMathJax(strong("Fault disp \\(*\\) area, \\(m^3\\)"))), value = 1, min = 0)),

                            tags$div(title="controls the energy released by volume source:",
                                     numericInput("v", label = h5(withMathJax(strong("Volume explosion, \\(m^3\\)"))), value = 0, min = 0 )),
                            
                            tags$div(title="only matters when tensile opening angle is non-zero:",
                                     sliderInput("sigma", label = h5(withMathJax(strong("Poisson ratio"))), min = 0,
                                                 max = 0.49999, value = 0.25)),

                            tags$div(title="the fault orientation clockwise measured from North",
                                     sliderInput("phi", label = h5(withMathJax(strong("Strike (\\(^o\\))"))), min = 0,
                                        max = 360, value = 0)),
                            tags$div(title="the fault orientation positively downward from horizontal plane",
                                     sliderInput("delta", label = h5(withMathJax(strong("Dip (\\(^o\\))"))), min = -90,
                                                 max = 90, value = 45)),
                            tags$div(title="the tensile angle (0~90) deg between two fault planes",
                                     sliderInput("theta", label = h5(withMathJax(strong("Open angle (\\(^o\\))"))), min = 0,
                                                 max = 90, value = 0)),
                            tags$div(title="the fault moving angle relative to the strike.
                                     It takes 180 deg along the strike and -180 deg in the opposite direction of the strike",
                                     sliderInput("eta", label = h5(withMathJax(strong("Rake (\\(^o\\))"))), min = -180,
                                                 max = 180, value = 0)),

                            #tags$div(title="The name you want to put on your figure",
                                     #textInput("name", "Author's Name", "First and Family Name")),

                            tags$div(title="The date stamp you want to put on your figure",
                                     dateInput("date", label = h5(strong("Date stamp")), value = as.character(Sys.Date())))
                    ),
                   actionButton("plot", label = "Plot", width = '100%')
             ),
             column(9,
                    h4("The cartoon below illustrates the angle definitions."),
                    img(src='faultcartoon.png', align = "center", height = 400),
                    br(),
                    h5("If interactive figures are not shown below after clicking 'plot', 
                       please wait a few seconds and click your mouse in the blank area below. 
                       It is a known bug from 'shiny-rbl' package that your browser 
                       would complain about javascript. Please ignore it."),
                    webGLOutput("plotRP"),
                    webGLOutput("plotRS"),
                    webGLOutput("plotRSV"),
                    webGLOutput("plotRSH")
                    
             )
      )
  )

))