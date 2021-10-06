#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(

    # Application title
    titlePanel("Andersen Lab NIL Browser"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            shiny::radioButtons("chr_input", "Chrom:", choices = c("I", "II", "III", "IV", "V", "X"), selected = "V", inline = T),
            shiny::uiOutput("positions"),
            shiny::checkboxInput("eca_only", "ECA strains only?", value = FALSE), 
            downloadButton('downloadData', "Download data")
        ),

        # Show a plot of the generated distribution
        mainPanel(
            shiny::uiOutput("nil_output")
        )
    )
))
