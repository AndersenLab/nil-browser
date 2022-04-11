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
            shiny::radioButtons("select_option", "Do you want to:", choices = c("Explore NILs", "Plot specific NILs"), selected = "Explore NILs", inline = T),
            shiny::uiOutput("choose_type")
        ),

        # Show a plot of the generated distribution
        mainPanel(
            shiny::uiOutput("nil_output")
        )
    )
))
