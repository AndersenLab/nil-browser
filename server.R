#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)
library(DT)

# source functions
source("functions.R")

# Define server logic required to draw a histogram
shinyServer(function(input, output) {

    # read in data
    nilgeno <- readr::read_tsv("https://github.com/AndersenLab/andersenlab.github.io/raw/kse_update_nil_data/pages/gt_hmm.tsv") %>%
        dplyr::mutate(nil_genotype = ifelse(gt == 1, "N2", "CB4856"))
    
    # initialize reactive value for length of plot
    rv <- shiny::reactiveValues(h = "1000px")
    rv <- shiny::reactiveValues(pos = 21e6)
    
    # get nil data
    getNILdata <- shiny::reactive({

        # filter strains - ECA only
        if(input$eca_only) {
            test <- nilgeno %>%
                dplyr::filter(grepl("ECA", strain))
        } else {
            test <- nilgeno
        }
        
        right_pos <- input$end_pos*1e6
        left_pos <- input$start_pos*1e6
        
        # filter strains - breakpoint on chromosome
        # Filter nil genotype data:
        # - filter to chromosome of interest
        # - filter to keep only ECA strains
        # - filter the length of a segment to be > 10,000 bp
        # - filter to keep only strains with N2 and CB on chrom
        nils <- test  %>%
            dplyr::filter(!(strain %in% c("N2", "CB4856"))) %>%
            dplyr::mutate(length = as.numeric(end) - start) %>%
            dplyr::group_by(strain, gt) %>%
            dplyr::mutate(background = sum(length)) %>%
            tidyr::pivot_wider(names_from = nil_genotype, values_from = background) %>%
            dplyr::group_by(strain) %>%
            dplyr::mutate(CB4856 = mean(CB4856, na.rm = T),
                          N2 = mean(N2, na.rm = T),
                          RIAIL = N2/CB4856) %>%
            dplyr::filter(RIAIL > 5 | RIAIL < 0.25) %>% #0.167 is one chromosome (CSS) but better safer than sorry probably.
            dplyr::group_by(strain) %>%
            dplyr::mutate(backgroundtype = ifelse(N2 > CB4856, "N2", ifelse(CB4856 > N2, "CB4856", NA))) %>%
            dplyr::filter(chromosome == input$chr_input) %>%
            # dplyr::filter(chromosome == "V") %>%
            dplyr::filter(length > 100000) %>% # what should I use as the threshold? smallest NIL in chrIII nils is 241,230 bp
            dplyr::group_by(strain) %>%
            dplyr::mutate(num_genos = sum(plyr:::nunique(gt))) %>%
            dplyr::filter(num_genos == 2, strain != "ECA257") # this is the kammenga lab N2 
        
        # Separate into N2 and CB nils, then make sure the NIL region is within the boundary
        cbnils <- nils %>%
            dplyr::filter(backgroundtype == "N2", gt == 2, start < right_pos, end > left_pos) #1 = N2, 2 = CB
        
        n2nils <- nils %>%
            dplyr::filter(backgroundtype == "CB4856", gt == 1, start < right_pos, end > left_pos)
        
        allnils <- unique(rbind(cbnils, n2nils)$strain)
        
        rv$h <- ifelse(length(allnils) > 20, glue::glue("{length(allnils)}0px"), "200px")
        rv$pos <- test %>% dplyr::filter(chromosome == input$chr_input) %>% dplyr::pull(end) %>% max()
        
        # plot
        nil_plot(test %>% dplyr::filter(strain %in% allnils), input$chr_input, background = T, ci = c(left_pos/1e6, right_pos/1e6))
        
    })
    
    # plot nils
    output$nilplot <- shiny::renderPlot({
        
        nil_data <- getNILdata()
        nil_data[[1]]
        
    })
    
    # ui for nil output
    output$nil_output <- shiny::renderUI({
        
        # show datatable of genotypes
        output$niltable <- DT::renderDT(
            getNILdata()[[4]] %>%
                dplyr::select(chromosome, start, end, strain, nil_genotype), filter = list(
                    position = 'top', clear = FALSE)
        )
        
        tagList(
            
            # show output for plot
            shiny::plotOutput("nilplot", height = rv$h),
            DT::DTOutput("niltable")
        )
        
    })
    
    output$positions <- shiny::renderUI({
        tagList(
            shiny::numericInput("start_pos", "Left Pos (Mb):", value = 0, min = 0, max = 20),
            shiny::numericInput("end_pos", "Right Pos (Mb):", value = rv$pos/1e6, min = 1, max = 21)
        )
    })
    
    # handle download of dataset for nil geno
    output$downloadData <- shiny::downloadHandler(
        filename = "nil_genotypes.csv",
        content = function(file) {
            write.csv(getNILdata()[[4]], file, row.names = FALSE)
        }
    )
    
    # output$saveButton <- shiny::renderUI({
    #     req(input$file1)
    #     downloadButton('saveImage', 'Save plot')
    # })
    # 
    # output$saveImage <- shiny::downloadHandler(
    #     filename = function() { 'NILs.png' },
    #     content = function(file) {
    #         ggsave(file, plot = getNILdata()[[1]], device = "png", height = 10, width = 7) # need reactive height again.
    #     }
    # )
    
    
})


# nil_plot(test, "III", background = T, ci = c(5, 10))[[1]]


