#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#


# Define UI for application that draws a histogram
shinyUI(fluidPage(
    useShinyjs(),
    titlePanel("GemmaMGP"),
    
    sidebarLayout(
        sidebarPanel(
            textInput(inputId = 'study',label ="Study ID/Name",value = 'GSE7621'),
            actionButton(inputId = 'submit',label = 'Submit'),
            hidden(
                div(id = 'hiddenThings',
                    checkboxInput(inputId = 'abbreviate',label = 'Abbreviate Factors?', value = FALSE),
                    checkboxGroupInput(inputId = 'factors',label = 'Pick factors'),
                    selectInput(inputId = 'brainRegion',
                                label = 'Pick brain region',
                                choices = mouseRegionHierarchy %>% 
                                    unlist %>%
                                    names %>%
                                    str_split('\\.') %>%
                                    unlist %>% 
                                    unique,
                                selected = 'Cortex'),
                    checkboxGroupInput(inputId = 'cellTypes',
                                       label = 'Cell Types',
                                       choices = names(mouseMarkerGenes$Cortex),
                                       selected = names(mouseMarkerGenes$Cortex)),
                    actionButton(inputId='mgpCalc',label = 'Calculate MGPs')
                )
            )
        ),
        
        # Show a plot of the generated distribution
        mainPanel(
            plotOutput('mgpPlot'),
            dataTableOutput('qualityTable'),
            dataTableOutput('groupInfo'),
            plotOutput('groupComparison')
        )
    )
))


