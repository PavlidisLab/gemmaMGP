#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)

# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
    vals = reactiveValues(expression = NULL,
                          study = NULL,
                          query = NULL)
    
    
    observe({
        vals$querry = parseQueryString(session$clientData$url_search)
        print(vals$querry)
    })
    
    observe({
        if(input$submit >0){
            isolate({
                withProgress(message = 'Importing expression data',value = 0,
                             min=0,max = 10, expr = {
                                 vals$expression = datasetInfo(input$study,request= 'data',
                                                               IdColnames=TRUE,
                                                               memoised = TRUE)
                                 
                                 list[vals$metadata, vals$expression] = mem_gemmaPrep(input$study)
                                 
                                 setProgress(9, detail ="Compiling metadata")
                                 
                                 factors = 
                                     vals$metadata$sampleAnnotBroadCategory %>%
                                     str_split(pattern = '\\|') %>%
                                     unlist %>% unique
                                 
                                 updateCheckboxGroupInput(session,inputId = 'factors',
                                                          choices = factors,selected=factors[!factors %in% 'block'])
                                 show('factors')
                                 show('brainRegion')
                                 show('cellTypes')
                                 setProgress(10, detail ="Compiling metadata")
                                 
                             })
                
            })}
    })
    
    observe({
        updateCheckboxGroupInput(session,inputId= 'cellTypes',
                                 choices = names(mouseMarkerGenes[[input$brainRegion]]),
                                 selected= names(mouseMarkerGenes[[input$brainRegion]]))
    })
    
    output$mgpPlot = renderPlot({
        if(!is.null(vals$expression) & !is.null(input$factors)){
            markers = mouseMarkerGenes[[input$brainRegion]][input$cellTypes]
            
            species = vals$metadata$taxon %>% unique
            
            taxonData = taxonInfo(species,memoised = TRUE)
            speciesID = taxonData$`Homo sapiens`$ncbiId
            
            if(speciesID != 10090){
                geneTransform = function(x){
                    homologene(x,inTax = 10090, outTax  =speciesID)[[speciesID %>% as.character]]
                }
            } else {
                geneTransform = NULL
            }
            
            
            groups =
                getCategoryAnnotations(data = vals$metadata,
                                       category = input$factors,
                                       categoryColumn = 'sampleAnnotBroadCategory',
                                       annotationColumns = 'sampleAnnotation',
                                       split = '\\|',
                                       merge=FALSE) %>% sapply(function(x){
                                           x%>% unlist %>% sort %>% paste(collapse='|')
                                       })
            
            estimates = mgpEstimate(exprData = vals$expression,
                                    genes = markers,
                                    geneColName = "GeneSymbol",
                                    geneTransform = geneTransform,
                                    groups = groups)
            
            # estimates %<>% lapply(function(x){
            #     x = x[!names(x) %in% c('Microglia_activation',
            #                           'Microglia_deactivation')]
            # })
            
            estimates$estimates %<>% lapply(scale01)
            
            toPlot = estimates$estimates %>% melt
            names(toPlot) = c('mgp','cellType')
            # browser
            toPlot = data.frame(toPlot,groups = estimates$groups[[1]])
            p = toPlot %>% ggplot(aes(x = groups,y = mgp)) + 
                facet_wrap(~cellType) + 
                ogbox::geom_ogboxvio() + 
                theme(axis.text.x = element_text(angle = 90, size = 10))
            
            return(p)
    
        }
    })


})
