#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#


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

                                 list[vals$metadata, vals$expression] = mem_gemmaPrep(input$study)
                                 
                                 setProgress(9, detail ="Setting up UI")
                                 
                                 factors = 
                                     vals$metadata$sampleAnnotBroadCategory %>%
                                     str_split(pattern = '\\|') %>%
                                     unlist %>% unique
                                 
                                 updateCheckboxGroupInput(session,inputId = 'factors',
                                                          choices = factors,selected=factors[!factors %in% 'block'])
                                 show('factors')
                                 show('brainRegion')
                                 show('cellTypes')
                                 show('abbreviate')
                                 setProgress(10, detail ="DONE!")
                                 
                             })
                
            })}
    })
    
    observe({
        updateCheckboxGroupInput(session,inputId= 'cellTypes',
                                 choices = names(mouseMarkerGenes[[input$brainRegion]]),
                                 selected= names(mouseMarkerGenes[[input$brainRegion]]))
        print('update checkbox')
    })
    
    estimates = reactive({
        vals$expression
        input$factors
        input$cellTypes
        isolate({
            if(!is.null(vals$expression) & !is.null(input$factors)){
                markers = mouseMarkerGenes[[input$brainRegion]][input$cellTypes]
                
                species = vals$metadata$taxon %>% unique
                
                taxonData = taxonInfo(species,memoised = TRUE)
                speciesID = taxonData[[1]]$ncbiId
                
                if(speciesID != 10090){
                    geneTransform = function(x){
                        homologene(x,inTax = 10090, outTax = speciesID)[[speciesID %>% as.character]]
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
                
                estimates = mem_mgpEstimate(exprData = vals$expression,
                                            genes = markers,
                                            geneColName = "GeneSymbol",
                                            geneTransform = geneTransform,
                                            groups = groups)
                estimates$estimates %<>% lapply(scale01)
                
                
                keep = estimates$estimates %>% is.na() %>% not
                
                estimates %<>% lapply(function(x){
                    x[keep]
                })
                
                return(estimates)
            } else {
                return(NULL)
            }
        })
    })
    
    output$mgpPlot = renderPlot({
        if(!is.null(vals$expression) & !is.null(input$factors)){
            print('make plot')
            toPlot = estimates()$estimates %>% melt
            names(toPlot) = c('mgp','cellType')
            # browser
            toPlot = data.frame(toPlot,groups = estimates()$groups[[1]])
            textSize = 11
            if(input$abbreviate){
                toPlot$groups %<>% str_replace('\\|',' | ')  %>% abbreviate()
                textSize = 16
            }
            p = toPlot %>% ggplot(aes(x = groups,y = mgp)) + 
                facet_wrap(~cellType) + 
                ogbox::geom_ogboxvio() + geom_jitter() + 
                theme(axis.text.x = element_text(angle = 90, size = textSize))
            
            return(p)
    
        }
    })
    
    output$qualityTable = renderDataTable({
        data.frame(removedMarkerRatio = estimates()$removedMarkerRatios,
                   varianceExplained = estimates()$trimmedPCAs %>% 
                       sapply(function(x){
                           x %>% summary %$% importance %>% {.[2,1]}
                       }),stringsAsFactors = FALSE) %>% 
            DT::datatable(rownames = TRUE)
    })


})
