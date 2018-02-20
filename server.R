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

                                 show('hiddenThings')
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
        input$mgpCalc
        # vals$expression
        # input$factors
        # input$cellTypes
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
                
                if(input$abbreviate){
                    groups %<>% str_replace('\\|',' | ')  %>% abbreviate()
                }
                
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
        estimates()
        isolate({
            if(!is.null(vals$expression) & !is.null(input$factors) & !is.null(estimates())){
                print('make plot')
                toPlot = estimates()$estimates %>% melt
                names(toPlot) = c('mgp','cellType')
                # browser
                toPlot = data.frame(toPlot,groups = estimates()$groups[[1]])
                textSize = 11
                # if(input$abbreviate){
                #     toPlot$groups %<>% str_replace('\\|',' | ')  %>% abbreviate()
                #     textSize = 16
                # }
                p = toPlot %>% ggplot(aes(x = groups,y = mgp)) + 
                    facet_wrap(~cellType) + 
                    ogbox::geom_ogboxvio() + geom_jitter() + 
                    theme(axis.text.x = element_text(angle = 90, size = textSize))
                
                return(p)
                
            }
        })
    })
    
    output$qualityTable = renderDataTable({
        data.frame(removedMarkerRatio = estimates()$removedMarkerRatios,
                   varianceExplained = estimates()$trimmedPCAs %>% 
                       sapply(function(x){
                           x %>% summary %$% importance %>% {.[2,1]}
                       }),stringsAsFactors = FALSE) %>%datatable(selection = 'single')
    })
    
    output$groupInfo = renderDataTable({
        if(!is.null(estimates()) & !is.null(input$qualityTable_rows_selected)){
            groups = unique(estimates()$groups[[1]])
            cellType = input$qualityTable_rows_selected
            groupInfo = groups %>% lapply(function(x){
                group = estimates()$groups[[cellType]] %in% x
                meanEstimate = estimates()$estimates[[cellType]][group] %>% mean
                meanMarkerExpression = estimates()$meanUsedMarkerExpression[[cellType]][group] %>% mean
                meanScaledMarkerExpression = estimates()$simpleScaledEstimation[[cellType]][group] %>% mean
                return(c(meanEstimate = meanEstimate,
                         meanMarkerExpression = meanMarkerExpression,
                         meanScaledMarkerExpression = meanScaledMarkerExpression))
            }) %>% as.data.frame
            names(groupInfo) = groups
            return(groupInfo %>% t)
        }
    })
    
    output$groupComparison = renderPlot({
        if(!is.null(estimates()) & !is.null(input$qualityTable_rows_selected)){
            groups = unique(estimates()$groups[[1]])
            cellType = input$qualityTable_rows_selected
            pairwise = combn(groups,2)
            groupComparisons = seq_len(ncol(pairwise)) %>% sapply(function(k){
                
                pair = pairwise[,k]
                group1 = estimates()$groups[[cellType]] %in% pair[1]
                group2 = estimates()$groups[[cellType]] %in% pair[2]
                
                pVal = wilcox.test(estimates()$estimates[[cellType]][group1],
                                   estimates()$estimates[[cellType]][group2]) %$% p.value
                
            })
            sigMark = signifMarker(groupComparisons)
            
            toPlot = data.frame(g1 = pairwise[1,],
                                g2 = pairwise[2,],
                                pVal = groupComparisons,
                                sigMark = sigMark)
            
            ggplot(toPlot, aes(x = g1,y = g2)) + xlab('') + ylab('') +
                geom_tile(aes(fill = pVal))+ cowplot::theme_cowplot() +
                geom_text(aes(label = sigMark),size = 18) + scale_fill_viridis(limits = c(0,1),direction = -1)
            
        }
    })


})
