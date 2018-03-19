library(gemmaAPI)
library(shinyjs)
library(stringr)
library(ogbox)
library(memoise)
library(markerGeneProfile)
library(homologene)
library(reshape2)
library(ggplot2)
library(shiny)
library(DT)
library(magrittr)
library(viridis)
library(cowplot)
data("mouseMarkerGenes", package = 'markerGeneProfile')
data('mouseMarkerGenesNCBI',package = 'markerGeneProfile')

data("mouseMarkerGenesCombined", package = 'markerGeneProfile')
data('mouseMarkerGenesCombinedNCBI',package = 'markerGeneProfile')

data("mouseRegionHierarchy", package='markerGeneProfile')
mouseMarkerGenes = mouseMarkerGenesCombined
mouseMarkerGenesNCBI = mouseMarkerGenesCombinedNCBI



mem_median = memoise(median)

gemmaPrep = function(study){
    if(exists('session')){
        setProgress(0, detail ="Compiling metadata")
    }
    
    meta = compileMetadata(study,memoised = TRUE)
    species = meta$taxon %>% unique
    
    neededGenes = mouseMarkerGenesNCBI %>% unlist %>% unique
    neededGenes = neededGenes[!grepl('\\|',neededGenes)]
    taxonData = taxonInfo(species,memoised = TRUE)
    speciesID = taxonData[[1]]$ncbiId
    
    if(speciesID != 10090){
        neededGenes = homologene(neededGenes,inTax = 10090, outTax = speciesID)[[paste0(speciesID,'_ID')]]
    }
    
    # expression = datasetInfo(study,request= 'data',
    #                               IdColnames=TRUE,
    #                               memoised = TRUE)
    if(exists('session')){
        setProgress(4, detail ="Getting expression data")
    }

    # expression = datasetInfo(study,request= 'data',
    #                          IdColnames=TRUE,
    #                          memoised = TRUE)
    expression = gemmaAPI::expressionSubset(study,neededGenes,consolidate = 'pickvar')
    
    if(exists('session')){
        setProgress(7, detail ="Filtering expression data")
    }
    
    list[gene,exp] = expression %>% sepExpr()
    # some samples have outliers. they are NaNs in the entire row. remove them
    outlier = exp %>% apply(2,function(y){all(is.nan(y))})
    # some rows have NaNs remove them too
    keep = exp[!outlier] %>% apply(1,function(y){!any(is.nan(y))})
    gene = gene[keep,]
    exp = exp[!outlier][keep,]
    # note that there can be a difference between samples in metadata and
    # samples in filtered exp now
    exp = exp[,match(meta$sampleName,colnames(exp))]
    
    
    expression = cbind(gene,exp)
    expression = expression[!expression$GeneSymbol=='',]
    
    medExp = mem_median(exp %>% unlist)
    
    expression = mostVariable(expression,
                                   genes = "GeneSymbol",
                                   threshFun = max,
                                   threshold = medExp)
    
    
    
    
    return(list(meta,expression))
}

mem_gemmaPrep = memoise(gemmaPrep)



getCategoryAnnotations = function(data,
                                  category,
                                  categoryColumn,
                                  annotationColumns,
                                  split = '(?<=[^\\s])\\|(?=[^\\s])',
                                  merge=FALSE){
    categories = data[[categoryColumn]] %>% stringr::str_split(split)
    categoryMatch = categories %>% lapply(function(x){
        x %in% category
    })
    
    sepAnnots = data[annotationColumns] %>% apply(2,function(x){stringr::str_split(x,split)})
    relevantAnnots = lapply(seq_len(length(sepAnnots[[1]])),function(i){
        out = lapply(seq_along(sepAnnots), function(j){
            sepAnnots[[j]][[i]][categoryMatch[[i]]]
        })
        
        names(out) = annotationColumns
        return(out)
    })
    
    if(merge){
        relevantAnnots %<>% lapply(function(x){
            x %>% lapply(function(y){paste(y,collapse= '|')})
        })
    }
    
    return(relevantAnnots)
}

mem_mgpEstimate = memoise(mgpEstimate)