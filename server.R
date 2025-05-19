# A shiny app for RNAseq data analysis


shinyServer(function(input, output, session) {

  
############Panel 1: Asscoations between RNAseq and preteomics########
  
  #A reactive object for preprocessing proteomics data
  protNorm <- reactive({
    protCLL <- readRDS(paste0("data/",input$seleSet,".rds"))
    
    
    #filter out proteins on sex chromsomes
    if (input$ifNoSex) {
      protCLL <- protCLL[!rowData(protCLL)$chromosome_name %in% c("X","Y"),]
    }
    
    protCLL.new <- protCLL
    assays(protCLL.new)[assayNames(protCLL.new)] <- NULL
    
    #imputation (if specified)
    if (input$ifImpute == "yes") {
        if(input$seleSet == "LUMOS") {
          assays(protCLL.new)[["count"]] <- assays(protCLL)[["QRILC_combat"]]
          assays(protCLL.new)[["count_raw"]] <- assays(protCLL)[["QRILC"]]
        } else {
          assays(protCLL.new)[["count"]] <- assays(protCLL)[["QRILC"]]
          assays(protCLL.new)[["count_raw"]] <- assays(protCLL)[["QRILC"]]
        }
    } else {
      if(input$seleSet == "LUMOS") {
        assays(protCLL.new)[["count"]] <- assays(protCLL)[["count_combat"]]
        assays(protCLL.new)[["count_raw"]] <- assays(protCLL)[["count"]]
      } else {
        assays(protCLL.new)[["count"]] <- assays(protCLL)[["count"]]
        assays(protCLL.new)[["count_raw"]] <- assays(protCLL)[["count"]]
      }
    }
    
    protCLL.new
 })
  
#  #text output to show how many samples dataset
  output$numAllSample <- renderText({
    paste0("Number of samples: ", ncol(protNorm()))
  })
  
  #text output to show how many proteins remain after filtering
  output$numProt <- renderText({
    paste0("Number of proteins: ", nrow(protNorm()))
  })
  

  #A reactive object to calculate correlations between RNAseq and proteomics
  rnaCorr <- reactive({
    #subset data
    protSub <- protNorm()
    rnaMat <- assay(ddsSub)
    overSample <- intersect(colnames(protSub),colnames(rnaMat))
    protSub <- protSub[!is.na(rowData(protSub)$ensembl_gene_id),overSample]
    rnaMat <- assay(ddsSub)[,overSample]
    protMat <- assay(protSub)
    rownames(protMat) <- rowData(protSub)$ensembl_gene_id
    overGene <- intersect(rownames(protMat), rownames(rnaMat))
    protMat <- protMat[overGene,]
    rnaMat <- rnaMat[overGene, ]
    

    #calculate associations using pearson correlation
    resTab <- data.frame(id = rownames(rnaMat), 
                         coef = rep(NA,nrow(rnaMat)),
                         P = rep(NA,nrow(rnaMat)),
                         stringsAsFactors = FALSE)
    
    withProgress(message = "Loading results, please wait...", value = NULL, {
      for (i in seq(nrow(rnaMat))) {
        res <- cor.test(rnaMat[i,],protMat[i,],use = "pairwise.complete.obs")
        resTab[i,2] <- res$estimate
        resTab[i,3] <- res$p.value
      }
      })
    
   #annotate
    resTab <- filter(resTab, !is.na(P), !is.na(coef)) %>%
      mutate(adj.P = p.adjust(P, method = "BH"),
                     symbol = rowData(ddsSub[id,])$symbol,
             chromosome  = rowData(ddsSub[id,])$chromosome)
    
    list(resTab=resTab, rnaMat = rnaMat, protMat = protMat)
  })
  
  # Filter the correlation data based on selection
  rnaCorrFilt <- reactive({
    resTab <- rnaCorr()$resTab
    coefRange <- as.numeric(input$coefRange)
    resTab.filt <- filter(resTab, coef >= coefRange[1] & coef <= coefRange[2],
                          adj.P <= as.numeric(input$fdrRange)) %>% 
      arrange(desc(abs(coef)))
    resTab.filt
  })
  
  # Table output to show correlations
  output$rnaTab <- DT::renderDataTable({
    if (!is.null(rnaCorrFilt())) {
      rnaCorrFilt() %>% 
        select(symbol, chromosome, P, adj.P, coef) %>%
        dplyr::rename(`P-value` = P, `adjusted P-value` = adj.P, coefficient = coef) %>%
        mutate_if(is.numeric, formatC, digits=2) %>%
      datatable(selection = 'single', rownames = FALSE, caption = "RNA-protein correlations (click row to show correlation plot)", 
                options = list(
                  paging =TRUE,
                  pageLength =  5 
                )) 
    }
  })
  
  #output of coefficient histogram
  output$coefPlot <- renderPlot({
    plotTab <- rnaCorr()$resTab
    coefRange <- as.numeric(input$coefRange)
    recTab <- data.frame(x1 = coefRange[1],x2=coefRange[2],y1=0, y2=Inf)
    ggplot() + 
      geom_histogram(data = plotTab, aes(x=coef), col = "grey50", bins =50, fill = "lightblue") +
      geom_rect(data = recTab, aes(xmin = x1, xmax=x2, ymin = y1, ymax =y2), fill  = "red", alpha =0.2) +
      xlab("Pearson correlation coefficient") + ggtitle("")+
      theme_bw() + theme(plot.title = element_text(hjust=0.5))
  })
  
  #Correlation plot
  output$rnaPlot <- renderPlotly({
    lastClicked <- input$rnaTab_row_last_clicked
    if (is.null(lastClicked)) lastClicked <- 1
    id <- rnaCorrFilt()[lastClicked,]$id
    symbol <- rnaCorrFilt()[lastClicked,]$symbol
    plotTab <- tibble(rna = rnaCorr()$rnaMat[id,],
                      prot = rnaCorr()$protMat[id,],
                      patID = colnames(rnaCorr()$rnaMat)) %>%
      filter(!is.na(prot))
    p <- ggplot(plotTab, aes(x=rna,y=prot,label = patID)) + geom_point() + geom_smooth(method = "lm") +
      ggtitle(symbol) + xlab("RNA expression") + ylab("Protein expression") +
      theme_bw() + theme(plot.title = element_text(face = "bold", hjust =0.5))
    ggplotly(p) %>% config(displayModeBar = F)
  })
  
############Panel 2: Association test#####################
  
  ####Widgets####
  
  #the selection box for choosing a binary feature
  output$sampleSubsetBox <- renderUI({
    availList <- colnames(geneMatAll())
    selectInput("seleSubsetFeature","Subset samples by", availList)
  })
  
  output$statusSubsetBox <- renderUI({
    availStatus <- sort(unique(geneMatAll()[[input$seleSubsetFeature]]))
    selectInput("seleSubsetStatus", "Sample sub-group to include", availStatus)
  })
  
  #the selection box for choosing a binary feature
  output$featureBox <- renderUI({
    availList <- colnames(geneMat())
    availList <- availList[!availList %in% input$seleSubsetFeature]
    selectInput("seleFeature","Select a feature for testing", availList)
  })
  

  
  #the selection box for choosing blocking in the model
  output$blockingBox <- renderUI({
    ighvTri12Mat <- geneMat()[,colnames(geneMat()) %in% c("IGHV.status","trisomy12"), drop=FALSE]
    keepCol <- apply(ighvTri12Mat, 2, function(x) length(unique(x[!is.na(x)]))>1)
    blockList <- colnames(ighvTri12Mat)[keepCol]
    blockList <- blockList[!blockList %in% input$seleFeature]
    if (length(blockList >0)) {
      checkboxGroupInput("seleBlocking", "Blocking for", blockList, NULL)
    }
  })
  
  
  #a button to download DE gene table or enrichment results as csv
  output$downloadTable <- downloadHandler(
    filename = function() { paste('CorTable', '.csv', sep='') },
    content = function(file) {
      write.csv2(filterDE(), file)
    }
  )
  
  ###Reactive Objects####
  
  # a reactive value to monitor if plot histogram or boxplot
  ifHistogram <- reactiveValues(value = TRUE)
  
  #if table is clicked, turn the ifHistogram to false
  observeEvent(input$DEtab_row_last_clicked,{
    ifHistogram$value <- FALSE
  })
  
  
  #subset samples by user selection
  protSub <- reactive({
    if (input$subsetBox) {
      subPat <- intersect(patBack$Patient.ID[patBack[[input$seleSubsetFeature]] %in% input$seleSubsetStatus], colnames(protNorm()))
      protNorm()[,colnames(protNorm()) %in% subPat]
    } else {
      protNorm()
    }
  })
  
  #a reactive object to get genetic background matrix of the selected dataset
  geneMatAll <- reactive({
    geneBack <- patBack[match(colnames(protNorm()), patBack$Patient.ID ),]
    geneBack <- geneBack[,!colnames(geneBack ) %in% c("Patient.ID", "project", "diagnosis")]
    geneBack <- mutate_all(geneBack, as.factor) #change all to factors
    geneBack <- data.frame(geneBack)
    rownames(geneBack) <- colnames(protNorm())
    #only use features with enough cases in each level 
    keep <- apply(geneBack, 2, function(x) {
      x <- factor(x)
      l <- levels(x)
      length(l) == 2 & (sum(x == l[1], na.rm = TRUE) >=3 & sum(x == l[2], na.rm = TRUE) >=3)
    })
    geneBack <- geneBack[,keep,drop=FALSE]
    
    #reorder by mutated cases
    mutSum <- apply(geneBack, 2, function(x) sum(x!=0, na.rm = TRUE))
    geneBack <- geneBack[,order(mutSum, decreasing = TRUE), drop = FALSE]
    geneBack
  })
  
  #a reactive object to get genetic background matrix on the subseted sample
  geneMat <- reactive({
    geneBack <- patBack[match(colnames(protSub()), patBack$Patient.ID),]
    geneBack <- geneBack[, ! colnames(geneBack) %in% c("Patient.ID", "project", "diagnosis")]
    geneBack <- mutate_all(geneBack, as.factor) #change all to factors
    geneBack <- data.frame(geneBack)
    rownames(geneBack) <- colnames(protSub())
    #only use features with enough cases in each level 
    keep <- apply(geneBack, 2, function(x) {
      x <- factor(x)
      l <- levels(x)
      length(l) == 2 & (sum(x == l[1], na.rm = TRUE) >=3 & sum(x == l[2], na.rm = TRUE) >=3)
    })
    
    geneBack <- geneBack[,keep,drop=FALSE]
    
    #reorder by mutated cases
    mutSum <- apply(geneBack, 2, function(x) sum(x!=0, na.rm = TRUE))
    geneBack <- geneBack[,order(mutSum, decreasing = TRUE), drop = FALSE]
    geneBack
    
  })
  

  corrRes <- reactiveVal(value = NULL)
  
  #reactive event for calculating differential expression
  observeEvent(input$calcCorr, {
    withProgress(message = "Running differential expression analysis, please wait...", value = NULL, {
      #prepare design 
      dMat <- data.frame(row.names = colnames(protSub()), batch = factor(protSub()$batch))
      dMat <- cbind(dMat,geneMat()[,c(input$seleBlocking, input$seleFeature),drop=FALSE])
      dMat <- model.matrix(~., dMat)
      coefName <- colnames(dMat)[ncol(dMat)]
      
      
      #prepare input protein expression matrix
      protMat <- assays(protSub())[["count_raw"]]
      protMat <- protMat[,rownames(dMat)]
      
      if (input$seleTestModel == "Limma") {
        #moderate test using limma
        fit <- lmFit(protMat,  dMat)
        fit2 <- eBayes(fit)
        
        corRes <- topTable(fit2, number ="all", coef = coefName, adjust.method = "BH") %>% rownames_to_column("id") %>%
          mutate(symbol = rowData(protNorm()[id,])$hgnc_symbol,
                 chr = rowData(protNorm()[id,])$chromosome_name) %>%
          select(id, symbol, chr, logFC, P.Value, adj.P.Val, t) %>%
          arrange(P.Value)
        
      } else if (input$seleTestModel == "proDA") {
        fit <- proDA(protMat, design = dMat)
        corRes <- test_diff(fit, coefName) %>%
          dplyr::rename(id = name, logFC = diff, t=t_statistic,
                        P.Value = pval, adj.P.Val = adj_pval) %>% 
          mutate(symbol = rowData(protNorm()[id,])$hgnc_symbol,
                 chr = rowData(protNorm()[id,])$chromosome_name) %>%
          select(id, symbol, chr, logFC, P.Value, adj.P.Val, t) %>%
          arrange(P.Value)
      }
      ifHistogram$value <- TRUE
      corrRes(corRes)
    })
  })
  
 filterDE <- reactive({
   if (!is.null(corrRes())) {
     DEtab <- corrRes()
     if(input$ifAdjusted) {
        DEtab <- filter(DEtab, abs(logFC) >= input$fcFilter, adj.P.Val <= as.numeric(input$pFilter))
     } else {
        DEtab <- filter(DEtab, abs(logFC) >= input$fcFilter, P.Value <= as.numeric(input$pFilter))
     }
     
     DEtab
   }
 })
  
 ####Output#####
 
  #table for differentially expressed genes
  output$DEtab <- DT::renderDataTable({
    if (!is.null(corrRes())) {
      showTab <- dplyr::select(filterDE(), symbol, chr, logFC, P.Value, adj.P.Val) %>%
        mutate(across(c(P.Value, adj.P.Val), formatC, digits=2)) %>%
        mutate(across(logFC, formatC, digits=1)) %>%
        dplyr::rename(`P-value` = P.Value, `adjusted P-value` = adj.P.Val, `log2(fold change)` = logFC, chromosome = chr)
      datatable(showTab,selection = 'single', rownames = FALSE, caption = "Differentially expressed proteins (click row to show boxplot)") 
    }
  })

  
  #histogram plot or boxplot on the first panel
  output$plot1 <- renderPlotly({
      if (ifHistogram$value == TRUE) {
        #plot the histogram of p values
        p <- ggplot(corrRes(), aes(x=P.Value)) + geom_histogram(fill = "blue",col = "red", alpha=0.5, bins = 50) +
          ggtitle("P value histogram") + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                                                            plot.margin = margin(15,0,15,0)) 
        ggplotly(p) %>% config(displayModeBar = F)
        
      } else {
        lastClicked <- input$DEtab_row_last_clicked
        geneID <- filterDE()[lastClicked,]$id
        geneSymbol <- filterDE()[lastClicked,]$symbol
        plotTab <- tibble(group = geneMat()[,input$seleFeature],
                              value = assay(protSub())[geneID, ], 
                              patID = colnames(protSub())) %>%
          filter(!is.na(group),!is.na(value))
        p <- ggplot(plotTab, aes(x=group, y = value, fill = group, label = patID)) + 
          geom_boxplot(width = 0.5, alpha = 0.5, outlier.shape = NA) + 
          geom_quasirandom(width = 0.2, method = "smiley")  + 
          ylab("Normalized expression") + xlab(input$seleFeature) + 
          ggtitle(geneSymbol) + theme_bw() + theme(text=element_text(size=15),plot.title = element_text(hjust = 0.5))
        
        #hide outliers in the boxplot
        p <- plotly_build(p)
        p$x$data[[1]]$marker$opacity <- 0
        p$x$data[[2]]$marker$opacity <- 0
        
        ggplotly(p) %>% config(displayModeBar = F)
      }
  })
  
  #show the sample size information when using binary features
  output$infoBinary <- renderTable({
    testVec <- geneMat()[,input$seleFeature]
    testVec <- testVec[!is.na(testVec)]
    #data.frame(total = length(testVec),
    #           group1 = sum(testVec == levels(testVec)[1]),
    #           group2 = sum(testVec == levels(testVec)[2]))
    tt <- data.frame(table(testVec))
    colnames(tt) <- c("group","sample number")
    tt
  })


########Panel 4: Enrichment analysis###########################
  
  #to store the color values used to show gene sets that contain a certain gene
  colorList <- reactiveValues(cols = c())
  
  #a reactive variable to store GSE result
  GSEres <- reactiveValues(resTab=NULL,resObj=NULL)
  
  #a value to check wehther enrichment tab is clicked
  clickRecord <- reactiveValues(enrich = FALSE, gene = FALSE)
  
  #function to run GSEA
  resGSEA <- observeEvent(input$Run, {
    
    if (!is.null(filterDE())) {
      output$errMsg <- renderText("")
      withProgress(message = "Running enrichment analysis, please wait...", {
        
        #set color list to empty
        colorList$cols <- NULL
        
        #readin geneset database
        inGMT <- loadGSC(paste0("msigs/",input$sigSet),type="gmt")
        
        #method
        gseMethod <- input$enrichMethod
        
        #parameters for GSEA
        if (gseMethod == "GSEA") nPerm <- input$permNum
        
        #proceccing differential expression
        corTab <- filterDE() %>% arrange(P.Value) %>% filter(!duplicated(symbol), !is.na(symbol)) %>% arrange(t)
        
        #enrichment analysis
        
        if(input$statType == "stat") {
          myCoef <- data.frame(row.names = corTab$symbol, stat = corTab$t, stringsAsFactors = FALSE)
        } else {
          myCoef <- data.frame(row.names = corTab$symbol, stat = corTab$logFC, stringsAsFactors = FALSE)
        }
        
        if (gseMethod == "PAGE") {
          res <- runGSA(geneLevelStats = myCoef,geneSetStat = "page",adjMethod = "fdr", gsc=inGMT, signifMethod = 'nullDist')
        } else if (gseMethod == "GSEA") {
          res <- runGSA(geneLevelStats = myCoef,geneSetStat = "gsea",adjMethod = "fdr", gsc=inGMT, signifMethod = 'geneSampling', nPerm = nPerm)
        }
        
        resTab <- GSAsummaryTable(res)
        colnames(resTab) <- c("Name","Gene Number","Stat","p.up","p.up.adj","p.down","p.down.adj","Number up","Number down")
        if(input$ifFDR) {
          resTab <- filter(resTab, p.up.adj <= input$sigLevel | p.down.adj <= input$sigLevel) %>% arrange(desc(Stat))
        } else {
          resTab <- filter(resTab, p.up <= input$sigLevel | p.down <= input$sigLevel) %>% arrange(desc(Stat))
        }
        
  
        GSEres$resTab <- resTab
        GSEres$resObj <- res
        clickRecord$enrich <- FALSE
        clickRecord$gene <- FALSE
      })
    } else {
      output$errMsg <- renderText("Need to perform differential test on the Genomic assocation tab first!")
    }
  })
  
  
  #show differential expressed genes or enrichment results on table 1 (up table)
  output$enrichTab <- DT::renderDataTable({
    
      if( !is.null (GSEres$resTab)) {
        resTab <- GSEres$resTab
        if (!is.null(colorList$cols)) {  
          #when the bottom table is clicks, color the pathways according to whether it contains the clicked gene 
          datatable(resTab,selection = 'single', caption = "Enriched gene sets (click a row to show the proteins in the geneset)") %>%
            formatStyle('Stat',background=styleInterval(c(0),c("lightblue","pink"))) %>%
            formatStyle('Name', color = styleEqual(resTab$Name,colorList$cols)) %>%
            formatRound(c('Stat','p.up','p.up.adj', 'p.down','p.down.adj'), digits=3)
        } else { 
          #when bottom table was not clicked, do not show colors
          datatable(resTab,selection = 'single', caption = "Enriched gene sets (click a row to show the proteins in the geneset)") %>% 
            formatStyle('Stat',background=styleInterval(c(0),c("lightblue","pink"))) %>%
            formatRound(c('Stat','p.up','p.up.adj', 'p.down','p.down.adj'), digits=3)
        }
      }
  })
  
  #find gene sets that contain a certain gene
  setGene <- reactive({
    resTab <- GSEres$resTab
    setList <- loadGSC(paste0("msigs/",input$sigSet),type="gmt")$gsc[resTab$Name]
    genes <- filterDE()$symbol
    allSets <- sapply(genes, function(geneName) names(setList)[sapply(setList, function(x) geneName %in% x)])
    allSets <- allSets[sapply(allSets, function(x) length(x) != 0)]
    allSets
  })
  
  #list genes that enriched in a certain gene set
  gseaList <- reactive({
    corGene <- filterDE()
    setName <- GSEres$resTab[as.integer(input$enrichTab_row_last_clicked),"Name"]
    geneList <- loadGSC(paste0("msigs/",input$sigSet),type="gmt")$gsc[[setName]]
    geneTab <- corGene[corGene$symbol %in% geneList,]
    geneTab$setNum <- sapply(geneTab$symbol, function(x) length(setGene()[[x]]))
    geneTab <- select(geneTab, id, symbol, chr, logFC, P.Value, adj.P.Val, setNum)
    geneTab
  })
  
  #if the enrichtab is clicked, cancel the color
  observeEvent(input$enrichTab_row_last_clicked, {
    colorList$cols <- NULL
    clickRecord$enrich <- TRUE
  })
  
  #get the gene ID as well as set color of the gene sets when click the bottom table
  observeEvent(input$geneTab_row_last_clicked, {
    clickSym <- gseaList()[as.integer(input$geneTab_row_last_clicked),"symbol"][[1]]
    colorList$cols <- sapply(GSEres$resTab$Name,function(x) ifelse(x %in% setGene()[[clickSym]],"red","black"))
    clickRecord$gene <- TRUE
  })
  
  #using the bottom right table to show genes enriched in the selected set
  output$geneTab <- DT::renderDataTable({
    if (clickRecord$enrich) {
      datatable(select(gseaList(), symbol, logFC, P.Value, adj.P.Val),
                selection = 'single', rownames = FALSE, 
                caption = "Genes in the selected set (click a row to show expression boxplot)") %>% 
        formatRound(c('logFC','adj.P.Val','P.Value'))
    }
  })
  
  output$plot2 <- renderPlotly({
    if (clickRecord$gene) {
      lastClicked <- input$geneTab_row_last_clicked
      geneID <- gseaList()[lastClicked,]$id
      geneSymbol <- gseaList()[lastClicked,]$symbol
      plotTab <- data.frame(group = geneMat()[,input$seleFeature],
                            value = assay(protSub())[geneID, ], patID = colnames(protSub())) %>%
        filter(!is.na(group),!is.na(value))
      p<- ggplot(plotTab, aes(x=group, y = value, fill = group, label = patID)) + 
        geom_boxplot(width = 0.5, alpha = 0.5, outlier.shape = NA) + 
        geom_quasirandom(width = 0.2, method = "smiley")  + 
        ylab("Normalized expression") + xlab(input$seleFeature) + 
        ggtitle(geneSymbol) + theme_bw() + theme(text=element_text(size=15),plot.title = element_text(hjust = 0.5))
      
      #hide outliers in the boxplot
      p <- plotly_build(p)
      p$x$data[[1]]$marker$opacity <- 0
      p$x$data[[2]]$marker$opacity <- 0
      ggplotly(p) %>% config(displayModeBar = F)
    } else {
      p <- ggplot(data.frame()) + theme_minimal()
      ggplotly(p) %>% config(displayModeBar = F)
    }
  })
  
  #link for download the enriched set list
  output$downloadUI1 <- renderUI({
    if( !is.null (GSEres$resTab)) {
      downloadLink("downloadSet","Download enriched set list")
    }
  })
  
  output$tableLegendUI <- renderUI({
    if( !is.null (GSEres$resTab)) {
      h6("In the table below, 'Gene Number' indiactes the number of genes in each set, 'Stat' is the enrichment statistics, 'p.up' and 'p.up.adj' indicate the enrichment p-values and adjusted p-values for up-regulated genes while 'p.down' and 'p.down.adj' are for down-regulted genes. 'Number up' and 'Number down' indicate the number of up-/down-regulated genes in the set.")
    }
  })
  
  
  output$downloadUI2 <- renderUI({
    if (clickRecord$enrich) {
      downloadLink("downloadGene","Download gene list")
    }
  })
  
  #a link to download enrichment results as csv
  output$downloadSet <- downloadHandler(
    filename = function() { paste('geneSetList', '.csv', sep='') },
    content = function(file) {
      write.csv2(GSEres$resTab, file)
    }
  )
  
  
  #a link to download genes in the clicked set as csv
  output$downloadGene <- downloadHandler(
    filename = function() { paste('geneList', '.csv', sep='') },
    content = function(file) {
      write.csv2(gseaList(), file)
    }
  )
  
  

###############################################################################
  
#######Panel 3: Heatmap of differentially expressed genes######################
  output$colAnnoBox <- renderUI({
    showFeatures <- colnames(geneMat())
    selectInput("colAnno","Select additional features", showFeatures, 
                multiple = TRUE, size = 5, selectize = FALSE, selected = NULL)
  })
  
  output$setListBox <- renderUI({
    #readin geneset database
    inGMT <- sort(names(loadGSC(paste0("msigs/",input$sigSet1),type="gmt")$gsc))
    selectInput("seleHeatmapSet","Select gene sets (multiple selection allowed)", inGMT, size =8, multiple = TRUE, selectize = FALSE)
  })
  
  genesInSet <- reactive({
    #to get genes in selected gene sets
    allSets <- loadGSC(paste0("msigs/",input$sigSet1),type="gmt")$gsc
    symbolList <- unique(unlist(allSets[input$seleHeatmapSet], use.names = FALSE))
    idList <- na.omit(rownames(protNorm())[match(symbolList,rowData(protNorm())$hgnc_symbol)])
    setName <- ifelse(length(input$seleHeatmapSet) ==1,input$seleHeatmapSet,"")
    return(list(idList = idList, setName = setName))
  })
  
  orderID <- reactive({
    exprMat <- assay(protSub())
    sds <- apply(exprMat,1,sd)
    return(names(sort(sds, decreasing = TRUE)))
  })
  
  plotMap <- reactive({
    
    if (input$chooseType == "Top variant") { #plot top variant genes
      geneIDs <- orderID()[seq(1,as.integer(input$numGenes))]
      if (input$ifSet) geneIDs <- geneIDs[geneIDs %in% genesInSet()$idList]
      setName <- sprintf("Top %s most variant proteins", input$numGenes)
      if (!is.null(input$seleFeature)) 
        feature <- geneMat()[,input$seleFeature] 
      else feature <- geneMat()[,"IGHV.status"]
      names(feature) <- rownames(geneMat())
      exprMat <- assay(protSub()[geneIDs,])
    } else {  #plot differentially expressed genes
      if(!is.null(corrRes())) {
        
        #prepare the data matrix
        geneIDs <- arrange(filterDE(), t)$id
        setName <- "Differentially expressed proteins"
        if (input$ifSet) {
          geneIDs <- geneIDs[geneIDs %in% genesInSet()$idList]
          setName <- genesInSet()$setName
        }
          
        exprMat <- assay(protSub()[geneIDs,])
        feature <- geneMat()[,input$seleFeature]
        names(feature) <- rownames(geneMat())
      } else return(NULL)
    }

    #hierarchical clustering on colmns
    if (input$ifHCcol) {
      #hierarchical clustering on column
      hc <- hclust(as.dist(1-cor(exprMat, method = "spearman", use = "pairwise.complete.obs")), method = "ward.D2")
    } else {
      #order columns by feature levels
      feature <- sort(feature)
      exprMat <- exprMat[,names(feature)]
    }
    
    #manual row normalization
    exprMat <- t(scale(t(exprMat)))
    exprMat[exprMat > 4] <- 4
    exprMat[exprMat < -4] <- -4 
    
    #prepare column annotations
    annoCol <- data.frame(row.names = names(feature), test = feature)
    colnames(annoCol) <- input$seleFeature
    annoAdd <- geneMat()[names(feature), input$colAnno, drop = FALSE]
    annoCol <- cbind(annoCol, annoAdd)
    
    #prepare color scale
    breaks <- seq(-4,4,length.out = 100)
    colorBar <- colorRampPalette(c("navy","white","firebrick"))(length(breaks))
    
    #get gene symbol
    geneSymbol <- rowData(protNorm()[rownames(exprMat),])$hgnc_symbol
    
    #plot heatmap, 4 possibilities
    if (input$ifHCcol & input$ifSymbol){
      #hc on columns and show gene symbol
      p <- pheatmap(exprMat, cluster_cols = hc, cluster_rows = TRUE, treeheight_row = 0, 
                    main = setName, color = colorBar, breaks = breaks,
                    scale = "none", labels_row = geneSymbol, annotation_col = annoCol, silent = TRUE)$gtable
    } else if (input$ifHCcol & !input$ifSymbol){
      p <- pheatmap(exprMat, cluster_cols = hc, cluster_rows = TRUE, treeheight_row = 0, 
                    main = setName, color = colorBar,breaks = breaks,
                    scale = "none", show_rownames = FALSE, annotation_col = annoCol, silent=TRUE)$gtable
    } else if (!input$ifHCcol & (!input$ifSymbol)){
      p <- pheatmap(exprMat, cluster_cols = FALSE, cluster_rows = TRUE, treeheight_row = 0, 
                    main = setName, color = colorBar,breaks = breaks,
                    scale = "none", show_rownames = FALSE, annotation_col = annoCol, silent = TRUE)$gtable
    } else if ((!input$ifHCcol) & input$ifSymbol){
      p <- pheatmap(exprMat, cluster_cols = FALSE, cluster_rows = TRUE, treeheight_row = 0, 
                    main = setName, color = colorBar,breaks = breaks,
                    scale = "none", labels_row = geneSymbol, annotation_col = annoCol, silent = TRUE)$gtable
    }
    p
  })
  
  #plot the heatmap
  output$plot3 <- renderPlot({
    if (!is.null(plotMap())) {
      output$errMsg1 <- renderText("")
      withProgress(message = "Plotting heatmap, please wait...", value = NULL, {
        grid.draw(plotMap())
      })
    } else {
      output$errMsg1 <- renderText("Please perform differential expression analysis first or load a previous result!")
    }
  })
  
  
  #download the heatmap plot
  output$downHeatmap <- downloadHandler(
    filename = function() { paste0("heatmap", '.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = plotMap(), 
             device = "pdf", width = input$figWidth, height = input$figHeight,
             limitsize = FALSE)
    }
  )
  




  ########Panel 3: Assocation with drug responses################
  
  
  #the selection box for choosing a binary feature
  output$sampleSubsetBox2 <- renderUI({
    availList <- colnames(geneMatAll())
    selectInput("seleSubsetFeature2","Subset samples by", availList)
  })
  
  output$statusSubsetBox2 <- renderUI({
    availStatus <- sort(unique(geneMatAll()[[input$seleSubsetFeature2]]))
    selectInput("seleSubsetStatus2", "Sample sub-group to include", availStatus)
  })
  
  
  
  
  #the selection box for choosing blocking in the model
  output$blockingBox2 <- renderUI({
    ighvTri12Mat <- geneMat()[,colnames(geneMat()) %in% c("IGHV.status","trisomy12"), drop=FALSE]
    ighvTri12Mat <- ighvTri12Mat[rownames(ighvTri12Mat) %in% drugData()$patientID,]
    keepCol <- apply(ighvTri12Mat, 2, function(x) length(unique(x[!is.na(x)]))>1)
    blockList <- colnames(ighvTri12Mat)[keepCol]
    if (length(blockList >0)) {
      checkboxGroupInput("seleBlocking2", "Blocking for", blockList, NULL)
    }
  })
  
  # reactive object to store the selected drug screening data
  drugData <- reactive({
    if (input$subsetBox2) {
      subPat <- intersect(patBack$Patient.ID[patBack[[input$seleSubsetFeature2]] %in% input$seleSubsetStatus2], colnames(protNorm()))
      useSet <- filter(screenData, patientID %in% subPat) %>% mutate_if(is.factor, droplevels)
    } else {
      useSet <- screenData %>% 
        filter(patientID %in% colnames(protNorm())) %>% mutate_if(is.factor, droplevels)
    }
    useSet
  })
  
  #show how many samples are there in each set
  output$numSample <- renderText({
    num <- length(unique(drugData()$patientID))
    sprintf("Number of samples with proteomic data: %s",num)
  })
  
  #ui to select a drug
  output$seleDrugUI <- renderUI({
    drugList <- sort(unique(drugData()$Drug))
    selectInput("seleDrug","", drugList)
  })
  
  #ui to select a protein
  output$seleProtUI <- renderUI({
    protList <- structure(rownames(protNorm()),names = rowData(protNorm())$hgnc_symbol)
    selectInput("seleProtDrug","", protList)
  })
  
  #ui to select concentration range
  output$concRangeUI <- renderUI({
    allConc <- sort(as.integer(as.character(drugData()$concIndex)))
    sliderInput("concRange","Select concentration range", 1,max(allConc),value = c(1,max(allConc)), step = 1)
  })
  
  #ui to select color feature
  output$seleFeatureUI <- renderUI({
    selePat <- unique(drugData()$patientID)
    geneTab <- filter(patAnno, Patient.ID %in% selePat) %>%
      select(Patient.ID, IGHV.status:sex) %>%
      gather(key = "feature", value = "value",-Patient.ID) %>%
      filter(!is.na(value)) %>% group_by(feature) %>%
      summarise(n=length(unique(value))) %>% arrange(feature) %>%
      filter(n>1)
    selectInput("colorSample","Color sample by", geneTab$feature, selected = "trisomy12")
  })
  
  # a reactive value to monitor if plot histogram or boxplot
  ifHistogramDrug <- reactiveValues(value = TRUE)
  
  #if table is clicked, turn the ifHistogram to false
  observeEvent(input$DrugTab_row_last_clicked,{
    ifHistogramDrug$value <- FALSE
  })
  
  
  #prepare drug response data
  drugResponse <- reactive({
    viabData <- filter(drugData(), Drug %in% input$seleDrug, concIndex %in% as.numeric(input$concRange)) %>%
      group_by(patientID) %>% summarise(viab = mean(viab)) %>% filter(!is.na(viab)) %>%
      data.frame() %>%
      column_to_rownames("patientID")
  })
  
  #drug response matrix for testing of drugs correlated to a specific protein
  viabMat <- reactive({
    if (!input$ifDrugAUC) {
      drugMat <- mutate(drugData(), Drug = paste0(Drug,"_",concIndex))
    } else drugMat <- drugData()
    
    drugMat <- drugMat %>% filter(concIndex %in% as.numeric(input$concRange)) %>%
      group_by(patientID, Drug) %>%
      summarise(viab = mean(viab)) %>%
      spread(key = patientID, value = viab) %>% data.frame() %>%
      column_to_rownames("Drug") %>% as.matrix()
    
    drugMat
  })
  
  #reactive event for calculating drugs correlated with a specific protein
  corrProtDrug <- reactive({
    
    drugMat <- viabMat()
    
    #table to record concentrations
    concTab <- drugData() %>% mutate(id = paste0(Drug,"_",concIndex)) %>%
      distinct(id, Concentration) %>%
      bind_rows(tibble(id = unique(drugData())$Drug, Concentration = NA))
    featureVec <- assay(protNorm())[input$seleProtDrug,colnames(drugMat)]
    names(featureVec) <- colnames(drugMat)
    featureVec <- featureVec[!is.na(featureVec)]
    
    
    #prepare design 
    dMat <- data.frame(row.names = names(featureVec), prot = featureVec)
    dMat <- cbind(dMat,geneMat()[names(featureVec),input$seleBlocking2,drop=FALSE])
    dMat <- model.matrix(~., dMat)
    
 
    drugMat <- drugMat[,rownames(dMat)]
    #moderate test using limma
    fit <- lmFit(drugMat,  dMat)
    fit2 <- eBayes(fit)
    corRes <- topTable(fit2, number ="all", coef = "prot", adjust.method = "BH") %>% rownames_to_column("id") %>%
      mutate(Concentration = concTab[match(id, concTab$id),]$Concentration) %>%
      select(id, Concentration, logFC, P.Value, adj.P.Val, t) %>%
      arrange(P.Value)
    ifHistogramDrug$value <- TRUE
    corRes
    
  })
  
  #reactive event for calculating differential expression
  corrDrug <- reactive( {
    #prepare input and design matrix
    protMat <- assay(protNorm()[,rownames(drugResponse())])
    stopifnot(all(colnames(protMat) == rownames(drugResponse())))
    featureVec <- drugResponse()[,1]
    
    #prepare design 
    dMat <- data.frame(row.names = rownames(drugResponse()), drug = featureVec)
    dMat <- cbind(dMat,geneMat()[rownames(dMat),input$seleBlocking2,drop=FALSE])
    dMat <- model.matrix(~., dMat)
    
    protMat <- protMat[,rownames(dMat)]
    
    #moderate test using limma
    fit <- lmFit(protMat,  dMat)
    fit2 <- eBayes(fit)
    corRes <- topTable(fit2, number ="all", coef ="drug", adjust.method = "BH") %>% rownames_to_column("id") %>%
      mutate(symbol = rowData(protNorm()[id,])$hgnc_symbol,
             chr = rowData(protNorm()[id,])$chromosome_name) %>%
      select(id, symbol, chr, logFC, P.Value, adj.P.Val, t) %>%
      arrange(P.Value)
    ifHistogramDrug$value <- TRUE
    corRes
  })  
  
  
  filterDrug <- reactive({
    
    if (!is.null(corrDrug()) && input$seleProtOrDrug == "Select a drug") {
      DEtab <- corrDrug()
    } else if (!is.null(corrProtDrug()) && input$seleProtOrDrug == "Select a protein") {
      DEtab <- corrProtDrug()
    }
    
    if(input$ifAdjustedDrug) {
      DEtab <- filter(DEtab, abs(logFC) >= input$fcFilterDrug, adj.P.Val <= as.numeric(input$pFilterDrug))
    } else {
      DEtab <- filter(DEtab, abs(logFC) >= input$fcFilterDrug, P.Value <= as.numeric(input$pFilterDrug))
    }
    
    DEtab
  })
  
  ####Output#####
  
  #table for differentially expressed genes
  output$DrugTab <- DT::renderDataTable({
    if (!is.null(filterDrug())) {
      datatable(filterDrug(),selection = 'single', rownames = FALSE, caption = "List of associations (click row to show correlation plot)") %>% 
        formatRound(c('P.Value',"adj.P.Val"),digits=3) %>% 
        formatRound(c('logFC','t'),digits=2)
    }
  })
  
  
  #histogram plot or boxplot on the first panel
  output$plotDrug <- renderPlotly({
    if (ifHistogramDrug$value == TRUE) {
      #plot the histogram of p values
      p <- ggplot(corrDrug(), aes(x=P.Value)) + geom_histogram(fill = "blue",col = "red", alpha=0.5, bins = 50) +
        ggtitle("P value histogram") + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                                                          plot.margin = margin(30,0,30,0)) 
    } else {
      lastClicked <- input$DrugTab_row_last_clicked
      if (input$seleProtOrDrug == "Select a drug") {
        geneID <- filterDrug()[lastClicked,]$id
        geneSymbol <- filterDrug()[lastClicked,]$symbol
        plotTab <- tibble(viab = drugResponse()[,1],
                          value = assay(protNorm())[geneID,rownames(drugResponse())], 
                          patID = rownames(drugResponse())) %>%
          filter(!is.na(viab),!is.na(value))
      } else {
        drugID <- filterDrug()[lastClicked,]$id
        geneID <- input$seleProtDrug
        geneSymbol <- rowData(protNorm()[geneID,])$hgnc_symbol
        plotTab <- tibble(viab = viabMat()[drugID,],
                          value = assay(protNorm())[geneID, colnames(viabMat())],
                          patID = colnames(viabMat()))
      }
      plotTab[[input$colorSample]] <- patAnno[match(plotTab$patID, patAnno$Patient.ID),][[input$colorSample]]
      p <- ggplot(plotTab, aes(x=viab, y = value, label = patID)) + 
        geom_point(aes_string(col = input$colorSample )) + 
        geom_smooth(method = "lm") +
        ylab("Normalized proetin expression") + xlab("Viability after treatment") + 
        ggtitle(ifelse(input$seleProtOrDrug == "Select a drug", geneSymbol, drugID)) + 
        theme_bw() + theme(text=element_text(size=15),
                           plot.title = element_text(hjust = 0.5))
      
      
    }
   ggplotly(p) %>% config(displayModeBar = F)
  })


})