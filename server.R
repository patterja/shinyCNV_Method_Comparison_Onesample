shinyServer(function(input, output) {

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Depth of Coverage 
output$plot_doc1 <- renderPlot({
  
  ggplot(doc1, aes(x=coverage)) +
    geom_histogram(binwidth=1) +
    geom_rect(data=data.frame(xmin=min(input$depth1), xmax=max(input$depth1)), 
              aes(xmin=xmin, xmax=xmax, ymin=0, ymax=Inf),
              color="red", alpha=0.5, inherit.aes = FALSE) +
    theme_bw() +
    theme(axis.text.x = element_text(size=15),axis.text.y = element_text(size=15),
          axis.title.x = element_text(size=16), axis.title.y=element_text(size=16)) +
    labs(x="coverage", y="frequency")
})

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ED abbrevation text for help
output$text1 <- renderUI({
  HTML(paste("ED= ExomeDepth","EC=ExomeCNV",sep="<br/>"))
})
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Data for between Tumors
doc1_filt<- reactive({
  makeGRangesFromDataFrame(filter(
    doc1, coverage >= min(input$depth1) & coverage <= max(input$depth1)), 
    keep.extra.columns = TRUE)})


ed1_filt<- reactive({ 
  ed1_int <- ed1  %>% 
    filter(reads.ratio > (1+as.numeric(input$sigfilter))|reads.ratio <= (1-as.numeric(input$sigfilter))) %>%
    filter(!grepl('chrX|chrY', chromosome)) 
  subsetByOverlaps(query=makeGRangesFromDataFrame(ed1_int,keep.extra.columns = TRUE),subject=doc1_filt())})


ec1_filt<- reactive({
  ec1_int <- ec1  %>% 
    filter(ratio > (1+as.numeric(input$sigfilter))|ratio <= (1-as.numeric(input$sigfilter))) %>%
    filter(!grepl('chrX|chrY', chr)) 
  subsetByOverlaps(query=makeGRangesFromDataFrame(ec1_int,keep.extra.columns = TRUE),subject=doc1_filt())})


ovfilt_method <- reactive({
  ov<-olRanges(query=ed1_filt(),subject=ec1_filt(), output="df")
  filter(ov, OLtype %in% input$ovfilter)})

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Circos plots

output$circos <- renderPlot({
  if(input$circdata=="overlap"){
    ed1gr <- ed1_filt()[ovfilt_method()$Qindex]
    seqlevels(ed1gr,force=TRUE)<-seqlevels(hg19sub)
    seqinfo(ed1gr,force=TRUE)<-seqinfo(hg19sub)
    ec1gr <- ec1_filt()[ovfilt_method()$Sindex]
    seqlevels(ec1gr,force=TRUE)<-seqlevels(hg19sub)
    seqinfo(ec1gr,force=TRUE)<-seqinfo(hg19sub)
    
  }else{
    ed1gr <- ed1_filt()
    seqlevels(ed1gr,force=TRUE)<-seqlevels(hg19sub)
    seqinfo(ed1gr,force=TRUE)<-seqinfo(hg19sub)
    ec1gr <- ec1_filt()
    seqlevels(ec1gr,force=TRUE)<-seqlevels(hg19sub)
    seqinfo(ec1gr,force=TRUE)<-seqinfo(hg19sub)}
  
  ggbio(radius = 120) + 
    circle(hg19sub,geom ="text",aes(label=seqnames,angle=90, hjust=1), vjust=1,size=4, buffer=5) +
    circle(hg19sub,geom ="ideo",fill="black",trackWidth=2, buffer=0) +
    circle(ec1gr,geom="rect",stat="identity",aes(y="type"),fill="darkmagenta", color="darkmagenta", trackWidth=30,radius=140) +
    circle(ed1gr,geom="rect",stat="identity",aes(y="type"),fill="forestgreen",color="forestgreen",trackWidth=30, radius=170)
    
#     circle(ec1gr,geom="rect",fill="darkmagenta",color="darkmagenta", trackWidth=30,radius=170, buffer=0) +
#     circle(ed1gr,geom="rect",fill="forestgreen",color="forestgreen",trackWidth=30, radius=140)
},height = 600,width = 600)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Venn Diagrams
output$depth_venn<-renderD3vennR({
  
    venn_tooltip2(
      d3vennR(
        data = list(
          list( sets = list(0), 
                label=paste("Tumor1", "ExomeDepth", sep="\n"), 
                size=length(ed1_filt()), gene=""),
          list( sets = list(1), 
                label=paste("Tumor1","ExomeCNV", sep="\n"), 
                size=length(ec1_filt()), gene=""), 
          list( sets = list(0,1), 
                size=nrow(ovfilt_method()), 
                gene="")
        )
      )
    )
})

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Stats Tables
output$depthvalues <-renderTable({
  lened1<-length(ed1_filt())
   ov<- olRanges(ed1_filt(), ec1_filt(), output = "gr") 
  
  data.frame(        
    Metric = c("CNV within depth range",
               "LOF within depth range",
               "GOF within depth range",
               "Average size within depth range",
               "Std. deviation within depth range",
               "Min size within depth range",
               "Median size within depth range",
               "Max size within depth range",
               "Average Number of bases overlap",
               "Average % overlap"),
    
    ExomeDepth = c(paste0(lened1," / ",nrow(ed1)),
               paste0(length(ed1_filt()$reads.ratio[ed1_filt()$type=="deletion"])," / ", 
                      length(ed1$reads.ratio[ed1$type=="deletion"])),
               paste0(length(ed1_filt()$reads.ratio[ed1_filt()$type=="duplication"])," / ",
                      length(ed1$reads.ratio[ed1$type=="duplication"])),
               mean(width(ed1_filt())),
               sd(width(ed1_filt())),
               min(width(ed1_filt())),
               median(width(ed1_filt())),
               max(width(ed1_filt())),
               mean(ov$OLlength),
               mean(ov$OLpercQ)),
    ExomeCNV = c(paste0(length(ec1_filt())," / ",nrow(ec1)),
               paste0(length(ec1_filt()$ratio[ec1_filt()$type=="deletion"])," / ", 
                      length(ec1$ratio[ec1$type=="deletion"])),
               paste0(length(ec1_filt()$ratio[ec1_filt()$type=="duplication"])," / ",
                      length(ec1$ratio[ec1$type=="duplication"])),
               mean(width(ec1_filt())),
               sd(width(ec1_filt())),
               min(width(ec1_filt())),
               median(width(ec1_filt())),
               max(width(ec1_filt())),
               mean(ov$OLlength),
               mean(ov$OLpercS)))
  }, digits=3)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Data tables

output$overlap_methtab <- renderDataTable({
   options(scipen=6, digits=2)

dafr<-cbind.data.frame(overlapType=ovfilt_method()$OLtype,
                       overlapLength=ovfilt_method()$OLlength,
                       ED_overlap=as.numeric(ovfilt_method()$OLpercQ),
                       EC_overlap=as.numeric(ovfilt_method()$OLpercS),
                       ED_chrom=ovfilt_method()$space, ED_start=ovfilt_method()$Qstart, ED_end=ovfilt_method()$Qend,
                       EC_chrom=ovfilt_method()$space,EC_start=ovfilt_method()$Sstart, ECend=ovfilt_method()$Send,
                       ED_type=ed1_filt()$type[ovfilt_method()$Qindex],
                       EC_type=ec1_filt()$type[ovfilt_method()$Sindex],
                       ED_reads.ratio=ed1_filt()$reads.ratio[ovfilt_method()$Qindex],
                       EC_ratio=ec1_filt()$ratio[ovfilt_method()$Sindex], 
                       ED_gene=ed1_filt()$gene[ovfilt_method()$Qindex],
                       ED_BF=ed1_filt()$BF[ovfilt_method()$Qindex],
                       ED_cytoband=ed1_filt()$cytoband[ovfilt_method()$Qindex])
})

output$ed1_methtab <- renderDataTable({
  select(as.data.frame(ed1_filt()),"Tumor1_chrom"=seqnames,"Tumor1_start"=start,
         "Tumor1_end"=end,"Tumor1_type"=type,"Tumor1_nexons"=nexons,"Tumor1_BF"=BF,
         "Tumor1_reads.ratio"=reads.ratio,"Tumor1_gene"=gene, "Tumor1_cytoband"=cytoband)})

output$ec1_methtab <- renderDataTable({
  select(as.data.frame(ec1_filt()),"Tumor1_chrom"=seqnames,"Tumor1_start"=start,"Tumor1_end"=end,"Tumor1_logR"=logR,
  "Tumor1_ratio"=ratio,"Tumor1_avg.coverage"=average.coverage, "cytoband"=ovCNVec.cytoband)})



})