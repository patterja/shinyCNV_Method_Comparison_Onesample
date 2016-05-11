
shinyUI(fluidPage(
h3("Method Comparison with One Sample"),
p("Shiny interface compares two copy number",
  span("ExomeDepth [Plagnol,et al, 2012]",style="color:darkgreen"), "and",
  span("ExomeCNV [Sathirapongsasuti,etal, 2011]",style="color:purple"),
  "variant callers performance in terms of coverage and other parameters for a single whole exome sequencing sample."),
#----------------------------------------------------------------Between Tumors
tabsetPanel(
  tabPanel("Tumor1 Method Comparison", 
     
     fluidRow(  
       column(7,
          plotOutput("circos", width=550, height=550)),
       column(4,offset = 1,
          
          h4("Depth of Coverage Tumor1"),
          plotOutput("plot_doc1",height = 200),
          sliderInput('depth1', 'Minimum Depth of Coverage threshold', 
                      min = min(doc1$coverage, na.rm=TRUE), 
                      max = max(doc1$coverage, na.rm=TRUE), 
                      value = c(min(doc1$coverage, na.rm=TRUE),max(doc1$coverage, na.rm=TRUE))),
          submitButton("UPDATE"),
          column(6, 
            radioButtons("sigfilter", label = h4("Filter reads ratio"),
                         choices = list("reads ratio < 1 or > 1" = 0.5, "No Filter" = 0), 
                         selected = 0),
            submitButton("UPDATE")),
          column(6,
            radioButtons("circdata", label = h4("View Overlaps"),
                       choices = list("All Results"="allresults", "Overlaps Only" = "overlap"), 
                       selected= "allresults"),
            htmlOutput("text1"),
            conditionalPanel(condition = "input.circdata == 'overlap'",
               checkboxGroupInput("ovfilter",label = h4("Overlaps by Type"), 
                                  choices=list("ED contained in EC"="contained",
                                               "EC inside ED"="inside",
                                               "overlapdown of ED"="oldown",
                                               "overlapup of ED"="olup"),
                                  selected=c("contained","inside","oldown", "olup"))),
            submitButton("UPDATE")))),
     fluidRow(
       column(5,h4("ExomeDepth and ExomeCNV Stats"),tableOutput("depthvalues")),
       column(7,h4("ExomeDepth and ExomeCNV Venn"),h5("venn overlap if any ranges overlap"),
              d3vennROutput("depth_venn")),
       hr()),
     
     
     tabsetPanel(
       tabPanel("Overlap",dataTableOutput("overlap_methtab")),
       tabPanel("ExomeDepth Tumor1", dataTableOutput("ed1_methtab")),
       tabPanel("ExomeCNV Tumor 1", dataTableOutput("ec1_methtab"))
     ))

)))