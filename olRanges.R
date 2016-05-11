
#############################
## Identify Range Overlaps ##
#############################
## Author: Thomas Girke
## Last update: 8-Feb-11
## Details on usage and use cases are available here:
## http://manuals.bioinformatics.ucr.edu/home/ht-seq#TOC-Analysis-Routines-with-IRanges-Geno

## Utility: identify overlaps in range data sets, such as annotation or alignment positions defined
## by two IRanges/GRanges objects.  

## Overlap types
## olup: startup & endin
## Q --------------
## S       -------------

## oldown: startin & enddown
## Q       -------------
## S --------------

## inside: startin & endin 
## Q     -----
## S --------------

## contained: startup & enddown
## Q --------------
## S     -----

###########################################################
## (A) olRanges Function for IRanges and GRanges Objects ##
###########################################################
olRanges <- function(query, subject, output="gr", ...) {
  require(GenomicRanges); require(IRanges)
  
  ## Input check
  if(!((class(query)=="GRanges" & class(subject)=="GRanges") | (class(query)=="IRanges" & class(subject)=="IRanges"))) {
    stop("Query and subject need to be of same class, either GRanges or IRanges!")
  }
  
  ## Find overlapping ranges
  if(class(query)=="GRanges") {
    seqlengths(query) <- rep(NA, length(seqlengths(query)))
    seqlengths(subject) <- rep(NA, length(seqlengths(subject)))
  }
  olindex <- as.matrix(findOverlaps(query, subject, ...))
  query <- query[olindex[,1]]
  subject <- subject[olindex[,2]]
  olma <- cbind(Qstart=start(query), Qend=end(query), Sstart=start(subject), Send=end(subject))
  
  ## Pre-queries for overlaps
  startup <- olma[,"Sstart"] < olma[,"Qstart"]
  enddown <- olma[,"Send"] > olma[,"Qend"]
  startin <- olma[,"Sstart"] >= olma[,"Qstart"] & olma[,"Sstart"] <= olma[,"Qend"]
  endin <- olma[,"Send"] >= olma[,"Qstart"] & olma[,"Send"] <=  olma[,"Qend"]
  
  ## Overlap types
  olup <- startup & endin
  oldown <- startin & enddown
  inside <- startin & endin 
  contained <- startup & enddown
  
  ## Overlap types in one vector
  OLtype <- rep("", length(olma[,"Qstart"]))
  OLtype[olup] <- "olup"
  OLtype[oldown] <- "oldown"
  OLtype[inside] <- "inside" 
  OLtype[contained] <- "contained"
  
  ## Overlap positions
  OLstart <- rep(0, length(olma[,"Qstart"]))
  OLend <- rep(0, length(olma[,"Qstart"]))
  OLstart[olup] <- olma[,"Qstart"][olup]
  OLend[olup] <- olma[,"Send"][olup]
  OLstart[oldown] <- olma[,"Sstart"][oldown]
  OLend[oldown] <- olma[,"Qend"][oldown]
  OLstart[inside] <- olma[,"Sstart"][inside]
  OLend[inside] <- olma[,"Send"][inside]
  OLstart[contained] <- olma[,"Qstart"][contained]
  OLend[contained] <- olma[,"Qend"][contained]
  
  ## Absolute and relative length of overlaps
  OLlength <- (OLend - OLstart) + 1
  OLpercQ <- OLlength/width(query)*100
  OLpercS <- OLlength/width(subject)*100
  
  ## Output type
  oldf <- data.frame(Qindex=olindex[,1], Sindex=olindex[,2], olma, OLstart, OLend, OLlength, OLpercQ, OLpercS, OLtype)
  if(class(query) == "GRanges") {
    oldf <- cbind(space=as.character(seqnames(query)), oldf)
  }
  if(output=="df") {
    return(oldf)
  }
  if(output=="gr") {
    if(class(query)=="GRanges") {
      elementMetadata(query) <- cbind(as.data.frame(elementMetadata(query)), oldf)
    }
    if(class(query)=="IRanges") {
      query <- GRanges(seqnames = Rle(rep("dummy", length(query))), ranges = IRanges(start=oldf[,"Qstart"], end=oldf[,"Qend"]), strand = Rle(strand(rep("+", length(query)))), oldf)  
    }
    return(query)
  }
}

#------------------------------------------------------------------------------------
venn_tooltip2<- function( venn ){
  venn$x$tasks[length(venn$x$tasks)+1] <- list(
    htmlwidgets::JS('
                    function(){
                    console.log("here");
                    var div = d3.select(this);
                    
                    // add a tooltip
                    var tooltip = d3.select("body").append("div")
                    .attr("class", "venntooltip")
                    .style("position", "absolute")
                    .style("text-align", "center")
                    .style("width", 250)
                    .style("height", 500)
                    .style("background", "#333")
                    .style("color","#ddd")
                    .style("padding","2px")
                    .style("border","0px")
                    .style("border-radius","8px")
                    .style("opacity",0);
                    
                    div.selectAll("path")
                    .style("stroke-opacity", 0)
                    .style("stroke", "#fff")
                    .style("stroke-width", 0)
                    
                    // add listeners to all the groups to display tooltip on mousover
                    div.selectAll("g")
                    .on("mouseover", function(d, i) {
                    
                    // sort all the areas relative to the current item
                    venn.sortAreas(div, d);
                    
                    // Display a tooltip with the current size
                    tooltip.transition().duration(400).style("opacity", .9);
                    tooltip.html(d.size+"<br>"+d.gene);
                    
                    // highlight the current path
                    var selection = d3.select(this).transition("tooltip").duration(400);
                    selection.select("path")
                    .style("stroke-width", 3)
                    .style("fill-opacity", d.sets.length == 1 ? .4 : .1)
                    .style("stroke-opacity", 1);
                    })
                    
                    .on("mousemove", function() {
                    tooltip.style("left", (d3.event.pageX) + "px")
                    .style("top", (d3.event.pageY - 28) + "px");
                    })
                    
                    .on("mouseout", function(d, i) {
                    tooltip.transition().duration(400).style("opacity", 0);
                    var selection = d3.select(this).transition("tooltip").duration(400);
                    selection.select("path")
                    .style("stroke-width", 0)
                    .style("fill-opacity", d.sets.length == 1 ? .25 : .0)
                    .style("stroke-opacity", 0);
                    });
                    }
                    ')
    )
  venn
  }
