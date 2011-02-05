createProtovishtmlforCV <- function(numPlots, numsamples,dates,htmlfile,jsfile,graphWidth,graphHeight){
	

	
n = numPlots
panelWidth = n*(graphWidth + 30) + 15*n;
panelHeight = 1000
offsetFromTop = 15


  cat(file = htmlfile, append = "FALSE", "\n")
  cat(file = htmlfile, append = "TRUE", "<html", "\n")
  cat(file = htmlfile, append = "TRUE", "<head", "\n")
  cat(file = htmlfile, append = "TRUE", "<title> Dependency Graph </title>", "\n")
  cat(file = htmlfile, append = "TRUE", " <link type=\"text/css\" rel=\"stylesheet\" href=\"../extdata/ex.css\"/>","\n")
  cat(file = htmlfile, append = "TRUE", "  <script type=\"text/javascript\" src=\"../extdata/protovis-r3.2.js\"></script> ", "\n")
  cat(file = htmlfile, append = "TRUE", " <script type=\"text/javascript\" src=\"",jsfile,"\"></script>", "\n")
  cat(file = htmlfile, append = "TRUE", "<style type=\"text/css\">", "\n")
  cat(file = htmlfile, append = "TRUE","#fig {", "\n")
  cat(file = htmlfile, append = "TRUE", "width:",panelWidth,"px;", "\n", sep="")
  cat(file = htmlfile, append = "TRUE", " height:",panelHeight,"px;}", "\n",sep ="")
  cat(file = htmlfile, append = "TRUE", " </style>", "\n")
  cat(file = htmlfile, append = "TRUE", "</head>", "\n")
  cat(file = htmlfile, append = "TRUE", " <body><div id=\"center\"><div id=\"fig\">", "\n")
  cat(file = htmlfile, append = "TRUE", " <script type=\"text/javascript+protovis\">", "\n")
  cat(file = htmlfile, append = "TRUE", "function average() {", "\n")
  cat(file = htmlfile, append = "TRUE", "var items = average.arguments.length", "\n")
  cat(file = htmlfile, append = "TRUE", " var sum = 0", "\n")
  cat(file = htmlfile, append = "TRUE", " for (i = 0; i < items;i++){", "\n")
  cat(file = htmlfile, append = "TRUE", " sum += average.arguments[i]}", "\n")
  cat(file = htmlfile, append = "TRUE", "return (sum/items)}", "\n")
  cat(file = htmlfile, append = "TRUE", "var  c = pv.Colors.category20(),", "\n")
  cat(file = htmlfile, append = "TRUE", " w =",graphWidth,",", "\n",sep="");
  cat(file = htmlfile, append = "TRUE", "h =",graphHeight,";", "\n",sep="")
  cat(file = htmlfile, append = "TRUE", "var origPanel = new pv.Panel()", "\n")
  cat(file = htmlfile, append = "TRUE", ".width(",panelWidth,")", "\n")
  cat(file = htmlfile, append = "TRUE", ".height(",panelHeight,")", "\n")
  cat(file = htmlfile, append = "TRUE", " .fillStyle(\"white\")", "\n")
  cat(file = htmlfile, append = "TRUE", ".event(\"mousedown\", pv.Behavior.pan())", "\n")
  cat(file = htmlfile, append = "TRUE", " .event(\"mousewheel\", pv.Behavior.zoom());", "\n")
  
  for ( i in 1:n) {
  	
  	 cat(file = htmlfile, append = "TRUE", "var vis",i, "= origPanel.add(pv.Panel)", "\n", sep = "")
    cat(file = htmlfile, append = "TRUE", ".left(", i-1, "*w+", i, "*15)","\n",sep ="")
    cat(file = htmlfile, append = "TRUE", ".width(w)", "\n")
    cat(file = htmlfile, append = "TRUE", " .height(h)", "\n")
    cat(file = htmlfile, append = "TRUE", " .top(15)", "\n")
    cat(file = htmlfile, append = "TRUE", " .fillStyle(\"white\")", "\n")
    cat(file = htmlfile, append = "TRUE", ".event(\"mousedown\", pv.Behavior.pan())", "\n")
    cat(file = htmlfile, append = "TRUE", " .event(\"mousewheel\", pv.Behavior.zoom());", "\n")
    cat(file = htmlfile, append = "TRUE", "var border",i, "= vis", i,".add(pv.Bar)", "\n", sep = "")
    cat(file = htmlfile, append = "TRUE", " .fillStyle(\"white\")", "\n")
    cat(file = htmlfile, append = "TRUE", " .strokeStyle(\"black\");", "\n")
  	
  	
    if ( i ==n) 
    {
       startdt =1
       enddt = numsamples
    }
    else 
    {
        startdt = dates[i]
        enddt = dates[i+1] -1 
    
	}
	
	if ( i == 1) 
	{
		cat(file = htmlfile, append = "TRUE", " var layout= vis",i,".add(pv.Layout.Force)", "\n", sep = "")
      cat(file = htmlfile, append = "TRUE", ".nodes(data",i,".nodes)", "\n", sep = "")
      cat(file = htmlfile, append = "TRUE", ".links(data",i,".links);", "\n", sep = "")
      cat(file = htmlfile, append = "TRUE", "layout.springLength(100)", "\n")
      cat(file = htmlfile, append = "TRUE", ".dragConstant(0.15)", "\n")
      cat(file = htmlfile, append = "TRUE", ".bound(true)", "\n")
      cat(file = htmlfile, append = "TRUE", ".chargeConstant(-500)", "\n")
      cat(file = htmlfile, append = "TRUE", ".chargeMaxDistance(80);", "\n")
      cat(file = htmlfile, append = "TRUE", "layout.link.add(pv.Line)", "\n")
      cat(file = htmlfile, append = "TRUE", ".lineWidth(function(d) 2*Math.pow(Math.abs(d.condcorr),0.2))", "\n")
      cat(file = htmlfile, append = "TRUE", ".strokeStyle(function(p,d) d.condcorr<0 ? \"red\":\"green\");", "\n")
      cat(file = htmlfile, append = "TRUE", "layout.node.add(pv.Dot)", "\n")
      cat(file = htmlfile, append = "TRUE", " .fillStyle(function(d) c(d.group))", "\n")
      cat(file = htmlfile, append = "TRUE", " .strokeStyle(function() this.fillStyle().darker())", "\n")
      cat(file = htmlfile, append = "TRUE", " .size(160)", "\n")
      cat(file = htmlfile, append = "TRUE", ".title(function(d) d.nodeName)", "\n")
      cat(file = htmlfile, append = "TRUE", ".event(\"mousedown\", pv.Behavior.drag())", "\n")
      cat(file = htmlfile, append = "TRUE", " .event(\"drag\", layout);", "\n")
      cat(file = htmlfile, append = "TRUE", "layout.node.add(pv.Label)", "\n")
      cat(file = htmlfile, append = "TRUE", ".text(function(d) d.nodeName)", "\n")
      cat(file = htmlfile, append = "TRUE", ".font(\"10px sans-serif\")", "\n")
      cat(file = htmlfile, append = "TRUE", ".textAlign(\"center\")", "\n")
      cat(file = htmlfile, append = "TRUE", ".textBaseline(\"middle\");", "\n")
      
      
            cat(file = htmlfile, append = "TRUE", "layout", "\n")

   
      cat(file = htmlfile, append = "TRUE", ".add(pv.Label)", "\n")

      cat(file = htmlfile, append = "TRUE", ".text(\"Dates:",startdt, "-", enddt , "\")", "\n", sep ="")

          cat(file = htmlfile, append = "TRUE", ".font(\"20px sans-serif\")", "\n")
      cat(file = htmlfile, append = "TRUE", ".textAlign(\"left\")", "\n")

      cat(file = htmlfile, append = "TRUE", ".left(", i-1,"*w+",i,"*15)", "\n", sep ="")


      cat(file = htmlfile, append = "TRUE", ".top(", graphHeight-15, ")", "\n", sep ="")

      cat(file = htmlfile, append = "TRUE", ".textBaseline(\"middle\");", "\n")
    
		}
	else {
		
		cat(file = htmlfile, append = "TRUE", " var layout",i,"= vis",i,".add(pv.Layout.Force)", "\n", sep = "")
      cat(file = htmlfile, append = "TRUE", ".nodes(data",i,".nodes)", "\n", sep = "")
      cat(file = htmlfile, append = "TRUE", ".links(data",i,".links);", "\n", sep = "")
      cat(file = htmlfile, append = "TRUE", "layout",i,".springLength(100)", "\n",sep="")
      cat(file = htmlfile, append = "TRUE", ".dragConstant(0.15)", "\n")
      cat(file = htmlfile, append = "TRUE", ".bound(true)", "\n")
      cat(file = htmlfile, append = "TRUE", ".chargeConstant(-500)", "\n")
      cat(file = htmlfile, append = "TRUE", ".chargeMaxDistance(80);", "\n")
      cat(file = htmlfile, append = "TRUE", "layout",i,".link.add(pv.Line)", "\n",sep="")
      cat(file = htmlfile, append = "TRUE", ".lineWidth(function(d) 2*Math.pow(Math.abs(d.condcorr),0.2))", "\n")
      cat(file = htmlfile, append = "TRUE", ".strokeStyle(function(p,d) d.condcorr<0 ? \"red\":\"green\");", "\n")
      cat(file = htmlfile, append = "TRUE", "layout",i,".node.add(pv.Dot)", "\n",sep="")
      cat(file = htmlfile, append = "TRUE", " .fillStyle(function(d) c(d.group))", "\n")
      cat(file = htmlfile, append = "TRUE", " .strokeStyle(function() this.fillStyle().darker())", "\n")
      cat(file = htmlfile, append = "TRUE", " .size(160)", "\n")
      cat(file = htmlfile, append = "TRUE", ".title(function(d) d.nodeName)", "\n")
      cat(file = htmlfile, append = "TRUE", ".event(\"mousedown\", pv.Behavior.drag())", "\n")
      cat(file = htmlfile, append = "TRUE", " .event(\"drag\", layout);", "\n")
      cat(file = htmlfile, append = "TRUE", "layout",i,".node.add(pv.Label)", "\n",sep="")
      cat(file = htmlfile, append = "TRUE", ".text(function(d) d.nodeName)", "\n")
      cat(file = htmlfile, append = "TRUE", ".font(\"10px sans-serif\")", "\n")
      cat(file = htmlfile, append = "TRUE", ".textAlign(\"center\")", "\n")
      cat(file = htmlfile, append = "TRUE", ".textBaseline(\"middle\");", "\n")
      
      
            cat(file = htmlfile, append = "TRUE", "layout",i, "\n", sep ="")

   
      cat(file = htmlfile, append = "TRUE", ".add(pv.Label)", "\n")

      cat(file = htmlfile, append = "TRUE", ".text(\"Dates:",startdt, "-", enddt , "\")", "\n", sep ="")

          cat(file = htmlfile, append = "TRUE", ".font(\"20px sans-serif\")", "\n")
      cat(file = htmlfile, append = "TRUE", ".textAlign(\"left\")", "\n")


      cat(file = htmlfile, append = "TRUE", ".top(", graphHeight-15, ");", "\n", sep ="")

    

		
		}  	
  	
  	}
  	   
    cat(file = htmlfile, append = "TRUE", "origPanel.render();", "\n")
    cat(file = htmlfile, append = "TRUE", "</script>", "\n")
    cat(file = htmlfile, append = "TRUE", "<div class=\"caption\">", "\n")
    cat(file = htmlfile, append = "TRUE", "</div>", "\n")
    cat(file = htmlfile, append = "TRUE", " </div></div></body>", "\n")
    cat(file = htmlfile, append = "TRUE", "</html>", "\n")
    return(0)

  	return(0)
  	}
	
	
	
	
	
	
