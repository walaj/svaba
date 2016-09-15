
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#
# 
library(shiny)
library(plotly)
library(knitr)
source("helper.R")

shinyServer(function(input, output) {

  output$dummmy <- renderPlot({plot(1:10)})
    
  output$markdown <- renderUI({
    HTML(markdown::markdownToHTML(knit('README.Rmd', quiet = TRUE)))
  })
  
  output$recipes <- renderUI({
    HTML(markdown::markdownToHTML(knit('RECIPES.Rmd', quiet = TRUE)))
  })
  
  output$snowman_NA12878_table <- renderDataTable(snowman12878$dels, options=list(pageLength = 10))
  output$truth_NA12878_table <- renderDataTable(truth.NA12878$dt, options=list(pageLength = 10))
  output$truth_NA12878_table2 <- renderDataTable(truth2.NA12878$dt, options=list(pageLength = 10))
  output$pindel_NA12878_table <- renderDataTable(pindel12878$dels, options=list(pageLength = 10))
  output$lumpy_NA12878_table <- renderDataTable(lumpy12878$dels, options=list(pageLength = 10))
  
  output$simulated_events_table <- renderDataTable(events.d1, options=list(pageLength=10, autoWidth=TRUE))
  
  getDellySimData <- reactive({
    
    if (input$sim_select == "DELLY 2X") {
      ff <- readRDS("data/delly_sim.rds")[[3]]
    } else if (input$sim_select == "DELLY 5X") {
      ff <- readRDS("data/delly_sim.rds")[[4]]
    } else if (input$sim_select == "DELLY 10X") {
      ff <- readRDS("data/delly_sim.rds")[[2]]
    }
    
    if (input$tum_support_sim == "PASS-ONLY") {
      talt = 0;
      nalt = 0;
    } else {
      talt = as.numeric(input$tum_support_sim)
      nalt = 0;
    }

    mcols(ff)$TALT <- mcols(ff)$TSPLIT + mcols(ff)$TDISC
    ff <- ff[mcols(ff)$TALT >= talt]
    
    return(ff)
    
  })
  
  getStrelkaSimData <- reactive({
    
    if (input$sim_select == "Strelka 2X") {
      ff <- readRDS("data/strelka_sim.rds")[[3]]
    } else if (input$sim_select == "Strelka 5X") {
      ff <- readRDS("data/strelka_sim.rds")[[4]]
    } else if (input$sim_select == "Strelka 10X") {
      ff <- readRDS("data/strelka_sim.rds")[[2]]
    }
    
    if (input$tum_support_sim == "PASS-ONLY") {
      talt = 0;
      nalt = 0;
    } else {
      talt = as.numeric(input$tum_support_sim)
      nalt = 0;
    }
    
    ff$SPAN <- ff$CALC_SPAN
    
    ff <- ff[ff$TALT >= talt]
    
    return(ff)
    
  })
  
  getPlatypusSimData <- reactive({
    print("...getting platypus data")
    
    if (input$sim_select == "Platypus 2X") {
      ff <- readRDS("data/platypus_sim.rds")[[3]]
    } else if (input$sim_select == "Platypus 5X") {
      ff <- readRDS("data/platypus_sim.rds")[[4]]
    } else if (input$sim_select == "Platypus 10X") {
      ff <- readRDS("data/platypus_sim.rds")[[2]]
    }
    
    if (input$tum_support_sim == "PASS-ONLY") {
      talt = 0;
      nalt = 0;
    } else {
      talt = as.numeric(input$tum_support_sim)
      nalt = 0;
    }
    
    ff <- ff[ff$NV_t >= talt]
    
    return(ff)
    
  })
  
  getLumpySimData <- reactive({
    print("...getting lumpy data")
    
    if (input$sim_select == "Lumpy 2X") {
      ff <- readRDS("data/lumpy_sim.rds")[[3]]
    } else if (input$sim_select == "Lumpy 5X") {
      ff <- readRDS("data/lumpy_sim.rds")[[4]]
    } else if (input$sim_select == "Lumpy 10X") {
      ff <- readRDS("data/lumpy_sim.rds")[[2]]
    }
    
    if (input$tum_support_sim == "PASS-ONLY") {
      talt = 0;
      nalt = 0;
    } else {
      talt = as.numeric(input$tum_support_sim)
      nalt = 0;
    }
    
    mcols(ff)$TUMALT <- mcols(ff)$SR + mcols(ff)$PE
    ff <- ff[mcols(ff)$TUMALT >= talt]
    mcols(ff)$SPAN <- mcols(ff)$span
    
    print(paste(c("TALT", talt)))
    print(paste("LUMPY LEN " , length(ff)))
    #snowi = ff[ff$TUMALT >= talt & ff$NORMAL <= nalt]
    #print(length(snowi))
    #snow  = GRangesList()
    
    return(ff)
    
    
  })
  
  getPindelSimData <- reactive({
    
    print("...getting pindel data")
    if (input$sim_select == "Pindel 2X") {
      ff <- readRDS("data/pindel_sim.rds")[[3]]
    } else if (input$sim_select == "Pindel 5X") {
      ff <- readRDS("data/pindel_sim.rds")[[4]]
    } else if (input$sim_select == "Pindel 10X") {
      ff <- readRDS("data/pindel_sim.rds")[[2]]
    }
    
    if (input$tum_support_sim == "PASS-ONLY") {
      talt = 0;
      nalt = 0;
    } else {
      talt = as.numeric(input$tum_support_sim)
      nalt = 0;
    }
    ff = ff[ff$TUMALT >= talt & ff$NORMAL <= nalt]
    print(paste("Filtered PINDEL length",length(ff)))
    
    return(ff)
    
  })
  
  getSnowmanSimData <- reactive({

    if (input$sim_select=="Snowman 10X") {
      ff <- readRDS("data/160413_10x_sim.snowman.bps.rds")
      LOD=8
      #snow  <- .load_snowman_vcf_datatable("data/full88_10xsim.snowman.somatic.sv.vcf")
      #snowi <- load_indel("data/full88_10xsim.snowman.somatic.indel.vcf")
    } else if (input$sim_select == "Snowman 5X") {
      ff <- readRDS("data/160413_5x_sim.snowman.bps.rds")
      LOD=5
      #snow  <- .load_snowman_vcf_datatable("data/full88_5xsim.snowman.somatic.sv.vcf")
      #snowi <- load_indel("data/full88_5xsim.snowman.somatic.indel.vcf")
    } else if (input$sim_select == "Snowman 2X") {
      ff <- readRDS("data/160413_2x_sim.snowman.bps.rds")
      LOD = 2
      #snow  <- .load_snowman_vcf_datatable("data/160413_2x_sim.snowman.somatic.sv.vcf")
      #snowi <- load_indel("data/160413_2x_sim.snowman.somatic.indel.vcf")
    } 
   # cix <- evidence != "COMPL"
    if (input$tum_support_sim == "PASS-ONLY") {
      snowi = ff[(confidence == "PASS"| true_lod > LOD) & somatic_score & evidence == "INDEL"]
      snow  = ff[confidence == "PASS" & somatic_score & evidence != "INDEL" & evidence != "COMPL"]
    } else {
      snowi = ff[( (T_ALT >= as.integer(input$tum_support_sim) & confidence %in% c("LOWLOD","PASS") & N_ALT < 2)) & evidence == "INDEL"]
      snow <- ff[((T_ALT >= as.integer(input$tum_support_sim) & confidence %in% c("WEAKASSEMBLY", "LOWSUPPORT", "PASS") & N_ALT < 2)) & evidence != "INDEL"]
    }
    
    if (input$sim_plot_evdnc %in% c("ASSMB","ASDIS","DSCRD","COMPL"))
      snow  = snow[evidence == input$sim_plot_evdnc]
    
    snowi <- with(snowi, GRanges(chr1, IRanges(pos1, pos2), SPAN=span, ALT=T_ALT))
    snow$id <- seq(nrow(snow))
    grl.snow <- with(snow, GRanges(c(chr1,chr2), IRanges(c(pos1,pos2), width=1), strand=c(strand1,strand2), id=rep(id, 2)))
    grl.snow <- split(grl.snow, grl.snow$id)
    mcols(grl.snow)$SPAN <- snow$span
    mcols(grl.snow)$ALT <- snow$T_ALT

    ## load the snowman data
    #grl.snow <- with(snow, GRanges(chr, IRanges(pos,pos), strand=strand, id=RARID))
    #grl.snow <- split(grl.snow, grl.snow$id)
    #mcols(grl.snow) <- snow[!duplicated(RARID),.(SPAN,RARID)]

    return(list(sv=grl.snow, indel=snowi))
  })

  getFlagPlotOutput <- reactive({
    
    if (grepl("Snowman", input$sim_select)) {
      dat <- getSnowmanSimData()
      if (input$sim_plot_evdnc %in% c("ALL-SV","ASSMB","ASDIS","DSCRD","COMPL"))
        flagp <- flag.plot(xsv=dat$sv, e=gr.events)
      else if (input$sim_plot_evdnc == "INDEL")
        flagp <- flag.plot(xindel=dat$indel, e=gr.events)
      else
        flagp <- flag.plot(xsv=dat$sv, xindel=dat$indel, e=gr.events)
    } else if (grepl("Pindel", input$sim_select)) {
      flagp <- flag.plot(xindel=getPindelSimData(), e=gr.events)
    } else if (grepl("Lumpy", input$sim_select)) {
      flagp <- flag.plot(xsv=getLumpySimData(), e=gr.events)
    } else if (grepl("Platypus", input$sim_select)) {
      flagp <- flag.plot(xindel=getPlatypusSimData(), e=gr.events)
    } else if (grepl("Strelka", input$sim_select)) {
      flagp <- flag.plot(xindel=getPlatypusSimData(), e=gr.events)
    } else if (grepl("DELLY", input$sim_select)) {
      flagp <- flag.plot(xsv=getDellySimData(), e=gr.events)
    }
    
    return (flagp)
    
  })
  
  output$sim_plot <- renderPlot({
    
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Sim Plot", value = 0)

    progress$inc(1/4, detail = paste("...loading data"))
    dat <- getFlagPlotOutput()

    progress$inc(3/4, detail = paste("...rendering plot"))
    dat$g

  })
  
  output$na12878_set <- renderPlot({
    
    df <- data.frame(TP=c(snowman12878$TP1,
                          snowman12878$TP2,
                          pindel12878$TP1,
                          pindel12878$TP2,
                          lumpy12878$TP1,
                          lumpy12878$TP2,
                          nrow(truth.NA12878$dt),
                          nrow(truth2.NA12878$dt)),
                     CALLER=c("Snowman","Snowman","Pindel","Pindel","LUMPY", "LUMPY", "Truth","Truth"),
                     SET=rep(c("Deletion Set 1","Deletion Set 2"),4),
                     FP=c(snowman12878$FP1,snowman12878$FP2, pindel12878$FP1, pindel12878$FP2, lumpy12878$FP1, lumpy12878$FP2, 0,0))
    
    g <- ggplot(data=melt(df, id.vars=c("SET","CALLER"))) + geom_bar(position='dodge',aes(x=CALLER, fill=variable, y=value), stat='identity') + facet_wrap(~ SET) + 
      ylab("Event count") + xlab("") + scale_y_continuous(breaks=seq(0,12000, by=1000)) + scale_fill_manual(name="", values=c("TP"="red","FP"="grey"), labels=c("TP"="Detected deletion","FP"="Unvalidated call")) +
      geom_text(aes(x=CALLER, y=value, label=value),position=position_dodge(width=1), stat='identity') + coord_flip()
    print(g)
    #p <- ggplotly(g)
    #p
    
  })
  
  output$sim_type_plot <- renderPlot({
    
    d <- getFlagPlotOutput()

    df <- as.data.frame(gr.events)
    df$TP <- d$TP
    df$type[df$type=="RAR" & nchar(df$ins_seq) >= 10] <- "RARI"
    g <- ggplot(data=df, aes(x=type, fill=TP)) + geom_bar(position='dodge') +
      scale_fill_manual(values=c("TRUE"="darkgreen","FALSE" ="grey"), labels=c("TRUE"="TP","FALSE"="FN"), name='') +
      xlab("") + ylab("Event Count") + scale_y_continuous(expand = c(0, 0)) + 
      scale_x_discrete(expand = c(0, 0), labels=c("RARI"="Rearrangement (>= 10 bp ins)", "RAR"="Rearrangement","INT"="Translocation","ins"="insertion < 100bp", "del"="deletion < 100 bp", "DEL"="Deletion > 100 bp", "DUP"="Tandem duplication")) + coord_flip()
    print(g)
    #p <- ggplotly(g)
    #p
    
  })

  output$dum <- renderPlot({
    plot(1:10)
  })

})