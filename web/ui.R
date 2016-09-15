
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#
# 
library(shiny)
library(plotly)
library(knitr)
library(markdown)
#source("helper.R")

shinyUI(navbarPage(inverse=TRUE,
  title="Snowman - Structural variation detection by genome-wide local assembly",
  # tabPanel(title="sfd",
  #          includeHTML("README.html")
  #          ),
  tabPanel(title="Home", icon = icon("fa fa-home"),
           #includeCSS("style.css"),
           fluidPage(
             h1("Snowman", align="center", color="#0000FF"),
             h3("Structural variation and indel detection by genome-wide local assembly", align="center"),
             br(),br(),
             fluidRow(column(width=3,align="right",
             a(img(src="http://www.broadinstitute.org/collaboration/gcc/wp-content/uploads/2011/07/broad-logo.jpg", width=200),href="http://www.broadinstitute.org/", target="_blank")
             ), column(width=1,
               p("")
             ),
             column(width=3,align="left",
               a(img(src="https://upload.wikimedia.org/wikipedia/en/thumb/5/52/Dana-Farber_Cancer_Institute_logo.svg/1280px-Dana-Farber_Cancer_Institute_logo.svg.png", width=200),href="http://www.dana-farber.org/", target="_blank")
             ),column(width=1,
                      p("")
             ),
             column(width=3,align="left",
                    a(img(src="https://assets-cdn.github.com/images/modules/logos_page/GitHub-Logo.png", width=200),href="https://github.com/broadinstitute/SnowmanSV", target="_blank")
             )
             ),
             br(),
             h5("Snowman was developed by", align="center"),
             fluidRow(
               column(width=1,p("")),
               column(width=2, 
                      div(
                        img(src="https://avatars1.githubusercontent.com/u/6922120?v=3&s=460", height = 100, width=100, style='border-radius:50px', align='center'),
                        br(),
                        a("Jeremiah Wala", href="https://github.com/jwalabroad", target="_blank"),
                        style = 'float:center; text-align:center'
                        )
                      ),
               column(width=2,
                      div(
                        img(src="http://www.nygenome.org/wp-content/uploads/2015/10/Imielinski-headshot.jpeg", height = 100, width=100, style='border-radius:50px', align='center'),
                        br(),
                        a("Marcin Imielinski", href="http://www.nygenome.org/lab-groups-overview/imielinski-lab/", target="_blank"),
                        style = 'float:center; text-align:center'
                      )
                  ),
               column(width=2,
                      div(
                        img(src="http://springerlab.tch.harvard.edu/springer/uploads/Alumni/cheng-zhongzhang.jpg", height = 100, width=100, style='border-radius:50px', align='center'),
                        br(),
                        a("Cheng-Zhong Zhang", href="http://www.ncbi.nlm.nih.gov/pubmed/?term=Zhang%20CZ%5Bauth%5D", target="_blank"),
                        style = 'float:center; text-align:center'
                      )
               ),
               column(width=2,
                      div(
                        img(src="http://meyersonlab.dana-farber.org/uploads/5/2/2/7/52274861/6071117_orig.jpg", height = 100, width=100, style='border-radius:50px', align='center'),
                        br(),
                        a("Matthew Meyerson (PI)", href="http://meyersonlab.dana-farber.org/", target="_blank"),
                        style = 'float:center; text-align:center'
                      )
               ),
               column(width=2,
                      div(
                        img(src="rameen_pic.jpg", height = 100, width=100, style='border-radius:50px', align='center'),
                        br(),
                        a("Rameen Beroukhim (PI)", href="http://beroukhimlab.dfci.harvard.edu/Beroukhim_Laboratory.html", target="_blank"),
                        style = 'float:center; text-align:center'
                      )
               ),
               column(width=1,p(""))
             ),
             br(),
             br(),br(),
             img(src='schematic_snowman.png',style="display: block; margin-left: auto; margin-right: auto;", width="50%")
           ) ## fluid Page
  ),
   navbarMenu(title="Documentation",
     tabPanel(title="README", icon = icon(lib="font-awesome", "file-text-o"),
              #p("")
              includeHTML("README.html")
              ##rmarkdown::render('RECIPES.Rmd') ## THEN HAVE TO REMOVE <DOCTYPE> tag pairs, and everything after README and until first <style type="text/css">
              ##rmarkdown::render('RECIPES.Rmd')
              #uiOutput('markdown')
     ),
  tabPanel(title="RECIPES", icon = icon(lib="font-awesome", "file-text-o"),
          # p("")
          includeHTML("RECIPES.html")
  )
),
   navbarMenu("NA12878",
   tabPanel(title="Data tables",icon = icon("fa fa-th-list"),
          h3("NA12878 deletions from 1000 Genomes (Mills et al)"),
            div(dataTableOutput("truth_NA12878_table"), style="font-size:70%"),
            h3("NA12878 deletions from 1000 Genomes + Validated LUMPY+DELLY+GASVPro+Pindel (Layer et al)"),
            div(dataTableOutput("truth_NA12878_table2"), style="font-size:70%"),
            h3("NA12878 deletions from Snowman"),
            div(dataTableOutput("snowman_NA12878_table"), style="font-size:70%"),
            h3("NA12878 deletions from Lumpy"),
            div(dataTableOutput("lumpy_NA12878_table"), style="font-size:70%"),
            h3("NA12878 deletions from Pindel"),
            div(dataTableOutput("pindel_NA12878_table"), style="font-size:70%")
          ),
   tabPanel(title="Comparison figures", icon = icon(lib="font-awesome", "bar-chart"),
            h3("Events detected"),
            plotOutput("na12878_set")
            )
   ),
   navbarMenu("Simulated tumor",
              tabPanel(title="Events table", icon = icon("fa fa-th-list"),
                       h3("Simulated SVs"),
                     div(dataTableOutput("simulated_events_table"), style="font-size:70%")
                       ),
             tabPanel(title="Caller plots",icon= icon(lib="font-awesome", "bar-chart"),
                      selectInput("sim_select", label="Select data to plot", 
                                  choices=c("Snowman 10X", "Snowman 5X", "Snowman 2X", "Pindel 10X", "Pindel 5X", "Pindel 2X",
                                            "Lumpy 10X", "Lumpy 5X", "Lumpy 2X", "Platypus 2X", "Platypus 5X", "Platypus 10X",
                                            "Strelka 2X", "Strelka 5X", "Strelka 10X", "DELLY 2X","DELLY 5X", "DELLY 10X"), selected=""),
                      selectInput("tum_support_sim", label="Select minimum T_ALT", choices=c("PASS-ONLY",2:9), selected="PASS-ONLY"),
                      selectInput("sim_plot_evdnc", label="Select EVDNC to plot", choices=c("ALL", "INDEL", "ALL-SV", "ASSMB","ASDIS","DSCRD","COMPL"), selected="ALL"),
                   
                      plotOutput("sim_plot"),
                      plotOutput("sim_type_plot", width="80%")
                      )
   ),
   navbarMenu("HCC1143",
              tabPanel(title="Snowman", icon = icon("fa fa-th-list"),
                       h3("Snowman somatic events (101-bp Illumina)"),
                       #div(dataTableOutput("pindel_hcc1143_101"), style="font-size:70%"),
                       h3("Snowman somatic events (250-bp PCR-free Illumina)")
                       #div(dataTableOutput("pindel_hcc1143_250"), style="font-size:70%")
              ),
              tabPanel(title="LUMPY", icon = icon("fa fa-th-list"),
                       h3("LUMPY somatic events (101-bp Illumina)"),
                       #div(dataTableOutput("pindel_hcc1143_101"), style="font-size:70%"),
                       h3("LUMPY somatic events (250-bp PCR-free Illumina)")
                       #div(dataTableOutput("pindel_hcc1143_250"), style="font-size:70%")
              ),
              tabPanel(title="Pindel", icon = icon("fa fa-th-list"),
                       h3("Pindel somatic events (101-bp Illumina)"),
                       #div(dataTableOutput("pindel_hcc1143_101"), style="font-size:70%"),
                       h3("Pindel somatic events (250-bp PCR-free Illumina)")
                       #div(dataTableOutput("pindel_hcc1143_250"), style="font-size:70%")
              ),
              tabPanel(title="DELLY", icon = icon("fa fa-th-list"),
                       h3("Delly somatic events (101-bp Illumina)"),
                       #div(dataTableOutput("delly_hcc1143_101"), style="font-size:70%"),
                       h3("Delly somatic events (250-bp PCR-free Illumina)")
                       #div(dataTableOutput("delly_hcc1143_250"), style="font-size:70%")
              ),
              tabPanel(title="Strelka", icon = icon("fa fa-th-list"),
                       h3("Strelka somatic events (101-bp Illumina)"),
                       #div(dataTableOutput("strelka_hcc1143_101"), style="font-size:70%"),
                       h3("Strelka somatic events (250-bp PCR-free Illumina)")
                       #div(dataTableOutput("strelka_hcc1143_250"), style="font-size:70%")
              ),
              tabPanel(title="Platypus", icon = icon("fa fa-th-list"),
                       h3("Platypus somatic events (101-bp Illumina)"),
                       #div(dataTableOutput("platypus_hcc1143_101"), style="font-size:70%"),
                       h3("Platypus somatic events (250-bp PCR-free Illumina)")
                       #div(dataTableOutput("platypus_hcc1143_250"), style="font-size:70%")
              )
   )
   ))


