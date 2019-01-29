library(shiny)
library(shinydashboard)
library(DT)
library(shinyjs)
library(ggplot2)
options(shiny.fullstacktrace=TRUE)
ui<-dashboardPage(
  dashboardHeader(title="Mixed Model Power & Size",titleWidth=300),
  dashboardSidebar(width=300,sidebarMenu(
    id = "tabs",
    menuItem("Introduction",tabName = "intro",icon=icon("calculator")
    ),
    menuItem("RCRM",tabName = "rcrm",icon=icon("calculator")
    ),
    menuItem("Two-Stage Model",tabName = "tsm",icon=icon("calculator")
    )
  )),
  dashboardBody(
    shinyjs::useShinyjs(),
    tabItems(
      tabItem(tabName="intro",
              mainPanel(h1("Introduction"),
                        p("This R Shiny App calculates the sample size and power for the following 2 longitidinal models:"),
                        p("- Random Coefficient Regression Model (RCRM)"),
                        p("- Two-Stage Mixed Effects Model."),
                        br(),
                        br(),
                        p("The RCRM theoretical results are summarized in a manuscript (Hu, Mackey, and Thomas; in progress). The two-stage mixed effects model details can be found in Section 8.4 of Fitzmaurice, Laird, and Ware (2011)."),
                        br(),
                        br(),
                        p("This R Shiny App was co-developed by Nan Hu (Genentech Inc., Biostatistics) and a summer intern Zhe Qu (Ph.D. candidate at Tulane University). If you have any technical questions or any ideas for improving functionality, please contact Nan Hu (hu.nan@gene.com).", style = "color:blue"))
              
              
      ),
      tabItem(tabName="rcrm",
              fluidRow(
                box(width = 9,
                    title = "Parameters",
                    status = "primary",
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    fluidRow(
                      column(width = 4,
                             # div(style="height:5px;"),
                             radioButtons("level", 
                                          h5("Test significance level \U03B1"), 
                                          choices = list("1-sided" = 1, 
                                                         "2-sided" = 2),
                                          selected = 1),
                             uiOutput("side1"),
                             uiOutput("side2")
                             # numericInput("siglevel", "",value=0.025)
                             #numericInput("siglevel", div(style = "height=40px",h5("One-sided test significance level")),value=0.025)
                      ),
                      column(width = 4,
                             div(style="height:40px;",h5("Effect size (absolute difference in rate of decline between tx and placebo)")),
                             numericInput("effect", 
                                          "", 
                                          value = 0.5)
                      ),
                      column(width = 4,
                             #div(style="height:40px;",h5(paste("SD of random intercept: \U03B2"))),
                             div(style="height:40px;",h5(HTML(paste0("SD of random intercept: \U03C3",tags$sub("\U03B1"))))),
                             numericInput("sig0", 
                                          "", 
                                          value = 7))
                    ),
                    fluidRow(
                      column(width = 4,
                             numericInput("rho", 
                                          h5(HTML(paste0("Correlation between random intercept & slope: \U03C1"))), 
                                          value = 0.3, min = -1, max = 1, step = 0.1)
                      ),
                      column(width = 4,
                             numericInput("sig1", 
                                          h5(HTML(paste0("SD of random slope: \U03C3",tags$sub("\U03B2")))), 
                                          value = 1.5)
                      ),
                      column(width = 4,
                             numericInput("sig", 
                                          h5("SD of pure error: \U03C3"), 
                                          value = 3)
                      )
                    ),
                    fluidRow(
                      column(width = 4,
                      #        numericInput("h", 
                      #                     h5("Time length between any two consecutive assessments"), 
                      #                     value = 0.5)
                             radioButtons("time_length", 
                                          h5("Time length between any two consecutive assessments"), 
                                          choices = list("equal length" = 1, 
                                                         "unequal length" = 2),
                                          selected = 1),
                             uiOutput("length1"),
                             uiOutput("length2")
                      ),
                    
                      column(width = 4,
                             numericInput("n_t", 
                                          h5("Number of assessments (including baseline)"), 
                                          value = 5)
                      ),
                      column(width = 4,
                             numericInput("droprate", 
                                          h5("Exponential dropout rate"), 
                                          value = 0.01, min = 0, step = 0.1)
                      )
                    )
                ),
                box(
                  width = 3,
                  title = "Result Type",
                  status = "primary",
                  solidHeader = TRUE,
                  collapsible = TRUE,
                  radioButtons("choice", 
                               h5("What are you solving for"), 
                               choices = list("Given power, find n" = 1, 
                                              "Given n, find power" = 2),
                               selected = 2),
                  
                  uiOutput("conditionalInput"),
                  uiOutput("conditionalInput2"),
                  uiOutput("conditionalInput3"),
                  uiOutput("conditionalInput4"),
                  uiOutput("conditionalInput5"),
                  actionButton("subBtn", "Submit", width = "100%")
                  
                )
              ),
              fluidRow(
                box(
                  width = 12,
                  title = "Results",
                  status = "primary",
                  solidHeader = TRUE,
                  h4(textOutput("choice", container = span)),
                  #span(textOutput("power"),style="font-size:14pt;font-weigth:bold;"),
                  DT::dataTableOutput("rcrm_table"),
                  uiOutput("showTblDownBtn"),
                  uiOutput("showplot"),
                  uiOutput("showplotDownBtn")
                )
              )
      ),
      tabItem(tabName="tsm",
              fluidRow(
                box(width = 9,
                    title = "Parameters",
                    status = "primary",
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    fluidRow(
                      column(width = 4,
                             # div(style="height:5px;"),
                             radioButtons("tsm_level", 
                                          h5("Test significance level \U03B1"), 
                                          choices = list("1-sided" = 1, 
                                                         "2-sided" = 2),
                                          selected = 1),
                             uiOutput("tsm_side1"),
                             uiOutput("tsm_side2")
                             # numericInput("siglevel", "",value=0.025)
                             #numericInput("siglevel", div(style = "height=40px",h5("One-sided test significance level")),value=0.025)
                      ),
                      column(width = 4,
                             div(style="height:40px;",h5("Effect size (absolute difference in rate of decline between tx and placebo)")),
                             numericInput("tsm_effect", 
                                          "", 
                                          value = 0.5)
                      ),
                      column(width = 4,
                             div(style="height:40px;",h5("Duration of follow up:")),
                             numericInput("tsm_duration", 
                                          "", 
                                          value = 5))
                    ),
                    fluidRow(
                      column(width = 4,
                             numericInput("tsm_nt", 
                                          h5("Number of assessments (including baseline)"), 
                                          value = 5)
                      ),
                      column(width = 4,
                             numericInput("tsm_sig1", 
                                          h5(HTML(paste0("Between-subject SD: \U03C3",tags$sub("b")))), 
                                          value = 1.5)
                      ),
                      
                      column(width = 4,
                             numericInput("tsm_sig", 
                                          h5(HTML(paste0("Within-subject SD: \U03C3",tags$sub("w")))), 
                                          value = 3)
                      )
                    )
                    
                ),
                box(
                  width = 3,
                  title = "Result Type",
                  status = "primary",
                  solidHeader = TRUE,
                  collapsible = TRUE,
                  radioButtons("tsm_choice", 
                               h5("What are you solving for"), 
                               choices = list("Given power,find n" = 1, 
                                              "Given n,find power" = 2),
                               selected = 2),
                  
                  uiOutput("tsm_conditionalInput"),
                  uiOutput("tsm_conditionalInput2"),
                  uiOutput("tsm_conditionalInput3"),
                  uiOutput("tsm_conditionalInput4"),
                  uiOutput("tsm_conditionalInput5"),
                  actionButton("tsm_subBtn", "Submit", width = "100%")
                  
                )
              ),
              fluidRow(
                box(
                  width = 12,
                  title = "Results",
                  status = "primary",
                  solidHeader = TRUE,
                  h4(textOutput("tsm_choice", container = span)),
                  
                  DT::dataTableOutput("tsm_table"),
                  uiOutput("tsm_showTblDownBtn"),
                  uiOutput("tsm_showplot"),
                  uiOutput("tsm_showplotDownBtn")
                )
              )
          )
    )
  )
)

server<-function(input,output,session)
{
  rst_table<-data.frame(matrix(nrow=0,ncol=12))
  
  observeEvent(input$time_length,{
    if(input$time_length==2){
      shinyjs::disable("n_t")
    }else{
      shinyjs::enable("n_t")
    }
  })
######RCRM#########
  output$side1<-renderUI({
    if(input$level==1){
      numericInput("oneside", "",value=0.025)
    }
  })
  output$side2<-renderUI({
    if(input$level==2){
      numericInput("twoside", "",value=0.05)
    }
  })
  output$length1<-renderUI({
    if(input$time_length==1){
      numericInput("h", "",value=0.5)
    }
  })
  output$length2<-renderUI({
    if(input$time_length==2){
      textInput("gap_length",h5("Comma separated list of visit time points (including time 0)"),value="0,0.5,1,1.5,2")
    }
  })
  output$conditionalInput <- renderUI({
    if(input$choice==1){
      numericInput("power", 
                   h5("Expected power"), 
                   value = 0.8)
    } 
  })
  
  output$conditionalInput2<-renderUI({
    if(input$choice==2)
    {
      numericInput("n_tre", 
                   h5("Sample size of treatment group"), 
                   value = 100)
    }
  })
  
  output$conditionalInput3<-renderUI({
    if(input$choice==2)
    {
      numericInput("n_pla", 
                   h5("Sample size of placebo group"),value=100)
    }
  })
  
  output$conditionalInput4<-renderUI({
    if(input$choice==1)
    {
      numericInput("ratio", 
                   h5("Allocation ratio (treatment/placebo)"),value=1)
    }
  })
  
  output$conditionalInput5<-renderUI({
    if(input$choice==2){
    checkboxInput("checkbox", "Power Plot", value = FALSE)
    }
  })
  
  
  
  # f <- reactive({})       #update when input changes
  # f <- reactive(actionButton,{funct}) #update when action button from input is pressed
  # f <- reactive({         #update according to multiple input changes
  #   input$uielement1
  #   input$uielement2
  #   },{funct})
  
  ppr <- eventReactive({
    input$subBtn
  },{
    bet=input$effect
    sig0 = input$sig0
    rho = input$rho
    sig1=input$sig1
    sig01 = rho*sig0*sig1
    sig=input$sig 
    n_t=input$n_t
    r=input$droprate
    
    if (input$time_length==1){
    C=0
    h=input$h
    for (i in 1:(n_t-2))
    {
      ti<-seq(0,i*h,h)
      ti_mean<-mean(ti)
      ti2_mean<-mean(ti^2)
      Ci<- exp(-i*h*r)*(1-exp(-h*r))*(i+1)*(sig^2*ti2_mean+(i+1)*sig0^2*(ti2_mean-ti_mean^2))/(sig^4+(i+1)*sig^2*(sig1^2*ti2_mean+sig0^2+2*sig01*ti_mean)+(i+1)^2*(sig0^2*sig1^2-sig01^2)*(ti2_mean-ti_mean^2))
      C<-C+Ci
    }
    
    t<-seq(0,(n_t-1)*h,h)
    t_mean<-mean(t)
    t2_mean<-mean(t^2)
    
    C<-C+exp(-(n_t-1)*h*r)*n_t*(sig^2*t2_mean+n_t*sig0^2*(t2_mean-t_mean^2))/(sig^4+n_t*sig^2*(sig1^2*t2_mean+sig0^2+2*sig01*t_mean)+n_t^2*(sig0^2*sig1^2-sig01^2)*(t2_mean-t_mean^2))
    }
    
    if (input$time_length==2){
      C=0
      length_char=input$gap_length
      t=as.numeric(strsplit(length_char,split=",")[[1]])
      n_t=length(t)
      # t_mean=mean(t)
      # t2_mean=mean(t^2)
      for (i in 1:(n_t-2)){
        ti=t[1:(i+1)]
        ti_mean<-mean(ti)
        ti2_mean<-mean(ti^2)
        Ci<-(exp(-r*t[i+1])-exp(-r*t[i+2]))*(i+1)*(sig^2*ti2_mean+(i+1)*sig0^2*(ti2_mean-ti_mean^2))/(sig^4+(i+1)*sig^2*(sig1^2*ti2_mean+sig0^2+2*sig01*ti_mean)+(i+1)^2*(sig0^2*sig1^2-sig01^2)*(ti2_mean-ti_mean^2))
        C<-C+Ci
      }
      t_mean<-mean(t)
      t2_mean<-mean(t^2)
      
      C<-C+exp(-r*t[n_t])*n_t*(sig^2*t2_mean+n_t*sig0^2*(t2_mean-t_mean^2))/(sig^4+n_t*sig^2*(sig1^2*t2_mean+sig0^2+2*sig01*t_mean)+n_t^2*(sig0^2*sig1^2-sig01^2)*(t2_mean-t_mean^2))
  
    }
    
    if (input$level==1)
    {
      alpha=input$oneside
      if (input$choice==2)
      {
        m1=ifelse(is.null(input$n_tre),100,input$n_tre)
        m2=ifelse(is.null(input$n_pla),100,input$n_pla)
      
        var_btre<-(m1+m2)/(C*m1*m2)
        d<-bet/sqrt(var_btre)
        pr=1-pnorm(qnorm(1-alpha)-d)
        pr=round(pr,4)
      # pr<-pwr.t.test(d=d,n=(m1+m2),sig.level=alpha,type="one.sample",alternative="two.sided")
        # return(paste("The RCRM statistical power is ", as.numeric(format(round(pr, 3), nsmall = 2)), sep = ""))
        return(c(pr,C,bet,m1,m2))
        }
      
    
      else{
        ratio=input$ratio
        powerinput=ifelse(is.null(input$power),0.8,input$power)
      
        sp=(1+ratio)^2*(qnorm(1-alpha)+qnorm(powerinput))^2/(bet^2*C*ratio)
        spsize=ceiling(sp)
        # return(paste("The total sample size needed is ", spsize, sep = ""))
        return(c(spsize,ratio,powerinput))
      }
    }
    if (input$level==2)
    {
      alpha=input$twoside
      if (input$choice==2)
      {
        m1=ifelse(is.null(input$n_tre),100,input$n_tre)
        m2=ifelse(is.null(input$n_pla),100,input$n_pla)
        
        var_btre<-(m1+m2)/(C*m1*m2)
        d<-bet/sqrt(var_btre)
        pr=1-(pnorm(qnorm(1-alpha/2)-d)-pnorm(qnorm(alpha/2)-d))
        pr=round(pr,4)
        # pr<-pwr.t.test(d=d,n=(m1+m2),sig.level=alpha,type="one.sample",alternative="two.sided")
        # return(paste("The RCRM statistical power is ", as.numeric(format(round(pr, 3), nsmall = 2)), sep = ""))
        return(c(pr,C,bet,m1,m2))
          

        }
      
      
      else{
        ratio=input$ratio
        powerinput=ifelse(is.null(input$power),0.8,input$power)
        
        sp=(1+ratio)^2*(qnorm(1-alpha/2)+qnorm(powerinput))^2/(bet^2*C*ratio)
        spsize=ceiling(sp)
        # return(paste("The total sample size needed is ", spsize, sep = ""))
        return(c(spsize,ratio,powerinput))
      }
    }
})
    

  
  output$choice<-renderText({
    txt = ""
    if (input$choice==1){
      txt = "Given power, total sample size (treatment+placebo) is:"
    }else{
      txt = "Given sample size, power is:"
    }
    paste(txt,ppr()[1])
  })
  
  # output$power <- renderText({
  #   ppr()[1]
  # })
  
  
  
 
  
  powerplot<-eventReactive(input$subBtn,{
         m1=ppr()[4]
         m2=ppr()[5]
         C=ppr()[2]
         bet=ppr()[3]
         rate=m1/m2
        m2_seq=seq(1,500,0.1)
        m_seq=(1+rate)*m2_seq
        var_btre_seq<-(1+rate)/(C*rate*m2_seq)
        d_seq<-bet/sqrt(var_btre_seq)

        alpha=input$oneside
        pr_seq=1-pnorm(qnorm(1-alpha)-d_seq)
        labelname = "one-sided"
        
        if (input$level==2){
          alpha=input$twoside
          pr_seq=1-(pnorm(qnorm(1-alpha/2)-d_seq)-pnorm(qnorm(alpha/2)-d_seq))
          labelname = "two-sided"
          }
        
        ggplot2::ggplot(data =data.frame(m_seq,pr_seq),aes(x=m_seq,y=pr_seq,"l", colour='red'))+
          geom_point(show.legend=F)+
          geom_line(show.legend=F)+
          annotate("label",x = max(m_seq)-50, y = 0.125, label = paste(expression("\U03B1"),"=",alpha,"\n",labelname,"\ndropout=",input$droprate,sep=""))+
          labs(y='Power',
               x="Total sample size",
               #subtitle=paste(expression("\U03B1"),"=",alpha,"\none sided","\ndropout=",input$droprate,sep=""),
               title="Power Plot") +
          theme(plot.title=element_text(hjust=0.5))
      })
  
  output$rcrmTblDownBtn <- downloadHandler(
    filename = function() {
      "rcrmTbl.csv"
    },
    content = function(file) {
      write.csv(rst_table, file, row.names=FALSE)
    }
  )
  
  output$showTblDownBtn <- renderUI({
    if(input$subBtn > 0){
      downloadButton("rcrmTblDownBtn", "Download Table")
    }  
  })
  
  output$showplot<-renderUI({
    if (input$checkbox==TRUE){
      plotOutput("power_plot")
    }
  })
  
  output$showplotDownBtn <- renderUI({
    if (input$checkbox==TRUE){
      downloadButton("rcrm_downPlotBtn", label = "Download Plot")
    }
  })
  
  output$rcrm_downPlotBtn <- downloadHandler(
    filename = function() {"rcrmPlot.png"},
    content = function(file) {
      ggsave(plot = powerplot(), filename = file, device = "png", width = 11, height = 8.5)
    }
  )
  
  output$power_plot <- renderPlot({
    powerplot()
  })
  
  
  rcrm_newRow<-reactive({
    test_type=ifelse(input$level==1,"one-sided","two-sided")
    alpha=ifelse(input$level==1,input$oneside,input$twoside)
    power=ifelse(input$choice==1,input$power,ppr()[1])
    size=ifelse(input$choice==1,ppr()[1],(input$n_tre+input$n_pla))
    visit_time=ifelse(input$time_length==1,paste(seq(0,(input$n_t-1)*input$h,input$h),collapse=", "),input$gap_length)
    ratio=ifelse(input$choice==1,input$ratio,input$n_tre/input$n_pla)
    row=data.frame(test_type,alpha,power,size,ratio,input$effect,input$sig0,input$rho,input$sig1,input$sig,visit_time,input$droprate)

    return(row)
  })
  
  observeEvent(input$subBtn,{
    rst_table <<- rbind(rst_table,rcrm_newRow())
    colnames(rst_table)<-c("Test Type","\U03B1","Power","Total Sample Size","Allocation Ratio (treatment/placebo)",
                     "Effect Size",paste0("\U03C3","_","\U03B1"),"\U03C1",paste0("\U03C3","_","\U03B2"),"\U03C3","Assessment Time","Exponential Dropout Rate")
    output$rcrm_table <- DT::renderDataTable(
      datatable(rst_table,
                rownames = FALSE,
                extensions=c('Buttons','ColReorder'),
                options=list(scrollX=TRUE,dom='Bfrtip',pageLength = 10,
                             buttons=list('colvis'),
                                          # list(
                                          #   extend='csv',
                                          #   filename='resultsTbl'),
                                          # list(
                                          #   extend='pdf',
                                          #   pageSize = 'A4',
                                          #   orientation = 'landscape',
                                          #   filename='resultsTbl')),
                             colReorder = TRUE),
                filter="top"), server = FALSE
      )
  })
  
  

  
  # output$rcrmTbl <- DT::renderDataTable(
  #   datatable(get_rcrm_table(),
  #             rownames = FALSE,
  #             extensions=c('Buttons'),
  #             options=list(scrollX=TRUE, dom ='Bfrtip',
  #                          pageLength = 10,
  #                          buttons = list('colvis',
  #                                         list(extend='csf',filename ='results'))),
  #             filter = "top"),
  #   server = FALSE
  # )

  
  ###########Two-stage############
  tsm_rst_table<-data.frame(matrix(nrow=0,ncol=9))
  
  
  output$tsm_side1<-renderUI({
    if(input$tsm_level==1){
      numericInput("tsm_oneside", "",value=0.025)
    }
  })
  output$tsm_side2<-renderUI({
    if(input$tsm_level==2){
      numericInput("tsm_twoside", "",value=0.05)
    }
  })
  output$tsm_conditionalInput <- renderUI({
    if(input$tsm_choice==1){
      numericInput("tsm_power", 
                   h5("Expected power"), 
                   value = 0.8)
    } 
  })
  
  output$tsm_conditionalInput2<-renderUI({
    if(input$tsm_choice==2)
    {
      numericInput("tsm_n_tre", 
                   h5("Sample size of treatment group"), 
                   value = 100)
    }
  })
  
  output$tsm_conditionalInput3<-renderUI({
    if(input$tsm_choice==2)
    {
      numericInput("tsm_n_pla", 
                   h5("Sample size of placebo group"),value=100)
    }
  })
  
  output$tsm_conditionalInput4<-renderUI({
    if(input$tsm_choice==1)
    {
      numericInput("tsm_ratio", 
                   h5("Allocation ratio (treatment/placebo)"),value=1)
    }
  })
  
  output$tsm_conditionalInput5<-renderUI({
    if(input$tsm_choice==2){
      checkboxInput("tsm_checkbox", "Power Plot", value = FALSE)
    }
  })
  
  ppr_tsm<-eventReactive({
    input$tsm_subBtn
  },{
      bet_tsm=input$tsm_effect
      duration_tsm=input$tsm_duration
      M=input$tsm_nt
      sigb_tsm=input$tsm_sig1
      sigw_tsm=input$tsm_sig
      C_tsm=(sigw_tsm^2*12*(M-1)/(M*(M+1)*duration_tsm^2)+sigb_tsm^2)
      
      
      if (input$tsm_level==1){
        alpha_tsm=input$tsm_oneside
        if (input$tsm_choice==2){
          m1_tsm=ifelse(is.null(input$tsm_n_tre),100,input$tsm_n_tre)
          m2_tsm=ifelse(is.null(input$tsm_n_pla),100,input$tsm_n_pla)
          var_bet_tsm=C_tsm*(1/m1_tsm+1/m2_tsm)
          d_tsm=bet_tsm/sqrt(var_bet_tsm)
          
          pr_tsm=1-pnorm(qnorm(1-alpha_tsm)-d_tsm)
          pr_tsm=round(pr_tsm,4)
          return(c(pr_tsm,C_tsm,bet_tsm,m1_tsm,m2_tsm))
        }
        else{
          ratio_tsm=input$tsm_ratio
          powerinput_tsm=ifelse(is.null(input$tsm_power),0.8,input$tsm_power)
          
          sp_tsm=(1+ratio_tsm)^2*(qnorm(1-alpha_tsm)+qnorm(powerinput_tsm))^2*C_tsm/(bet_tsm^2*ratio_tsm)
          spsize_tsm=ceiling(sp_tsm)
          return(c(spsize_tsm,ratio_tsm,powerinput_tsm))
        }
      }
      
      if (input$tsm_level==2){
        alpha_tsm=input$tsm_twoside
        if (input$tsm_choice==2){
          m1_tsm=ifelse(is.null(input$tsm_n_tre),100,input$tsm_n_tre)
          m2_tsm=ifelse(is.null(input$tsm_n_pla),100,input$tsm_n_pla)
          var_bet_tsm=C_tsm*(1/m1_tsm+1/m2_tsm)
          d_tsm=bet_tsm/sqrt(var_bet_tsm)
          
          pr_tsm=1-(pnorm(qnorm(1-alpha_tsm/2)-d_tsm)-pnorm(qnorm(alpha_tsm/2)-d_tsm))
          pr_tsm=round(pr_tsm,4)
          
          return(c(pr_tsm,C_tsm,bet_tsm,m1_tsm,m2_tsm))
        }
        else{
          ratio_tsm=input$tsm_ratio
          powerinput_tsm=ifelse(is.null(input$tsm_power),0.8,input$tsm_power)
          
          sp_tsm=(1+ratio_tsm)^2*(qnorm(1-alpha_tsm/2)+qnorm(powerinput_tsm))^2*C_tsm/(bet_tsm^2*ratio_tsm)
          spsize_tsm=ceiling(sp_tsm)
          return(c(spsize_tsm,ratio_tsm,powerinput_tsm))
        }
      }
      
  })
  
  output$tsm_choice<-renderText({
    txt = ""
    if (input$tsm_choice==1){
      txt = "Given power, total sample size (treatment+placebo) is: "
    }else{
      txt = "Given sample size, power is: "
    }
    paste(txt,ppr_tsm()[1])
    
  })
  
  # output$tsm_power <- renderText({
  #   ppr_tsm()[1]
  # })
  
  
  powerplot_tsm<-eventReactive(input$tsm_subBtn,{
    m1_tsm=ppr_tsm()[4]
    m2_tsm=ppr_tsm()[5]
    C_tsm=ppr_tsm()[2]
    bet_tsm=ppr_tsm()[3]
    rate_tsm=m1_tsm/m2_tsm
    m2_tsm_seq=seq(1,500,0.1)
    m_tsm_seq=(1+rate_tsm)*m2_tsm_seq
    var_btre_tsm_seq<-C_tsm*(1+rate_tsm)/(rate_tsm*m2_tsm_seq)
    d_tsm_seq<-bet_tsm/sqrt(var_btre_tsm_seq)
    
    alpha_tsm=input$tsm_oneside
    pr_tsm_seq=1-pnorm(qnorm(1-alpha_tsm)-d_tsm_seq)
    labelname_tsm = "one-sided"
    
    if (input$tsm_level==2){
      alpha_tsm=input$tsm_twoside
      pr_tsm_seq=1-(pnorm(qnorm(1-alpha_tsm/2)-d_tsm_seq)-pnorm(qnorm(alpha_tsm/2)-d_tsm_seq))
      
      ggplot2::ggplot(data =data.frame(m_tsm_seq,pr_tsm_seq),aes(x=m_tsm_seq,y=pr_tsm_seq,"l", colour='red'))+
        geom_point(show.legend=F)+
        geom_line(show.legend=F)+
        annotate("label",x = max(m_tsm_seq)-50, y = 0.125, label = paste(expression("\U03B1"),"=",alpha_tsm,"\ntwo sided"))+
        labs(y='Power',
             x="Total sample size",
             title="Power Plot"
        ) +
        theme(plot.title=element_text(hjust=0.5))
      
    }
    
    ggplot2::ggplot(data =data.frame(m_tsm_seq,pr_tsm_seq),aes(x=m_tsm_seq,y=pr_tsm_seq,"l", colour='red'))+
      geom_point(show.legend=F)+
      geom_line(show.legend=F)+
      annotate("label",x = max(m_tsm_seq)-50, y = 0.125, label = paste(expression("\U03B1"),"=",alpha_tsm,"\n",labelname_tsm))+
      labs(y='Power',
           x="Total sample size",
           title="Power Plot"
      ) +
      theme(plot.title=element_text(hjust=0.5))
 
    
  })
  
  output$tsm_showplot<-renderUI({
    if (input$tsm_checkbox){
      plotOutput("tsm_power_plot")
    }
  })
  
  output$tsm_showplotDownBtn <- renderUI({
    if (input$tsm_checkbox){
      downloadButton("tsm_downPlotBtn", label = "Download Plot")
    }
  })
  
  output$tsm_downPlotBtn <- downloadHandler(
    filename = function() {"tsmPlot.png"},
    content = function(file) {
      ggsave(plot = powerplot_tsm(), filename = file, device = "png", width = 11, height = 8.5)
    }
  )
  
  output$tsm_power_plot <- renderPlot({
    powerplot_tsm()
  })
  
  output$tsmTblDownBtn <- downloadHandler(
    filename = function() {
      "tsmTbl.csv"
    },
    content = function(file) {
      write.csv(rst_table, file, row.names=FALSE)
    }
  )
  
  output$tsm_showTblDownBtn <- renderUI({
    if(input$tsm_subBtn > 0){
      downloadButton("tsmTblDownBtn", "Download Table")
    }  
  })
  
  tsm_newRow<-reactive({
    test_type=ifelse(input$tsm_level==1,"one-sided","two_sided")
    alpha=ifelse(input$tsm_level==1,input$tsm_oneside,input$tsm_twoside)
    power=ifelse(input$tsm_choice==1,input$tsm_power,ppr_tsm()[1])
    size=ifelse(input$tsm_choice==1,ppr_tsm()[1],(input$tsm_n_tre+input$tsm_n_pla))
    visit_time=paste(round(seq(0,input$tsm_duration,length.out = input$tsm_nt),3),collapse=", ")
    ratio=ifelse(input$tsm_choice==1,input$tsm_ratio,input$tsm_n_tre/input$tsm_n_pla)
    row=data.frame(test_type,alpha,power,size,ratio,input$tsm_effect,input$tsm_sig1,input$tsm_sig,visit_time)
    
    return(row)
  })
  
  observeEvent(input$tsm_subBtn,{
    tsm_rst_table <<- rbind(tsm_rst_table,tsm_newRow())
    colnames(tsm_rst_table)<-c("Test Type","\U03B1","Power","Total Sample Size","Allocation Ratio",
                           "Effect Size","Between suject SD","Within subject SD","Assessment Time")
    #output$rcrm_table <- DT::renderDataTable(rst_table)
    output$tsm_table <- DT::renderDataTable(
      datatable(tsm_rst_table,
                rownames = FALSE,
                extensions=c('Buttons','ColReorder'),
                options=list(scrollX=TRUE,dom='Bfrtip',pageLength = 10,
                             buttons=list('colvis'),
                             # list(
                             #   extend='csv',
                             #   filename='resultsTbl'),
                             # list(
                             #   extend='pdf',
                             #   pageSize = 'A4',
                             #   orientation = 'landscape',
                             #   filename='resultsTbl')),
                             colReorder = TRUE),
                filter="top"), server = FALSE
    )
  })

  
}


shinyApp(ui,server, options = list(launch.browser=T))


