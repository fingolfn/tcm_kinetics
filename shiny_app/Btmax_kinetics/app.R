



####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
##     #########  ###     ##   ########    ##     ##  ###########    #######              
##        ##      ####    ##   ##     ##   ##     ##      ##        ##                
##        ##      ## ##   ##   ##     ##   ##     ##      ##        ##                                     
##        ##      ##  ##  ##   ########    ##     ##      ##         ######          
##        ##      ##   ## ##   ##          ##     ##      ##              ##
##        ##      ##    ####   ##           ##   ##       ##              ##                
##     ########   ##     ###   ##            #####        ##         ######   
####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
library(shiny)
library(shinyWidgets)
library(deSolve)# loads the library

coop_list<-c("0.001","0.002","0.003","0.004","0.005","0.006","0.007","0.008","0.009",
             "0.01","0.02","0.03","0.04","0.05","0.06","0.07","0.08","0.09",
             "0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9",
             "1","2","3","4","5","6","7","8","9",
             "10","20","30","40","50","60","70","80","90",
             "100","200","300","400","500","600","700","800","900","1000")



numlist<-c("0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9",
           "1","2","3","4","5","6","7","8",
           "10","20","30","40","50","60","70","80","90",
           "100","200","300","400","500","600","700","800","900","1000",
           "1000","2000","3000","4000","5000","6000","7000","8000","9000","10000")

scale_list<-c("0.001","0.002","0.003","0.004","0.005","0.006","0.007","0.008","0.009",
              "0.01","0.02","0.03","0.04","0.05","0.06","0.07","0.08","0.09",
              "0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1","1","2","3","4","5","6","7","8","9","10")


slider_Kab <- sliderTextInput(inputId = "Kab",label = NULL,choices = numlist,
                              selected ="1000",grid = T,width="100%",pre="Keb = ",post=" nM",hide_min_max = T)

slider_At <- sliderTextInput(inputId = "At",label = NULL,choices = numlist,
                             selected ="1",grid = T,width="100%",pre="[E]t = ",post=" nM",hide_min_max = T)

slider_Kbc <- sliderTextInput(inputId = "Kbc",label = NULL,choices = numlist,
                              selected ="100",grid = T,width="100%",pre="Kbt = ",post=" nM",hide_min_max = T)

slider_Ct <- sliderTextInput(inputId = "Ct",label = NULL,choices = numlist,
                             selected ="0.1",grid = T,width="100%",pre="[T]t = ",post=" nM",hide_min_max = T)

slider_coop <- sliderTextInput(inputId = "alpha",label = NULL,choices = coop_list,
                               selected ="1",grid = T,width="100%",pre="Cooperativity = ",hide_min_max = T)

slider_scale <- sliderTextInput(inputId = "kcat",label = NULL,choices = scale_list,
                                selected ="0.1",grid = T,width="100%",pre="kcat = ",post=" s^-1",hide_min_max = T)


#slider_scale <- noUiSliderInput(
#  inputId = "scale",label=NULL, min = -3, max = 0,
#  value = 0, tooltips = F,
#  step = 0.1, orientation = "vertical",
#  color = "black", inline = F,
#  height = "330", width = "0",direction="rtl",margin=-100,
#)


####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
##     ########      ########      ########     ########   ########  ###     ##  #############                                                          
##     ##     ##     ##     ##     ##          ##          ##        ####    ##       ##                                      
##     ##     ##     ##    ##      ##          ##          ##        ## ##   ##       ##                                      
##     ########      #######       ########     #######    ########  ##  ##  ##       ##                                                 
##     ##            ##    ##      ##                 ##   ##        ##   ## ##       ##                                   
##     ##            ##     ##     ##                 ##   ##        ##    ####       ##                                      
##     ##            ##      ##    ########     #######    ########  ##     ###       ##                                                 
####################################################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################################################### 

ui <- fluidPage(

  
  
  chooseSliderSkin("Square"),
  
  setSliderColor(c("red","red","blue","blue","purple","black"), c(1,2,3,4,5,6)),
  
  
  
  sidebarLayout(
    sidebarPanel(width=4,
                 div(style = "display:inline",slider_Kab),
                 div(style = "display:inline",slider_At),
                 div(style = "display:inline",slider_Kbc),
                 div(style = "display:inline",slider_Ct),
                 div(style = "display:inline",slider_coop),
                 div(style = "display:inline",slider_scale),

                 downloadButton('downloadPlot', 'Download Plot'),div(style="display:inline",downloadButton('downloadData', 'Download Data'))
                 ),
    
    
    mainPanel(plotOutput("hist",height="600px"),width=8)
  )
  
)

#input & output are list-like objects
server <- function(input, output) {
  ####################################################################################################################################################################################################################################################################
  ####################################################################################################################################################################################################################################################################
  ##    ########    ##          ########    ##########                                     
  ##    ##     ##   ##         ##      ##       ##           
  ##    ##     ##   ##         ##      ##       ##                              
  ##    ########    ##         ##      ##       ##      
  ##    ##          ##         ##      ##       ##      
  ##    ##          ##         ##      ##       ##                         
  ##    ##          ########    ########        ##      
  ####################################################################################################################################################################################################################################################################
  ####################################################################################################################################################################################################################################################################
  

  output$hist <- renderPlot({
    
    At1<-as.numeric(input$At)
    
    Ct1<-as.numeric(input$Ct)
    Kab1<-as.numeric(input$Kab)
    Kbc1<-as.numeric(input$Kbc)
    
    alpha1=as.numeric(input$alpha)
    kcat1=as.numeric(input$kcat)
    
    ################################################################################################################
    #TIMESCALE APPROXIMATION:
    ################################################################################################################
    #time_scale_approx<-log(2)/(alpha1*kcat1*At1/(max(Kbc1,Kab1)))
    
    if (max(At1,Ct1)<(max(Kab1,Kbc1)/alpha1)){
      #PSEUDO-FIRST ORDER:   [X]t << Kweak
      time_scale_approx<-log(2)/(alpha1*kcat1*At1/(2*sqrt(Kbc1*Kab1)+Kbc1+Kab1))
    }else{
      if (At1>=Ct1){
        #SATD C:
        time_scale_approx<-log(2)/kcat1
      }
      if (At1<Ct1){
        #SATD A:
        time_scale_approx<-(Ct1/2)/(kcat1*At1)
      }
    }
    
    ################################################################################################################
    #SIMULATION:
    ################################################################################################################
    ## time sequence
    max=300000
    
    time <- seq(from=0, to=max, by = max/100)
    time_plus_approx<-sort(c(time,time_scale_approx))
    time<-time_plus_approx
    
    # parameters: a named vector
    parameters <- c(alpha=alpha1,At=At1,Kab=Kab1,Kbc=Kbc1,kcat=kcat1)
    
    # initial conditions: also a named vector
    state <- c(x = Ct1)
    
    logistic <- function(t, state, parameters){
      with(
        as.list(c(state, parameters)),{
          dx <- -kcat*(     x+At+ (sqrt(Kab)+sqrt(Kbc))^2/alpha   -sqrt( (    x+At+ (sqrt(Kab)+sqrt(Kbc))^2/alpha   )^2-4*x*At)   )   /2
          return(list(dx))
        }
      )
    }
    
    
    out <- ode(y = state, times = time, func = logistic, parms = parameters)
    
    
    ternaryMAXcoop_full_sim<-function(ode_soln,alpha=alpha1,At=At1,Kab=Kab1,Kbc=Kbc1,kcat=kcat1){
      #CALC each time course from At time course
      x=ode_soln[,2]
      
      ABC<-(     x+At+ (sqrt(Kab)+sqrt(Kbc))^2/alpha   -sqrt( (    x+At+ (sqrt(Kab)+sqrt(Kbc))^2/alpha   )^2-4*x*At)   )   /2
      B<-sqrt(Kab*Kbc)
      
      C<-(x-ABC)/(1+B/Kbc)
      A<-(At-ABC)/(1+B/Kab)
      
      BC<-B*C/Kbc
      AB<-A*B/Kab
      
      
      #ADD to matrix:
      
      all_sims<-cbind(ode_soln,ABC,AB,A,B,C,BC)
      colnames(all_sims)<-c("time","Ct","ABC","AB","A","B","C","BC")
      
      return(all_sims)
    }
    
    
    
    
    ################################################################################################################
    #ANNOTATED CURVE:
    ################################################################################################################
    full_maxCOOPsim<-ternaryMAXcoop_full_sim(out)
    
    par(mar=c(5,4,0,2))
    plot(full_maxCOOPsim[,c("time")]/3600,full_maxCOOPsim[,c("Ct")]/Ct1,ylim=c(0,1 ),
         type="l",lwd=2,xlab="hours",ylab=NA,cex.lab=1.5,cex.axis=1.5,las=1)
    
    
    
    #ADD TIMESCALE APPROXIMATION:
    t_half_idx<-which(time %in% time_scale_approx)
    t_half_MAG<-full_maxCOOPsim[t_half_idx,"Ct"]/Ct1
    t_half_TIME<-full_maxCOOPsim[t_half_idx,"time"]/3600
    
    
    points(t_half_TIME,t_half_MAG,pch=16,cex=1.5,col="black")
    
    segments(x0=t_half_TIME,y0=t_half_MAG,x1=t_half_TIME,y1=-1,lty=2,col="black")
    segments(x0=t_half_TIME,y0=t_half_MAG,x1=-10,y1=t_half_MAG,lty=2,col="black")
    
    
    
    
    

    
    
    
    
    ###################################################
    #T1/2 ALGORITHM:
    ###################################################
    if (max(At1,Ct1)<(max(Kab1,Kbc1)/alpha1)){
      ################################################################################################################
      #PSEUDO-FIRST ORDER:   [X]t << Kweak
      ################################################################################################################
      if (time_scale_approx<60){
        ###################
        #Seconds:
        ###################
        time_format<-round(time_scale_approx,1)
        legend("top",legend=bquote(t[1/2]  %~~% frac(ln(2),k[cat] * "[E]"[t])  %*% frac((sqrt(K[EB])+sqrt(K[BT]))^2,alpha) ~ "=" ~ .(time_format) ~ "seconds"),
               bty="n",cex=1.4)
        
      }
      if (time_scale_approx>=60 & time_scale_approx<3600){
        ###################
        #Minutes:
        ###################
        time_format<-round(time_scale_approx/60,1)
        legend("top",legend=bquote(t[1/2]  %~~% frac(ln(2),k[cat] * "[E]"[t])  %*% frac((sqrt(K[EB])+sqrt(K[BT]))^2,alpha) ~ "=" ~ .(time_format) ~ "minutes"),
               bty="n",cex=1.4)
      }
      if (time_scale_approx>=3600 & time_scale_approx<86400){
        ###################
        #Hours:
        ###################
        time_format<-round(time_scale_approx/3600,1)
        legend("top",legend=bquote(t[1/2]  %~~% frac(ln(2),k[cat] * "[E]"[t])  %*% frac((sqrt(K[EB])+sqrt(K[BT]))^2,alpha) ~ "=" ~ .(time_format) ~ "hours"),
               bty="n",cex=1.4)
      }
      if (time_scale_approx>86400){
        ###################
        #Hours:
        ###################
        time_format<-round(time_scale_approx/86400,1)
        legend("top",legend=bquote(t[1/2]  %~~% frac(ln(2),k[cat] * "[E]"[t])  %*% frac((sqrt(K[EB])+sqrt(K[BT]))^2,alpha) ~ "=" ~ .(time_format) ~ "days"),
               bty="n",cex=1.4)
      }
    }else{
      if (At1>=Ct1){
        ################################################################################################################
        #SATD C:
        ################################################################################################################
        if (time_scale_approx<60){
          ###################
          #Seconds:
          ###################
          time_format<-round(time_scale_approx,1)
          legend("top",legend=bquote(t[1/2]  %~~% frac(ln(2),k[cat]) ~ "=" ~ .(time_format) ~ "seconds"),
                 bty="n",cex=1.4)
          
        }
        if (time_scale_approx>=60 & time_scale_approx<3600){
          ###################
          #Minutes:
          ###################
          time_format<-round(time_scale_approx/60,1)
          legend("top",legend=bquote(t[1/2]  %~~% frac(ln(2),k[cat]) ~ "=" ~ .(time_format) ~ "minutes"),
                 bty="n",cex=1.4)
        }
        if (time_scale_approx>=3600 & time_scale_approx<86400){
          ###################
          #Hours:
          ###################
          time_format<-round(time_scale_approx/3600,1)
          legend("top",legend=bquote(t[1/2]  %~~% frac(ln(2),k[cat])  ~ "=" ~ .(time_format) ~ "hours"),
                 bty="n",cex=1.4)
        }
        if (time_scale_approx>86400){
          ###################
          #Hours:
          ###################
          time_format<-round(time_scale_approx/86400,1)
          legend("top",legend=bquote(t[1/2]  %~~% frac(ln(2),k[cat])  ~ "=" ~ .(time_format) ~ "days"),
                 bty="n",cex=1.4)
        }
      }
      if (At1<Ct1){
        ################################################################################################################
        #SATD A:
        ################################################################################################################
        if (time_scale_approx<60){
          ###################
          #Seconds:
          ###################
          time_format<-round(time_scale_approx,1)
          legend("top",legend=bquote(t[1/2]  %~~% frac("[T]"[t]/2,k[cat]* "[E]"[t]) ~ "=" ~ .(time_format) ~ "seconds"),
                 bty="n",cex=1.4)
          
        }
        if (time_scale_approx>=60 & time_scale_approx<3600){
          ###################
          #Minutes:
          ###################
          time_format<-round(time_scale_approx/60,1)
          legend("top",legend=bquote(t[1/2]  %~~% frac("[T]"[t]/2,k[cat]* "[E]"[t]) ~ "=" ~ .(time_format) ~ "minutes"),
                 bty="n",cex=1.4)
        }
        if (time_scale_approx>=3600 & time_scale_approx<86400){
          ###################
          #Hours:
          ###################
          time_format<-round(time_scale_approx/3600,1)
          legend("top",legend=bquote(t[1/2]  %~~% frac("[T]"[t]/2,k[cat]* "[E]"[t])  ~ "=" ~ .(time_format) ~ "hours"),
                 bty="n",cex=1.4)
        }
        if (time_scale_approx>86400){
          ###################
          #Hours:
          ###################
          time_format<-round(time_scale_approx/86400,1)
          legend("top",legend=bquote(t[1/2]  %~~% frac("[T]"[t]/2,k[cat]* "[E]"[t])  ~ "=" ~ .(time_format) ~ "days"),
                 bty="n",cex=1.4)
        }
      }
    }



    
    
    
    
    
    
    
    
    
  })
  ####################################################################################################################################################################################################################################################################
  ####################################################################################################################################################################################################################################################################
  ##    ########    ##          ########    ##########                                     
  ##    ##     ##   ##         ##      ##       ##           
  ##    ##     ##   ##         ##      ##       ##                              
  ##    ########    ##         ##      ##       ##      
  ##    ##          ##         ##      ##       ##      
  ##    ##          ##         ##      ##       ##                         
  ##    ##          ########    ########        ##      
  ####################################################################################################################################################################################################################################################################
  ####################################################################################################################################################################################################################################################################
  
  
  plotInput <- function(){
    
    At1<-as.numeric(input$At)
    
    Ct1<-as.numeric(input$Ct)
    Kab1<-as.numeric(input$Kab)
    Kbc1<-as.numeric(input$Kbc)
    
    alpha1=as.numeric(input$alpha)
    kcat1=as.numeric(input$kcat)
    
    ################################################################################################################
    #TIMESCALE APPROXIMATION:
    ################################################################################################################
    #time_scale_approx<-log(2)/(alpha1*kcat1*At1/(max(Kbc1,Kab1)))
    
    if (max(At1,Ct1)<(max(Kab1,Kbc1)/alpha1)){
      #PSEUDO-FIRST ORDER:   [X]t << Kweak
      time_scale_approx<-log(2)/(alpha1*kcat1*At1/(2*sqrt(Kbc1*Kab1)+Kbc1+Kab1))
    }else{
      if (At1>=Ct1){
        #SATD C:
        time_scale_approx<-log(2)/kcat1
      }
      if (At1<Ct1){
        #SATD A:
        time_scale_approx<-(Ct1/2)/(kcat1*At1)
      }
    }
    
    ################################################################################################################
    #SIMULATION:
    ################################################################################################################
    ## time sequence
    max=300000
    
    time <- seq(from=0, to=max, by = max/100)
    time_plus_approx<-sort(c(time,time_scale_approx))
    time<-time_plus_approx
    
    # parameters: a named vector
    parameters <- c(alpha=alpha1,At=At1,Kab=Kab1,Kbc=Kbc1,kcat=kcat1)
    
    # initial conditions: also a named vector
    state <- c(x = Ct1)
    
    logistic <- function(t, state, parameters){
      with(
        as.list(c(state, parameters)),{
          dx <- -kcat*(     x+At+ (sqrt(Kab)+sqrt(Kbc))^2/alpha   -sqrt( (    x+At+ (sqrt(Kab)+sqrt(Kbc))^2/alpha   )^2-4*x*At)   )   /2
          return(list(dx))
        }
      )
    }
    
    
    out <- ode(y = state, times = time, func = logistic, parms = parameters)
    
    
    ternaryMAXcoop_full_sim<-function(ode_soln,alpha=alpha1,At=At1,Kab=Kab1,Kbc=Kbc1,kcat=kcat1){
      #CALC each time course from At time course
      x=ode_soln[,2]
      
      ABC<-(     x+At+ (sqrt(Kab)+sqrt(Kbc))^2/alpha   -sqrt( (    x+At+ (sqrt(Kab)+sqrt(Kbc))^2/alpha   )^2-4*x*At)   )   /2
      B<-sqrt(Kab*Kbc)
      
      C<-(x-ABC)/(1+B/Kbc)
      A<-(At-ABC)/(1+B/Kab)
      
      BC<-B*C/Kbc
      AB<-A*B/Kab
      
      
      #ADD to matrix:
      
      all_sims<-cbind(ode_soln,ABC,AB,A,B,C,BC)
      colnames(all_sims)<-c("time","Ct","ABC","AB","A","B","C","BC")
      
      return(all_sims)
    }
    
    
    
    
    ################################################################################################################
    #ANNOTATED CURVE:
    ################################################################################################################
    full_maxCOOPsim<-ternaryMAXcoop_full_sim(out)
    
    par(mar=c(5,4,0,2))
    plot(full_maxCOOPsim[,c("time")]/3600,full_maxCOOPsim[,c("Ct")]/Ct1,ylim=c(0,1 ),
         type="l",lwd=2,xlab="hours",ylab=NA,cex.lab=1.5,cex.axis=1.5,las=1)
    
    
    
    #ADD TIMESCALE APPROXIMATION:
    t_half_idx<-which(time %in% time_scale_approx)
    t_half_MAG<-full_maxCOOPsim[t_half_idx,"Ct"]/Ct1
    t_half_TIME<-full_maxCOOPsim[t_half_idx,"time"]/3600
    
    
    points(t_half_TIME,t_half_MAG,pch=16,cex=1.5,col="black")
    
    segments(x0=t_half_TIME,y0=t_half_MAG,x1=t_half_TIME,y1=-1,lty=2,col="black")
    segments(x0=t_half_TIME,y0=t_half_MAG,x1=-10,y1=t_half_MAG,lty=2,col="black")
    
    
    
    
    
    
    
    
    
    
    ###################################################
    #T1/2 ALGORITHM:
    ###################################################
    if (max(At1,Ct1)<(max(Kab1,Kbc1)/alpha1)){
      ################################################################################################################
      #PSEUDO-FIRST ORDER:   [X]t << Kweak
      ################################################################################################################
      if (time_scale_approx<60){
        ###################
        #Seconds:
        ###################
        time_format<-round(time_scale_approx,1)
        legend("top",legend=bquote(t[1/2]  %~~% frac(ln(2),k[cat] * "[E]"[t])  %*% frac((sqrt(K[EB])+sqrt(K[BT]))^2,alpha) ~ "=" ~ .(time_format) ~ "seconds"),
               bty="n",cex=1.4)
        
      }
      if (time_scale_approx>=60 & time_scale_approx<3600){
        ###################
        #Minutes:
        ###################
        time_format<-round(time_scale_approx/60,1)
        legend("top",legend=bquote(t[1/2]  %~~% frac(ln(2),k[cat] * "[E]"[t])  %*% frac((sqrt(K[EB])+sqrt(K[BT]))^2,alpha) ~ "=" ~ .(time_format) ~ "minutes"),
               bty="n",cex=1.4)
      }
      if (time_scale_approx>=3600 & time_scale_approx<86400){
        ###################
        #Hours:
        ###################
        time_format<-round(time_scale_approx/3600,1)
        legend("top",legend=bquote(t[1/2]  %~~% frac(ln(2),k[cat] * "[E]"[t])  %*% frac((sqrt(K[EB])+sqrt(K[BT]))^2,alpha) ~ "=" ~ .(time_format) ~ "hours"),
               bty="n",cex=1.4)
      }
      if (time_scale_approx>86400){
        ###################
        #Hours:
        ###################
        time_format<-round(time_scale_approx/86400,1)
        legend("top",legend=bquote(t[1/2]  %~~% frac(ln(2),k[cat] * "[E]"[t])  %*% frac((sqrt(K[EB])+sqrt(K[BT]))^2,alpha) ~ "=" ~ .(time_format) ~ "days"),
               bty="n",cex=1.4)
      }
    }else{
      if (At1>=Ct1){
        ################################################################################################################
        #SATD C:
        ################################################################################################################
        if (time_scale_approx<60){
          ###################
          #Seconds:
          ###################
          time_format<-round(time_scale_approx,1)
          legend("top",legend=bquote(t[1/2]  %~~% frac(ln(2),k[cat]) ~ "=" ~ .(time_format) ~ "seconds"),
                 bty="n",cex=1.4)
          
        }
        if (time_scale_approx>=60 & time_scale_approx<3600){
          ###################
          #Minutes:
          ###################
          time_format<-round(time_scale_approx/60,1)
          legend("top",legend=bquote(t[1/2]  %~~% frac(ln(2),k[cat]) ~ "=" ~ .(time_format) ~ "minutes"),
                 bty="n",cex=1.4)
        }
        if (time_scale_approx>=3600 & time_scale_approx<86400){
          ###################
          #Hours:
          ###################
          time_format<-round(time_scale_approx/3600,1)
          legend("top",legend=bquote(t[1/2]  %~~% frac(ln(2),k[cat])  ~ "=" ~ .(time_format) ~ "hours"),
                 bty="n",cex=1.4)
        }
        if (time_scale_approx>86400){
          ###################
          #Hours:
          ###################
          time_format<-round(time_scale_approx/86400,1)
          legend("top",legend=bquote(t[1/2]  %~~% frac(ln(2),k[cat])  ~ "=" ~ .(time_format) ~ "days"),
                 bty="n",cex=1.4)
        }
      }
      if (At1<Ct1){
        ################################################################################################################
        #SATD A:
        ################################################################################################################
        if (time_scale_approx<60){
          ###################
          #Seconds:
          ###################
          time_format<-round(time_scale_approx,1)
          legend("top",legend=bquote(t[1/2]  %~~% frac("[T]"[t]/2,k[cat]* "[E]"[t]) ~ "=" ~ .(time_format) ~ "seconds"),
                 bty="n",cex=1.4)
          
        }
        if (time_scale_approx>=60 & time_scale_approx<3600){
          ###################
          #Minutes:
          ###################
          time_format<-round(time_scale_approx/60,1)
          legend("top",legend=bquote(t[1/2]  %~~% frac("[T]"[t]/2,k[cat]* "[E]"[t]) ~ "=" ~ .(time_format) ~ "minutes"),
                 bty="n",cex=1.4)
        }
        if (time_scale_approx>=3600 & time_scale_approx<86400){
          ###################
          #Hours:
          ###################
          time_format<-round(time_scale_approx/3600,1)
          legend("top",legend=bquote(t[1/2]  %~~% frac("[T]"[t]/2,k[cat]* "[E]"[t])  ~ "=" ~ .(time_format) ~ "hours"),
                 bty="n",cex=1.4)
        }
        if (time_scale_approx>86400){
          ###################
          #Hours:
          ###################
          time_format<-round(time_scale_approx/86400,1)
          legend("top",legend=bquote(t[1/2]  %~~% frac("[T]"[t]/2,k[cat]* "[E]"[t])  ~ "=" ~ .(time_format) ~ "days"),
                 bty="n",cex=1.4)
        }
      }
    }
    
    
    
  }
  
  
  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste("plot-", Sys.Date(), ".pdf", sep="")
    },
    content = function(file) {
      pdf(file,width=10,height=8)
      plotInput()
      dev.off()
    },
    contentType = "image/pdf") 
  
  
  ####################################################################################################################################################################################################################################################################
  ####################################################################################################################################################################################################################################################################
  ##      #########           ##       ############      ##
  ##      ##      ##         ####           ##          ####
  ##      ##       ##       ##  ##          ##         ##  ##
  ##      ##        ##     ##    ##         ##        ##    ##
  ##      ##       ##     ##########        ##       ##########
  ##      ##      ##     ##        ##       ##      ##        ##
  ##      #########     ##          ##      ##     ##          ##
  ####################################################################################################################################################################################################################################################################
  #################################################################################################################################################################################################################################################################### 
  
  
  outputData <- function(){
    
    
    At1<-as.numeric(input$At)
    
    Ct1<-as.numeric(input$Ct)
    Kab1<-as.numeric(input$Kab)
    Kbc1<-as.numeric(input$Kbc)
    
    alpha1=as.numeric(input$alpha)
    kcat1=as.numeric(input$kcat)
    
    ################################################################################################################
    #TIMESCALE APPROXIMATION:
    ################################################################################################################
    #time_scale_approx<-log(2)/(alpha1*kcat1*At1/(max(Kbc1,Kab1)))
    
    if (max(At1,Ct1)<(max(Kab1,Kbc1)/alpha1)){
      #PSEUDO-FIRST ORDER:   [X]t << Kweak
      time_scale_approx<-log(2)/(alpha1*kcat1*At1/(2*sqrt(Kbc1*Kab1)+Kbc1+Kab1))
    }else{
      if (At1>=Ct1){
        #SATD C:
        time_scale_approx<-log(2)/kcat1
      }
      if (At1<Ct1){
        #SATD A:
        time_scale_approx<-(Ct1/2)/(kcat1*At1)
      }
    }
    
    ################################################################################################################
    #SIMULATION:
    ################################################################################################################
    ## time sequence
    max=300000
    
    time <- seq(from=0, to=max, by = max/100)
    time_plus_approx<-sort(c(time,time_scale_approx))
    time<-time_plus_approx
    
    # parameters: a named vector
    parameters <- c(alpha=alpha1,At=At1,Kab=Kab1,Kbc=Kbc1,kcat=kcat1)
    
    # initial conditions: also a named vector
    state <- c(x = Ct1)
    
    logistic <- function(t, state, parameters){
      with(
        as.list(c(state, parameters)),{
          dx <- -kcat*(     x+At+ (sqrt(Kab)+sqrt(Kbc))^2/alpha   -sqrt( (    x+At+ (sqrt(Kab)+sqrt(Kbc))^2/alpha   )^2-4*x*At)   )   /2
          return(list(dx))
        }
      )
    }
    
    
    out <- ode(y = state, times = time, func = logistic, parms = parameters)
    
    
    ternaryMAXcoop_full_sim<-function(ode_soln,alpha=alpha1,At=At1,Kab=Kab1,Kbc=Kbc1,kcat=kcat1){
      #CALC each time course from At time course
      x=ode_soln[,2]
      
      ABC<-(     x+At+ (sqrt(Kab)+sqrt(Kbc))^2/alpha   -sqrt( (    x+At+ (sqrt(Kab)+sqrt(Kbc))^2/alpha   )^2-4*x*At)   )   /2
      B<-sqrt(Kab*Kbc)
      
      C<-(x-ABC)/(1+B/Kbc)
      A<-(At-ABC)/(1+B/Kab)
      
      BC<-B*C/Kbc
      AB<-A*B/Kab
      
      
      #ADD to matrix:
      
      all_sims<-cbind(ode_soln,ABC,AB,A,B,C,BC)
      colnames(all_sims)<-c("time","Ct","ABC","AB","A","B","C","BC")
      
      return(all_sims)
    }
    
    
    
    
    ################################################################################################################
    #ANNOTATED CURVE:
    ################################################################################################################
    full_maxCOOPsim<-ternaryMAXcoop_full_sim(out)
    

    

    
    return(full_maxCOOPsim)
  }
  
  
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("data-", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      data<-outputData()
      write.csv(data, file)
    }
  )
  
}

shinyApp(ui = ui, server = server)

