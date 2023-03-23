################################################################################.

library(shiny)
library(MASS)
library(spatstat)
library(ggplot2)

################################################################################.
#random variate generation
rdirichlet<-function (n = 1, alpha) 
{
  Gam <- matrix(0, n, length(alpha))
  for (i in 1:length(alpha)) Gam[, i] <- rgamma(n, shape = alpha[i])
  Gam/rowSums(Gam)
}

#random variate generation (avoid NAs)
rdirichlet1<-function (n = 1, alpha) 
{
  rv=rdirichlet(n,alpha)
  if(n==1)
  {
    rv=t(as.matrix(rv))
  }
  #
  ids=which(is.na(rv[,1]))
  if(length(ids)>0)
  {
    for(i in 1:length(ids))
    {
      rv[ids[i],]=rep(0,1,length(alpha))
      rv[ids[i],sample(1:length(alpha),1)]=1
    }#i
  }#
  #
  return(rv)
}

#probability mass (or density) function of the multinomial-Dirichlet distribution
dMDir<-function(x,alpha,log=FALSE)
{
  x=as.vector(x)
  if(length(x)==length(alpha))
  {
    num=length(x)
    n=sum(x)
    #
    p1=lgamma(n+1)
    p2=lgamma(sum(alpha))
    p3=lgamma(n+sum(alpha))
    #
    p4=sum(lgamma(x+alpha))
    p5=sum(lgamma(x+1))
    p6=sum(lgamma(alpha))
    #
    f=p1+p2-p3+p4-p5-p6
    #
    if(log==FALSE)
    {
      return(exp(f))
    }else
    {
      return(f)
    }
    #
  }else
  {
    return(NA)
  }
}

#fitting of the SDM model
#mat is a species-site matrix
#if it is a vector, convert it into a matrix
fit<-function(mat)
{
  if(is.vector(mat))
  {
    mat=t(as.matrix(mat))
  }#
  ##################
  n=dim(mat)[2] #number of sites
  likelihood<-function(pars)
  {
    alpha=pars[1]
    alpha=rep(alpha,1,n)
    #
    one=function(v)
    {
      return(dMDir(v,alpha,log=TRUE))
    }#
    #
    all=apply(mat,1,one)
    return(-sum(all))
  }#
  ##################
  res=nlminb(runif(1),likelihood,lower=1e-20,upper=1e+20)
  #
  return(res)
}


#fitting of the independent NBD model
#mat is a species-site matrix
#if it is a vector, convert it into a matrix
fitNBD<-function(mat)
{
  x=as.vector(mat) #use vector data for fitting
  likelihood<-function(pars)
  {
    k=pars[1]
    mu=pars[2]
    #
    all=dnbinom(x,size=k,mu=mu,log=TRUE)
    return(-sum(all,na.rm=TRUE))
  }#
  ##################
  res=nlminb(runif(2),likelihood,lower=rep(1e-20,1,2),upper=rep(1e+20,1,2))
  #
  return(res)
}


#simulation of species distribution using Poisson cluster process
#generate a matrix with cnum and rnum, species presence in each cell was recorded as 1
#using bivariate probability function to generate offsprings
#Poisson distribution to generate number of parental points
#geometric distribution to generate number of offsprings for each parental point
#uniform distribution to obtain locations of parental points
#cnum determines number of columns, while rnum determines number of rows
#lambda is the parameter for Poisson distribution, prob is the parameter for geometric distribution
species.distribution<-function(delta,lambda,prob,cnum,rnum)
{
  pnum<-rpois(1,lambda)+1 # avoiding zeros, Poission distribution,lambda =  parameter intensity ρ
  offnum<-rgeom(pnum,prob)+1 # avoiding zeros,The Geometric Distribution, prob = 0.01
  # generate the locations of parental points
  x<-sample(cnum,pnum,replace=FALSE)
  y<-sample(rnum,pnum,replace=FALSE)
  #
  #now generating offspring locations
  sig<-matrix(c(delta,0,0,delta),2,2)
  pp<-vector()
  for(i in 1:pnum)
  {
    this<-mvrnorm(offnum[i],c(x[i],y[i]),sig)#bivariate normal distribution,delta link to σ2 
    this<-rbind(this,c(x[i],y[i])) #including parental point is tricky, otherwise sometimes this is NULL
    #exluding out-of-bound points
    idd<-which(this[,1]>cnum | this[,2]>rnum | this[,1]<1 | this[,2]<1)
    if(length(idd)!=0)
    {
      this<-this[-idd,]
    }
    pp<-rbind(pp,this)
  }
  #
  return(pp)
}

################################################################################.
ui <- fluidPage(
  theme = shiny::shinytheme("cerulean"),
  titlePanel(h1("Symmetric Dirichlet-multinomial (SDM) distribution simulator", style = "font-size: 24px;")),
  
  # 添加引文
  p("Please cite the following paper when using this application:",
    br(),
    strong("Liao, Z., Zhou, J., Shen, T. and Chen, Y. (2023)."),
    "Inferring single- and multi-species distributional aggregation using quadrat sampling.",
    em("Ecography, Volume(Issue), Page range."),
    "DOI: xxxxxxxxxx"),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      numericInput("Spnum", "Number of Species:", 2, width = "100%"),
      numericInput("cnum", "Simulation Range - X Coordinate:", 200, width = "100%"),
      numericInput("rnum", "Simulation Range - Y Coordinate:", 200, width = "100%"),
      numericInput("x", "Sample Area Size - Length:", 50, width = "100%"),
      numericInput("y", "Sample Area Size - Width:", 50, width = "100%"),
      actionButton("simulate", "Simulate", width = "100%"),
      hr()
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        tabPanel("Spatial Distribution", plotOutput("spatial_plot")),
        tabPanel("SDM Results", verbatimTextOutput("sdm_results")),
        tabPanel("Actual vs Predicted", plotOutput("bar_plot"))
      )
    )
  )
)


server <- function(input, output) {
  dat <- q1 <- q2 <- list()
  res1 <- NULL
  sdm_num <- NULL
  
  observeEvent(input$simulate, {
    # 模拟物种分布
    for (u in 1:input$Spnum) {
      dat[[u]] <- species.distribution(delta = sample(5:5000, 1), 
                                       lambda = sample(1:10, 1), 
                                       prob = sample(0.01:0.1, 1), 
                                       cnum = input$cnum, 
                                       rnum = input$rnum)
      if (is.null(dat[[u]])) {
        showNotification("Error generating species distribution. Please try again.", type = "error")
        return()
      }
      
      xy_ppp <- as.ppp(dat[[u]], W = c(0, input$cnum, 0, input$rnum))
      q1[[u]] <- quadratcount(xy_ppp, nx = input$cnum/input$x, ny = input$rnum/input$y)
      q2[[u]] <- matrix(q1[[u]], nrow = 1, ncol = input$cnum/input$x * input$rnum/input$y)
    }
    
    dat_final <- Reduce(rbind, q2)
    
    # 计算 SDM
    res1 <<- fit(dat_final)
    p <- rdirichlet1(1, alpha = rep(res1$par, length(dat_final)))
    sdm_num <<- rmultinom(1, size = sum(dat_final), prob = p)
    
    res2 <<- fitNBD(dat_final)
    nbd_num  = rnbinom(length(dat_final),size=res2$par[1],mu=res2$par[2])
    
    # 
    output$spatial_plot <- renderPlot({
      plot(dat[[1]], xlim = c(0, input$cnum), ylim = c(0, input$rnum), col = 2, pch = 1, cex = 1, xlab = "Coordinate X", ylab = "Coordinate Y")
      for (u in 1:input$Spnum) {
        points(dat[[u]], col = u, pch = u, cex = 1)
      }
    })
    
    output$sdm_results <- renderPrint({
      input$simulate_button
      isolate({
        res1_alpha <- res1$par
        res1_AIC <- res1$objective * 2 + 2 * 1
        res1_BIC <- res1$objective * 2 + log(length(dat_final)) * 1
        
        res2_k <- res2$par[1]
        res2_miu <- res2$par[2]
        res2_AIC <-  res2$objective * 2 + 2 * 1
        res2_BIC <- res2$objective * 2 + log(length(dat_final)) * 1
        
        cat("SDM Model Performance:\n")
        cat("α:", res1_alpha, "\n")
        cat("AIC:", res1_AIC, "\n")
        cat("BIC:", res1_BIC, "\n")
        cat("\n")
        
        cat("NBD Model Performance:\n")
        cat("k:", res2_k, "\n")
        cat("μ",res2_miu,"\n")
        cat("AIC:", res2_AIC, "\n")
        cat("BIC:", res2_BIC, "\n")
        cat("\n")
        
      })
    })
    
    output$bar_plot <- renderPlot({
      y_actual <- sort(dat_final)
      y_predicted_sdm <- sort(sdm_num)
      y_predicted_nbd <- sort(nbd_num)
      
      df <- data.frame(value = c(y_actual,y_predicted_sdm,y_predicted_nbd),
                       type = rep(c("Actual", "Predicted by SDM", "Predicted by NBD"), 
                                  each = length(y_actual)),
                       quadrat = rep(seq_along(y_actual), times = 3))
      ggplot(df, aes(x = quadrat, y = value, fill = type)) + 
        geom_bar(stat = "identity", position = "dodge") +
        scale_fill_manual(values = c("darkgreen", "#56B4E9","#E69F00")) +
        
        labs(title = "Actual vs predicted quadrat abundances",
             x = "Quadrats", y = "Abundance") +
        
        scale_x_continuous(breaks = seq_along(y_actual)) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        theme(axis.title = element_text(size = 12)) +
        
        guides(fill = guide_legend(title = NULL)) +
        theme(legend.title = element_text(size = 12), legend.text = element_text(size = 10))
    })
  })
}

shinyApp(ui = ui, server = server)
