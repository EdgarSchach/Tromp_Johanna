####Loading sources and packages####

source("./gmMLA-master/R/mlaClasses.R")
source("./gmMLA-master/R/mlaModals.R")
source("./gmMLA-master/R/mlaSubsetting.R")
source("./gmMLA-master/R/mlaCurves.R")
source("./gmMLA-master/R/mlaXMLSession.R")
library("plyr")
library("dplyr")
library("reshape2")
library("ggplot2")
library("xlsx")
library("broom")
library("purrr")
library("pls")
library("ggfortify")

# load data 

load("./samples_all_map.RData")
nms <- names(samples)
samples <- samples[c(6:11)]
names(samples) <- nms[c(1:6)]

####calculating kerneldensity####

#Loading relevant Data:

x = samples$SC0_C1_grouped

relevants = function(x) data.frame(aspect_ratio = x$geomProps$MinRectMinLength/x$geomProps$MinRectMaxLength,
                                   density = x$geomProp$Mass/x$geomProp$ParticleArea/100,
                                   ECD = sqrt(x$geomProp$ParticleArea/pi)*2,
                                   mass = x$geomProp$Mass,
                                   mass_glass = x$minMassComp[,"A-Glass"]
                                      )
relevant.data <- lapply(samples,relevants)
relevant.data <- lapply(relevant.data,na.omit)

part.dat <- relevant.data$SC0_C1_grouped
count_function <- function(part.dat,min.ps=1,max){
part.dat$size_class <- cut(part.dat$ECD, breaks = c(0,1,10,1000),labels = c("0-1","1-10",">10"))
part.dat %>% group_by(size_class) %>% count()
part.dat %>% filter(aspect_ratio < 0.25) %>% count()
}

pre_processing <- function(part.data,max.ps = 10,min.ps = 1){
  n_unfiltered <- print(paste(nrow(part.data),"unfiltered particles"))
  res <- filter(part.data, ECD <= max.ps & ECD >= min.ps  & aspect_ratio > 0.25)  # delete particles bigger x µm with standard 10 µm
  n_filtered <- print(paste(nrow(res),"filtered particles"))
  return(res)
}


relevant.data <- lapply(relevant.data,pre_processing)


#Defining the ranges of important properties:

rf <- range(sapply(relevant.data,function(x)range(x$aspect_ratio)))
rs <- range(sapply(relevant.data,function(x)range(x$ECD)))

#Defining "smoothing parameter:

bws = c(0.05, 2)

#Calculating 2d Kernel Distribution:

mykde = function(x) kde2dWeighted(x$aspect_ratio, x$ECD, n=50, lims=c(rf,rs), w=x$mass_glass, h=bws)
kdes = lapply(relevant.data, mykde)

#Plotting the 2d-Kernel Distribution for Feed,Concentrate and Tailings:

image(kdes[[1]], main="Concentrate", xlab = "aspect ratio", ylab = "ECD in µm")

#defining weights of experiment: 

sample.weights = c(C1 = 0.6,C2 = 0.55,C3 = 0.89,C4 = 0.49,C5 = 0.24,T = 94.83) # in g

#Calculating 2d-Tromp-Curve:


twoDTromp = kdes[[1]] #defining dummy sample

twoDTromp$z = (sample.weights["C1"]*kdes[[1]]$z+sample.weights["C2"]*kdes[[2]]$z+sample.weights["C3"]*kdes[[3]]$z+sample.weights["C4"]*kdes[[4]]$z+sample.weights["C5"]*kdes[[5]]$z)/(sample.weights["C1"]*kdes[[1]]$z+sample.weights["C2"]*kdes[[2]]$z+sample.weights["C3"]*kdes[[3]]$z+sample.weights["C4"]*kdes[[4]]$z+sample.weights["C5"]*kdes[[5]]$z+sample.weights["T"]*kdes[[6]]$z)
image(twoDTromp, main="Prob to concentrate", xlab = "density in g/cm³", ylab = "ECD in µm" )

# Transform data for ggplot:

ggdat  <- as.data.frame(expand.grid(twoDTromp$x,twoDTromp$y))
ggdat[,"Recovery"] <- c(twoDTromp$z)

g <- ggplot()
g <- g + geom_raster(mapping = aes(x =Var2, y=Var1, fill=Recovery), data = ggdat, interpolate = FALSE)
g <- g + labs(x = "ECD in µm", y = "Aspect Ratio")
#g <- g + geom_point(data = relevant.data$SC0_C1_grouped,aes(x = ECD,y = aspect_ratio))
g <- g +scale_fill_gradientn(name=expression(paste(italic(r)," in wt.-%")),
                             limits=c(0,1),
                             breaks=seq(0, 1, by=0.2),
                             values = c(0,0.01,0.3,1),#!adjust to maximum
                             colours=c("white", "darkgreen","yellow","red"))
g <- g + scale_x_continuous(expand = c(0,0))
g <- g + scale_y_continuous(expand = c(0,0), limits = c(0.25,1))
print(g)

# prob distribution for concentrate für Konzentrat: 

kde.conc <- as.data.frame(expand.grid(twoDTromp$x,twoDTromp$y))
kde.conc$Prob <- c(sample.weights["C1"]*kdes[[1]]$z+sample.weights["C2"]*kdes[[2]]$z+sample.weights["C3"]*kdes[[3]]$z+sample.weights["C4"]*kdes[[4]]$z+sample.weights["C5"]*kdes[[5]]$z)/sum(sample.weights[c(1:5)])


g <- ggplot()
g <- g + geom_raster(mapping = aes(x =Var2, y=Var1, fill=Prob), data = kde.conc, interpolate = FALSE)
g <- g + labs(x = "ECD in µm", y = "Aspect Ratio")
#g <- g + geom_point(data = relevant.data$SC0_C1_grouped,aes(x = ECD,y = aspect_ratio))
g <- g +scale_fill_gradientn(name=expression(paste(italic(r)," in wt.-%")),
                             limits=c(0,5),
                             breaks=seq(0, 5, by=1),
                             values = c(0,0.01,0.3,1),#!adjust to maximum
                             colours=c("white", "darkgreen","yellow","red"))
g <- g + scale_x_continuous(expand = c(0,0))
g <- g + scale_y_continuous(expand = c(0,0), limits = c(0.25,1))
print(g)

# prob distribution for tailings: 

kde.tail <- as.data.frame(expand.grid(twoDTromp$x,twoDTromp$y))
kde.tail$Prob <- c(kdes[[6]]$z)


g <- ggplot()
g <- g + geom_raster(mapping = aes(x =Var2, y=Var1, fill=Prob), data = kde.tail, interpolate = FALSE)
g <- g + labs(x = "ECD in µm", y = "Aspect Ratio")
#g <- g + geom_point(data = relevant.data$SC0_C1_grouped,aes(x = ECD,y = aspect_ratio))
g <- g +scale_fill_gradientn(name=expression(paste(italic(r)," in wt.-%")),
                             limits=c(0,5),
                             breaks=seq(0, 5, by=1),
                             values = c(0,0.01,0.3,1),#!adjust to maximum
                             colours=c("white", "darkgreen","yellow","red"))
g <- g + scale_x_continuous(expand = c(0,0))
g <- g + scale_y_continuous(expand = c(0,0), limits = c(0.25,1))
print(g)

# pron distribution for feed recalculated:

kde.tail <- as.data.frame(expand.grid(twoDTromp$x,twoDTromp$y))
kde.tail$Prob <-  c(sample.weights["C1"]*kdes[[1]]$z+sample.weights["C2"]*kdes[[2]]$z+sample.weights["C3"]*kdes[[3]]$z+sample.weights["C4"]*kdes[[4]]$z+sample.weights["C5"]*kdes[[5]]$z+sample.weights["T"]*kdes[[6]]$z)/sum(sample.weights)


g <- ggplot()
g <- g + geom_raster(mapping = aes(x =Var2, y=Var1, fill=Prob), data = kde.tail, interpolate = FALSE)
g <- g + labs(x = "ECD in µm", y = "Aspect Ratio")
#g <- g + geom_point(data = relevant.data$SC0_C1_grouped,aes(x = ECD,y = aspect_ratio))
g <- g +scale_fill_gradientn(name=expression(paste(italic(r)," in wt.-%")),
                             limits=c(0,5),
                             breaks=seq(0, 5, by=1),
                             values = c(0,0.01,0.3,1),#!adjust to maximum
                             colours=c("white", "darkgreen","yellow","red"))
g <- g + scale_x_continuous(expand = c(0,0))
g <- g + scale_y_continuous(expand = c(0,0), limits = c(0.25,1))
print(g)

#Eigenschaftenverteilung aus Probe: 

part.dat <- bind_rows(relevant.data$SC0_C1_grouped,relevant.data$SC0_T_grouped,.id = "Sample")

g <- ggplot(part.dat,aes(x = ECD, y = aspect_ratio, fill = "Sample")) +
  geom_density_2d()
print(g)

