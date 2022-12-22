#WARNING: this work is in progress. The code and method has not yet been checked by statisticians... use at your own risk! 
#Please check the consistency of your results if you use it and stay informed of the further development of this script.

#Written by Marc Thomas marcthomas1@hotmail.fr
#last modification: 22th December 2022

#This script aims to test spatial associations of objects in three dimensions by using the Ripley's K function.
#loading packages
library(spatstat)
library(dplyr)
library(bioimagetools)
library(ggplot2)
library(cowplot)

#setRepositories(ind=c(1,2))
#install.packages("bioimagetools")

#function published on its github by volkerschmid (https://github.com/bioimaginggroup/bioimagetools/blob/master/R/Kcross3D.R)
#you don't need to edit anything in the following section
K.cross.3D.envelope<-function(X,Y,Z,X2,Y2,Z2,range,bw.int,psz,width,parallel=FALSE,N=100,use.intensity=TRUE)
{
  K.envelope<-c()
  if(use.intensity)intensity.red<-intensity3D(X+range,Y+range,Z+range,bw=bw.int,psz=psz)
  if(!use.intensity){intensity.red<-NULL;intensity.green<-NULL}
  shuffle<-function(i)
  {
    XX<-X2+runif(1,0,2*range)
    YY<-Y2+runif(1,0,2*range)
    ZZ<-Z2+runif(1,0,2*range)
    if(use.intensity)intensity.green<-intensity3D(XX,YY,ZZ,bw=bw.int,psz=psz)
    K.temp<-K.cross.3D(X,Y,Z,XX,YY,ZZ,psz=psz,width=width,intensity=intensity.red,intensity2=intensity.green,parallel=parallel)
    return(K.temp$y)
  }
  if(!parallel)envelope<-lapply(1:N,shuffle)
  if(parallel)envelope<-parallel::mclapply(1:N,shuffle)
  envelope<-unlist(envelope)
  envelope<-array(envelope,c(length(envelope)/N,N))
  
  if(use.intensity)intensity.green<-intensity3D(X2+range,Y2+range,Z2+range,bw=bw.int,psz=psz)
  K<-K.cross.3D(X+range,Y+range,Z+range,X2+range,Y2+range,Z2+range,psz=psz,width=width,intensity=intensity.red,intensity2=intensity.green,parallel=parallel)
  
  #for (i in 1:N)
  #envelope[,i]<-K$y-envelope[,i]
  
  e.max<-apply(envelope,1,max)
  e.min<-apply(envelope,1,min)
  return(list("y"=K$y,"x"=K$x,"max"=e.max,"min"=e.min,"envelope"=envelope))
}


#opening your data file containing x, y, z coordinates and one attribute (for example two types of objects: object1 and object2)
data <- read.csv(file = 'yourfile/file.csv', 
                 colClasses = c(x = "numeric", y = "numeric", z = "numeric", attribute = "character", ...), sep = ";")


#data preparation before the use of the K.cross.3D.envelope function
object1 <- dplyr::select(data, x,y,z,attribute)
object1 <- dplyr::filter(object1, attribute %in% c('object1')) #replace 'object1' with the name of the attribute you are testing

object2 <- dplyr::select(data, x,y,z,attribute)
object2 <- dplyr::filter(object2, attribute == 'object2') #replace 'object2' with the name of the attribute you are testing

#after you selected your two attributes (object1 and object2) above, you don't need to edit anything in the following sections

#use of the K.cross.3D.envelope to test the spatial association between object1 and object2
res_kcross3d <- K.cross.3D.envelope(abs(object1$x), abs(object1$y), abs(object1$z), abs(object2$x), abs(object2$y), abs(object2$z), range = 0.1, bw.int = 1, psz = 5, width = 0.35, use.intensity = FALSE, parallel = FALSE, N = 1000)
res_kcross3d <- data.frame(res_kcross3d)
res_kcross3d <- select(res_kcross3d, x, y)
res_kcross3d <- data.frame(res_kcross3d.x = res_kcross3d$x, res_kcross3d.y = res_kcross3d$y)


#then we want to compare it to a Poisson point pattern (spatial randomness)
#we have to create the Poisson point pattern in the same 3D boxes with the same numbers of points as the objects
#we keep the same 3D boxes for object1 and object2

object1_box <- pp3(object1$x,object1$y,object1$z, box3(c(0,1)))
object1_box <- as.box3(object1_box)

object2_box <- pp3(object2$x,object2$y,object2$z, box3(c(0,1)))
object2_box <- as.box3(object2_box)


#simulation of Poisson point patterns in the same boxes
poisson_object1 <- rpoispp3(nrow(object1), domain = object1_box)
poisson_object1 <- data.frame(x = poisson_object1$data$x, y = poisson_object1$data$y, z = poisson_object1$data$z)

poisson_object2 <- rpoispp3(23, domain = y)
poisson_object2 <- data.frame(x = poisson_object2$data$x, y = poisson_object2$data$y, z = poisson_object2$data$z)


#we test the spatial association of the two latter Poisson point pattern
object1v2_poisson <- K.cross.3D.envelope(abs(poisson_object1$x), abs(poisson_object1$y), abs(poisson_object1$z), abs(poisson_object2$x), abs(poisson_object2$y), abs(poisson_object2$z), range = 0.1, bw.int = 1, psz = 5, width = 0.35, use.intensity = FALSE, parallel = FALSE, N = 1000)
object1v2_poisson <- data.frame(object1v2_poisson)


#plot all the results with ggplot2
#first combining the results tables for both object1/object2 spatial associations and Poisson point patterns association tests
object1v2_poisson2 <-select(object1v2_poisson, x, y, min, max)
object1v2_poisson2 <- mutate(object1v2_poisson2, level = rep(1, times = 500))
object1v2_poisson_res_kcross3d <- cbind(object1v2_poisson2, res_kcross3d)


#second making the plot
kcross_plot_object1v2 <- ggplot(object1v2_poisson_res_kcross3d, aes(x = x, ymin = min, ymax = max, fill = "ligthgrey")) +
  geom_ribbon(alpha = 0.4) +
  scale_fill_manual(values = c("1" = "grey")) +
  theme_light() +
  geom_point(aes(res_kcross3d.x,res_kcross3d.y), size = 0.1) +
  geom_point(aes(x,y), col="red", size = 0.1) +
  theme(legend.position = "none") +
  xlab("r") + ylab("K(r)") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())


#third printing it
#the point pattern observed (spatial association of object1 and object2) corresponds to the black curve
#the grey shade area corresponds to the Poisson point pattern (no association)
#if the black curve is above the grey shade, object 1 and 2 are attracted to one other
#if the black curve is below the grey shade, object 1 and 2 are repulsive
ggdraw() +
  draw_plot(kcross_plot_object1v2, 0, 0, 1, 1)