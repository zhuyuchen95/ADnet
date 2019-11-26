library(KernSmooth)
library(akima)
denPre2D <- function(edge,DataFit,method=c("integers","bandwidth")[2]){
  index<-sample(1:nrow(DataFit),nrow(DataFit)/2)
  x <- DataFit[,edge[1]]
  y <- DataFit[,edge[2]]
  N <- max(c(x,y))

  
  Width <- c(bw.nrd0(x),bw.nrd0(y))

  
  gridSize <- switch(method,
                     integers  = c(N, N),
                     bandwidth = ceiling(N / c(min(Width[1]),min(Width[2]))))
  gridSize <- pmax(gridSize,10) # make sure there are at least 100 points in total
  

  
  #Compute a 2D Binned Kernel Density Estimate
  BSmooth <- bkde2D(x=cbind(x, y), bandwidth=Width, gridsize=gridSize)
  P <- BSmooth$fhat

  
  #Compute a Binned Kernel Density Estimate
  USmoothx<-bkde(x)
  USmoothy<-bkde(y)
  
  P <- pmax(P, 1e-10)
  
  USmoothx$y<-pmax(USmoothx$y, 1e-10)
  USmoothy$y<-pmax(USmoothy$y, 1e-10)

  
  fit <- bicubic(x=BSmooth$x1, y=BSmooth$x2, z=BSmooth$fhat, x0=x,y0=y)

  Ufitx <-aspline(USmoothx$x,USmoothx$y,x)
  Ufity <-aspline(USmoothy$x,USmoothy$y,y)

  
  Ufitx$y <- pmax(Ufitx$y, 1e-10)
  Ufity$y <- pmax(Ufity$y, 1e-10)
  fit$z <- pmax(fit$z, 1e-10)

  denpre<- log(fit$z/(Ufitx$y*Ufity$y))
  denpre
}


# make sure there are no zeros in the smooth function (since we will take a log of that)

