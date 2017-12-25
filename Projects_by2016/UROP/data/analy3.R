x <- c()
lambda <- c()
coll <- c("black","red","blue","green","yellow")

for( j in 2:61 ) {

  if ( j%%2 == 0 ) {
    k <- j/2 
    pol <- "Ex"
    next
  } else {
    k <- (j-1)/2
    pol <- "Ez"
    #next
  }

  for( i in 0:16 ) {
    string <- paste("./SPP_excit_coeff_work/raw/SPP_excit_coeff_" , as.character(i) , ".tmn", sep="")
    #string <- "./SPP_excit_coeff_work/raw/SPP_excit_coeff_0.tmn"
    data <- read.table(string,sep=" ",skip=3)
    #print(string)
    lambda[i] = 8.0+0.25*i
    #print(as.character(i))
    x[i]=data[390,j] 
  }

  print(j)
  #print(ramda)
  #print(x)

  l <- (k-1)%%5
  m <- as.integer((k-1)/5+1)
  xcord <- m*100
  zcord <- 4*l +2

  if ( l == 0 ) { 
  string2 <- paste("./SPP_excit_coeff_work/plot2/plot" ,pol,"_",xcord, ".png", sep="") 
  png(string2)
  plot(lambda,x,type="l",ann=F,col=coll[l],ylim=c(0,3))
  par(new=T)
  } else if ( l == 4 ) {
  plot(lambda,x,type="l",ylab=pol,col=coll[l],ylim=c(0,3))
  labels <- paste("Δ=", as.character(seq(2,18,by=4)),"μm", sep="")
  legend( "topleft", legend = labels , col = coll, lty = 1 )
  dev.off()
  } else {
  plot(lambda,x,type="l",ann=F,col=coll[l],ylim=c(0,3))
  par(new=T)
  }
 
}
