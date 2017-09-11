data=read.table("./RCWA_raw.dat", header=FALSE)

y <- seq(0.03, 1.4, length.out = nrow(data))

x <- seq(16.9, 18.0, length.out = ncol(data))

z1 <- matrix(as.numeric(unlist(data)),nrow(data))

z1 <- t(z1)

#z2 <- z1



#par <- 100

#div <- 0.6

#delta_y <- y[2]-y[1]



#for(i in 1:nrow(data)){

#  for(j in 1:ncol(data)){

#    sum <- 0.0

#    for(k in 1:par){

#      sum <- sum + dnorm(y[j]-(y[1]-delta_y*k),0,div)*z1[i,1]

#    }

#    for(k in 1:ncol(data)){

#      sum <- sum + dnorm(y[j]-y[k],0,div)*z1[i,k]

#    }

#    for(k in 1:par){

#      sum <- sum + dnorm(y[j]-(y[ncol(data)]+delta_y*k),0,div)*z1[i,ncol(data)]

#    }

#    z2[i,j] <- sum*delta_y

#  }

#}





png("RCWA_Rplot.png", height=1000, width=1441, res=216)

#par(mar = c(3.8, 3.8, 1, 0)) #  余白の広さを行数で指定．下，左，上，右の順．
par(mar = c(2, 2, 1, 1)) #  余白の広さを行数で指定．下，左，上，右の順．

#par(mgp = c(2.2, 0.7, 0)) #  余白の使い方．説明，ラベル，軸の位置を行で指定．
par(mgp = c(2.2, 0.7, 0)) #  余白の使い方．説明，ラベル，軸の位置を行で指定．

#filled.contour(x,y,z1, xlab="Incident angle [degree]", ylab=expression(paste("Grating depth [",mu,"m]")),cex.lab=1.1,cex.axis =0.5)
filled.contour(x,y,z1, xlab="", ylab="")


dev.off()



#png("RCWA_Rplot_convoluted.png")

#filled.contour(x,y,z2,xlim = range(x, finite = TRUE),ylim = range(y, finite = TRUE),zlim = range(z2, finite = TRUE))

#dev.off()
