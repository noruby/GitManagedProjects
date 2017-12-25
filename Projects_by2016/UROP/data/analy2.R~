  data <- read.table("coeff_base.tmn",sep=" ",skip=3)

  t <- data[,1]
  Ex_in <- data[,140]
  Ez_in <- data[,141]
  Ex_out <- data[,130]
  Ez_out <- data[,131]

  png("plotEx_base.png")
  plot(t,Ex_in,ann=F, type="l",ylim=c(0,4),col="blue")
  par(new=T)
  plot(t,Ex_out,ylab="Ex", type="l",ylim=c(0,4),col="red")
  dev.off()
