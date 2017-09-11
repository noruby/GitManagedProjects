n <- 50
a <- 3
x <- seq(-n, n, length=10000)
y <- (cos(2*pi*x/n)+a*1i*sin(2*pi*x/n))^n

png("superoscillations_sample.png", height=1300, width=1700, res=216)
par(mar = c(3, 3, 1, 1)) #  余白の広さを行数で指定．下，左，上，右の順．
par(mgp = c(1.5, 0.3, 0)) #  余白の使い方．説明，ラベル，軸の位置を行で指定．
plot(x, Arg(y),log="",type="p",lwd=2, tcl=0.5)
dev.off()