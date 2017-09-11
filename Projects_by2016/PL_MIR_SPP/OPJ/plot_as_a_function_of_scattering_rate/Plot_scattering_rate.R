speed_of_light <- 299792458 #m/s
permittivity_vacuum <- 8.854187817*10^(-12) #F/m
num<-50

scattering_rate <- seq(200,400, length=75) #cm^-1 (<380 cm^-1 for grain size of >56nm)
#scattering_time <- 10^{13}/ 2/ pi / speed_of_light/ scattering_rate #fs
plasma_frequency <-70000 #cm^-1
dielectric_background <- 8

operating_frequency <- 943 #cm^-1

epsilon <- dielectric_background - plasma_frequency^2 / ( operating_frequency^2 + 1i*operating_frequency*scattering_rate )

wavenumber_x <- 2*pi* operating_frequency* sqrt(epsilon/(1+epsilon)) #1/cm
propagation_length <- 1/(2*Im( wavenumber_x ))

wavenumber_z_air <- 2*pi* operating_frequency* sqrt(1/(1+epsilon)) #1/cm
penetration_length_air <- -1/(2*Im( wavenumber_z_air ))

wavenumber_z_metal <- 2*pi* operating_frequency* sqrt(epsilon^2/(1+epsilon)) #1/cm
penetration_length_metal <- 1/(2*Im( wavenumber_z_metal ))

#conductivity_0 <- plasma_frequency^2/scattering_rate*(100*speed_of_light)*permittivity_vacuum
conductivity <- (2*pi*100*plasma_frequency*speed_of_light)*permittivity_vacuum *plasma_frequency/(scattering_rate-1i*operating_frequency)

png("real_dielectric_constant.png", height=1000, width=1000, res=216)
par(mar = c(3, 3, 3, 1)) #  余白の広さを行数で指定．下，左，上，右の順．
par(mgp = c(1.5, 0.3, 0)) #  余白の使い方．説明，ラベル，軸の位置を行で指定．
plot(scattering_rate,abs(Re(epsilon)),type="l",lwd=2, tcl=0.5, col="red",xlim=c(min(scattering_rate),max(scattering_rate)),ylim=c(0,max(abs(Re(epsilon)))),xlab=expression("Scattering rate"~(cm^{-1})), ylab=expression("|Re("~epsilon~")|"))
par(new=T)
axis(side=3, xlim=c(min(scattering_rate),max(scattering_rate)), at=c(221,261,313,384), labels=c("bulk",170,80,60), tck=1, col="blue", lty="dashed")
mtext("Grain diameter (nm)*", side=3, line=1.5)
dev.off()

png("imaginary_dielectric_constant.png", height=1000, width=1000, res=216)
par(mar = c(3, 3, 3, 1)) #  余白の広さを行数で指定．下，左，上，右の順．
par(mgp = c(1.5, 0.3, 0)) #  余白の使い方．説明，ラベル，軸の位置を行で指定．
plot(scattering_rate,Im(epsilon),type="l",lwd=2, tcl=0.5, col="red",ylim=c(0,max(Im(epsilon))),xlab=expression("Scattering rate"~(cm^{-1})), ylab=expression("Im("~epsilon~")"))
par(new=T)
axis(side=3, xlim=c(min(scattering_rate),max(scattering_rate)), at=c(221,261,313,384), labels=c("bulk",170,80,60), tck=1, col="blue", lty="dashed")
mtext("Grain diameter (nm)*", side=3, line=1.5)
dev.off()

png("propagation_length.png", height=1300, width=1300, res=216)
par(mar = c(3, 3, 3, 1)) #  余白の広さを行数で指定．下，左，上，右の順．
par(mgp = c(1.5, 0.3, 0)) #  余白の使い方．説明，ラベル，軸の位置を行で指定．
plot(scattering_rate,10*propagation_length,type="l",lwd=2, tcl=0.5, col="red",ylim=c(0,10*max(propagation_length)),xlab=expression("Scattering rate"~(cm^{-1})), ylab="Propagation length (mm)")
par(new=T)
axis(side=3, xlim=c(min(scattering_rate),max(scattering_rate)), at=c(221,261,313,384), labels=c("bulk",170,80,60), tck=1, col="blue", lty="dashed")
mtext("Grain diameter (nm)*", side=3, line=1.5)
axis(side=4, xlim=c(min(propagation_length),max(propagation_length)), at=c(9.0,12.0,14.7), labels=c("","",""), tck=1, col="green", lty="dashed")
text(390,9.0, "As grown")
text(380,12.0, "Once annealed")
text(380,14.7, "Twice annealed")
dev.off()

png("penetration_length_air.png", height=1000, width=1000, res=216)
par(mar = c(3, 3, 3, 1)) #  余白の広さを行数で指定．下，左，上，右の順．
par(mgp = c(1.5, 0.3, 0)) #  余白の使い方．説明，ラベル，軸の位置を行で指定．
plot(scattering_rate,10*penetration_length_air,type="l",lwd=2, tcl=0.5, col="red",ylim=c(0,10*max(penetration_length_air)),xlab=expression("Scattering rate"~(cm^{-1})), ylab="Penetration length into air (mm)")
par(new=T)
axis(side=3, xlim=c(min(scattering_rate),max(scattering_rate)), at=c(221,261,313,384), labels=c("bulk",170,80,60), tck=1, col="blue", lty="dashed")
mtext("Grain diameter (nm)*", side=3, line=1.5)
dev.off()

png("penetration_length_metal.png", height=1000, width=1000, res=216)
par(mar = c(3, 3, 3, 1)) #  余白の広さを行数で指定．下，左，上，右の順．
par(mgp = c(1.5, 0.3, 0)) #  余白の使い方．説明，ラベル，軸の位置を行で指定．
plot(scattering_rate,10^7*penetration_length_metal,type="l",lwd=2, tcl=0.5, col="red",ylim=c(0,10^7*max(penetration_length_metal)),xlab=expression("Scattering rate"~(cm^{-1})), ylab="Penetration length into metal (nm)")
par(new=T)
axis(side=3, xlim=c(min(scattering_rate),max(scattering_rate)), at=c(221,261,313,384), labels=c("bulk",170,80,60), tck=1, col="blue", lty="dashed")
mtext("Grain diameter (nm)*", side=3, line=1.5)
dev.off()

png("conductivity.png", height=1000, width=1000, res=216)
par(mar = c(3, 3, 3, 1)) #  余白の広さを行数で指定．下，左，上，右の順．
par(mgp = c(1.5, 0.3, 0)) #  余白の使い方．説明，ラベル，軸の位置を行で指定．
plot(scattering_rate,Re(conductivity)/10^6,type="l",lwd=2, tcl=0.5, col="red",ylim=c(0,max(Re(conductivity))/10^6),xlab=expression("Scattering rate"~(cm^{-1})), ylab=expression("Electric conductivity"~(10^6~Omega^{-1}~m^{-1})))
par(new=T)
axis(side=3, xlim=c(min(scattering_rate),max(scattering_rate)), at=c(221,261,313,384), labels=c("bulk",170,80,60), tck=1, col="blue", lty="dashed")
mtext("Grain diameter (nm)*", side=3, line=1.5)
dev.off()

