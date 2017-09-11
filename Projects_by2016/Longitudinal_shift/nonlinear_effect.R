speed_of_light <- 299792458 #m/s
#calculation parameter
division_number_frequency <- 4000

#Drude parameter of gold
scattering_rate_wavenum <- 384 #cm^-1 (<380 cm^-1 for grain size of >56nm)
scattering_rate <- scattering_rate_wavenum * speed_of_light *100 #1/s
plasma_frequency_wavenum <-67900 #cm^-1
plasma_frequency <- plasma_frequency_wavenum  * speed_of_light *100 #1/s
dielectric_background <- 7.0

#Ti:saphire laser
central_wavelength <- 800 *10^(-9) #m
central_frequency <- speed_of_light/central_wavelength #1/s
pulse_duration <- 50 *10^(-15) #s 
frequency_deviation <- 1/(pulse_duration*2/2) #1/s
repetition_rate <- 1*10^3 #1/s
#operating_frequency <- seq(central_frequency-5/pulse_duration, central_frequency+5/pulse_duration, by=repetition_rate) #1/s
operating_frequency <- seq(central_frequency/100, central_frequency+10/pulse_duration, length=division_number_frequency) #1/s
operating_frequency_wavenum <- operating_frequency /speed_of_light /100 #cm^-1

dielectric_constant <- dielectric_background - plasma_frequency^2 / ( operating_frequency^2 + 1i*operating_frequency*scattering_rate )
k_free_light_wavenum <- operating_frequency_wavenum #cm^-1
k_SPP_wavenum <- operating_frequency_wavenum* sqrt(dielectric_constant/(dielectric_constant+1)) #cm^-1
k_SPP_Re_wavenum <- Re(k_SPP_wavenum)
intensity <- exp(-(operating_frequency-central_frequency)^2/(frequency_deviation^2)/2*2)

png("nonlinear.png", height=1300, width=1700, res=216)
par(mar = c(3, 3, 1, 1)) #  余白の広さを行数で指定．下，左，上，右の順．
par(mgp = c(1.5, 0.3, 0)) #  余白の使い方．説明，ラベル，軸の位置を行で指定．
plot(k_free_light_wavenum, operating_frequency_wavenum, xlim=c(0,max(operating_frequency_wavenum)), ylim=c(0,max(operating_frequency_wavenum)), type="l", lwd=2, tcl=0.5, col="black", xlab="", ylab="")
par(new=T)
plot(k_SPP_Re_wavenum, operating_frequency_wavenum, xlim=c(0,max(operating_frequency_wavenum)), ylim=c(0,max(operating_frequency_wavenum)),type="l", lwd=2, tcl=0.5, col="red", xlab="", ylab="")
par(new=T)
plot(intensity*max(operating_frequency_wavenum)*0.5, k_free_light_wavenum, xlim=c(0,max(operating_frequency_wavenum)), ylim=c(0,max(operating_frequency_wavenum)),type="l", lty="dashed", lwd=2, tcl=0.5, col="green",xlab="wavenumber (cm-1)", ylab="frequency (cm-1)")
par(new=T)
axis(side=2, xlim=c(0,max(operating_frequency_wavenum)), ylim=c(0,max(operating_frequency_wavenum)), at=c(500,1500), labels=c("", "FP"), tck=1, col="blue", lty="dashed")
dev.off()