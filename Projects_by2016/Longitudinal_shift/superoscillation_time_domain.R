speed_of_light <- 299792458 #m/s

#calculation parameter
division_number_frequency <- 4000
division_number_time <- 50000

#Ti:saphire laser
central_wavelength <- 800 *10^(-9) #m
central_frequency <- speed_of_light/central_wavelength #1/s
pulse_duration <- 50 *10^(-15) #s 
repetition_rate <- 1*10^3 #1/s
#operating_frequency <- seq(central_frequency-5/pulse_duration, central_frequency+5/pulse_duration, by=repetition_rate) #1/s
operating_frequency <- seq(central_frequency-10/pulse_duration, central_frequency+10/pulse_duration, length=division_number_frequency) #1/s

spectrum_incident <- exp(-(operating_frequency-central_frequency)^2*(pulse_duration^2)/16)

#Fourier transformation from spectrum to time
time_coordinate <- seq(-5*pulse_duration, 5*pulse_duration, length=division_number_time) #s
Et_incident <- numeric(division_number_time) #0
for (i in 1:division_number_frequency){
	#trapzoid
	Et_incident <- Et_incident + spectrum_incident[i]*exp(1i*operating_frequency[i]*time_coordinate)
	#Simpson
	#Et_incident <- Et_incident + i%2* spectrum_incident[i]*exp(1i*operating_frequency[i]*time_coordinate)
}

#arg
Arg_t <- numeric(division_number_time) #0
delta_Arg <- numeric(division_number_time) #0
for(j in 1:(division_number_time-1)){
	delta_Arg[j] <- Arg(Et_incident[j+1]-Et_incident[j])
}
for(j in 1:(division_number_time-2)){
	Arg_t[j] <- delta_Arg[j+1]-delta_Arg[j]
	if( Arg_t[j] < 0 ){Arg_t[j] <- Arg_t[j] +2*pi }
}


png("superoscillation_in_time.png", height=1300, width=1700, res=216)
par(mar = c(3, 3, 1, 1)) #  余白の広さを行数で指定．下，左，上，右の順．
par(mgp = c(1.5, 0.3, 0)) #  余白の使い方．説明，ラベル，軸の位置を行で指定．
plot(time_coordinate*10^15, abs(Et_incident)/max(abs(Et_incident)),log="y",type="l",lwd=2, tcl=0.5, col="red",xlab="fs", ylab="amplitude")
par(new=T)
plot(time_coordinate*10^15, Arg(Et_incident),type="l",lwd=2, tcl=0.5, col="blue",xlab="", ylab="")
#axis(side=3, xlim=c(min(scattering_rate),max(scattering_rate)), at=c(221,261,313,384), labels=c("bulk",170,80,60), tck=1, col="blue", lty="dashed")
#mtext("Grain diameter (nm)*", side=3, line=1.5)
#axis(side=4, xlim=c(min(propagation_length),max(propagation_length)), at=c(9.0,12.0,14.7), labels=c("","",""), tck=1, col="green", lty="dashed")
#text(390,9.0, "wavelength: 800 nm")
dev.off()