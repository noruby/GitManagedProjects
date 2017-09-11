speed_of_light <- 299792458 #m/s

#calculation parameter
division_number_frequency <- 4000
division_number_time <- 4000

#property of gold
scattering_rate <- 384 *speed_of_light *100 #1/s (<380 cm^-1 for grain size of >56nm)
plasma_frequency <- 67900 *speed_of_light *100 #1/s (67900 cm^-1)
dielectric_background <- 7.0

#parameter physics
jones_vector_initial <- c(1,1) #ps-basis (normalized)
jones_vector_final <- c(1,1) #ps-basis (normalized)
incident_angle <- 5/180*pi#radian

#Ti:saphire laser
central_wavelength <- 800 *10^(-9) #m
central_frequency <- speed_of_light/central_wavelength #1/s
pulse_duration <- 50 *10^(-15) #s 
repetition_rate <- 1*10^3 #1/s
#operating_frequency <- seq(central_frequency-5/pulse_duration, central_frequency+5/pulse_duration, by=repetition_rate) #1/s
operating_frequency <- seq(central_frequency-10/pulse_duration, central_frequency+10/pulse_duration, length=division_number_frequency) #1/s

Ep_spectrum_incident <- exp(-(operating_frequency-central_frequency)^2*(pulse_duration^2)/16)*jones_vector_initial[1]
Es_spectrum_incident <- exp(-(operating_frequency-central_frequency)^2*(pulse_duration^2)/16)*jones_vector_initial[2]

#Fourier transformation from spectrum to time
time_coordinate <- seq(-5*pulse_duration, 5*pulse_duration, length=division_number_time) #s
Ep_time_incident <- numeric(division_number_time) #0
Es_time_incident <- numeric(division_number_time) #0
for (i in 1:division_number_frequency){
	#trapzoid
	Ep_time_incident <- Ep_time_incident + Ep_spectrum_incident[i]*exp(1i*operating_frequency[i]*time_coordinate)
	Es_time_incident <- Es_time_incident + Es_spectrum_incident[i]*exp(1i*operating_frequency[i]*time_coordinate)
	#Simpson
	#Ep_time_incident <- Ep_time_incident + i%2* Ep_spectrum_incident[i]*exp(1i*operating_frequency[i]*time_coordinate)
	#Es_time_incident <- Es_time_incident + i%2* Es_spectrum_incident[i]*exp(1i*operating_frequency[i]*time_coordinate)
}

#calculating reflection coefficient
epsilon <- dielectric_background - plasma_frequency^2 / ( operating_frequency^2 + 1i*operating_frequency*scattering_rate )
refractive_index <- sqrt(epsilon)
refraction_angle <- asin(incident_angle/refractive_index)
reflection_s <- (cos(incident_angle)- refractive_index*cos(refraction_angle))/(cos(incident_angle)+ refractive_index*cos(refraction_angle))
reflection_p <- (refractive_index*cos(incident_angle)- cos(refraction_angle))/(refractive_index*cos(incident_angle)+ cos(refraction_angle))
print(epsilon)

 #reflection spectrum
Es_spectrum_refelection <- reflection_s *Es_spectrum_incident
Ep_spectrum_refelection <- reflection_p *Ep_spectrum_incident

#Fourier transformation from spectrum to time
Ep_time_refelection <- numeric(division_number_time) #0
Es_time_refelection <- numeric(division_number_time) #0
for (i in 1:division_number_frequency){
	#trapzoid
	Ep_time_refelection <- Ep_time_refelection + Ep_spectrum_refelection[i]*exp(1i*operating_frequency[i]*time_coordinate)
	Es_time_refelection <- Es_time_refelection + Es_spectrum_refelection[i]*exp(1i*operating_frequency[i]*time_coordinate)
	#Simpson
	#Ep_time_refelection <- Ep_time_refelection + i%2* Ep_spectrum_refelection[i]*exp(1i*operating_frequency[i]*time_coordinate)
	#Es_time_refelection <- Es_time_refelection + i%2* Es_spectrum_refelection[i]*exp(1i*operating_frequency[i]*time_coordinate)
}
print(Ep_spectrum_refelection[2000]/Es_spectrum_refelection[2000])

#analyzing
 detection_amplitude <- (jones_vector_final[1]*Ep_time_refelection+jones_vector_final[2]*Es_time_refelection)
 transition_amplitude <- 
 weak_value <- detection_amplitude/

png("longitudinal_shift_time_domain.png", height=1300, width=1700, res=216)
par(mar = c(3, 3, 1, 1)) #  余白の広さを行数で指定．下，左，上，右の順．
par(mgp = c(1.5, 0.3, 0)) #  余白の使い方．説明，ラベル，軸の位置を行で指定．
plot(time_coordinate*10^15, abs(Ep_time_incident)/max(abs(Ep_time_incident)),type="l",lwd=2, tcl=0.5, col="red",xlab="fs", ylab="amplitude")
par(new=T)
plot(time_coordinate*10^15, abs(detection_amplitude)/max(abs(detection_amplitude)),type="l",lwd=2, tcl=0.5, col="blue",xlab="", ylab="")
#axis(side=3, xlim=c(min(scattering_rate),max(scattering_rate)), at=c(221,261,313,384), labels=c("bulk",170,80,60), tck=1, col="blue", lty="dashed")
#mtext("Grain diameter (nm)*", side=3, line=1.5)
#axis(side=4, xlim=c(min(propagation_length),max(propagation_length)), at=c(9.0,12.0,14.7), labels=c("","",""), tck=1, col="green", lty="dashed")
#text(390,9.0, "wavelength: 800 nm")
dev.off()