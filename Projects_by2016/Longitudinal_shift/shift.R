division_number_angle <- 2000
division_number_x <- 2000

scattering_rate <- 384 #cm^-1 (<380 cm^-1 for grain size of >56nm)
plasma_frequency <-67900 #cm^-1
dielectric_background <- 7.0

jones_vector <- c(1,0) #xy-basis
wavelength <- 800 *10^(-9) #nm
operating_frequency <- 1/(wavelength*100) #cm^-1

incident_angle_mean <- 85/180*pi#radian
incident_angle_divergence <- 5/180*pi#radian
incident_angle <- seq(incident_angle_mean-5*incident_angle_divergence, incident_angle_mean+5*incident_angle_divergence, length=division_number_angle)
kx_angle <- 2*pi/wavelength*sin(incident_angle)
ky_angle <- 2*pi/wavelength*(-cos(incident_angle))
Ex_angle <- exp(-(incident_angle-incident_angle_mean)^2/(4*incident_angle_divergence^2))*cos(incident_angle)
Ey_angle <- exp(-(incident_angle-incident_angle_mean)^2/(4*incident_angle_divergence^2))*sin(incident_angle)
Ez_angle <- exp(-(incident_angle-incident_angle_mean)^2/(4*incident_angle_divergence^2))

spot_size <- wavelength/pi/incident_angle_divergence
x_coordinate <- seq(-2*spot_size/cos(incident_angle_mean), 2*spot_size/cos(incident_angle_mean), length=division_number_x)

Ex_x_incident_p <- numeric(division_number_x) #0
Ey_x_incident_p <- numeric(division_number_x) #0
Ez_x_incident_s <-  
for (i in 1:division_number_angle){
	Ex_x_incident <- Ex_x_incident + Ex_angle[i]*exp(1i*x_coordinate*kx_angle[i])
	Ey_x_incident <- Ey_x_incident + Ey_angle[i]*exp(1i*x_coordinate*kx_angle[i])
}
E_squared_incident <- abs(Ex_x_incident)^2 + abs(Ey_x_incident)^2

epsilon <- dielectric_background - plasma_frequency^2 / ( operating_frequency^2 + 1i*operating_frequency*scattering_rate )
refractive_index <- sqrt(epsilon)
refraction_angle <- asin(incident_angle/refractive_index)
reflection_s <- (cos(incident_angle)- refractive_index*cos(refraction_angle))/(cos(incident_angle)+ refractive_index*cos(refraction_angle))
reflection_p <- (refractive_index*cos(incident_angle)- cos(refraction_angle))/(refractive_index*cos(incident_angle)+ cos(refraction_angle))
print(epsilon)
#print(refraction_angle)
print(reflection_s)
#print(reflection_p)

Ex_angle_reflection_s <- reflection_s*Ex_angle
Ey_angle_reflection_s <- reflection_s*Ey_angle
Ex_angle_reflection_p <- reflection_p*Ex_angle
Ey_angle_reflection_p <- reflection_p*Ey_angle

Ex_x_reflection_s <- numeric(division_number_x) #0
Ey_x_reflection_s <- numeric(division_number_x) #0
Ex_x_reflection_p <- numeric(division_number_x) #0
Ey_x_reflection_p <- numeric(division_number_x) #0
for (i in 1:division_number_angle){
	Ex_x_reflection_s <- Ex_x_reflection_s + Ex_angle_reflection_s[i]*exp(1i*x_coordinate*kx_angle[i])
	Ey_x_reflection_s <- Ey_x_reflection_s + Ey_angle_reflection_s[i]*exp(1i*x_coordinate*kx_angle[i])
	Ex_x_reflection_p <- Ex_x_reflection_p + Ex_angle_reflection_p[i]*exp(1i*x_coordinate*kx_angle[i])
	Ey_x_reflection_p <- Ey_x_reflection_p + Ey_angle_reflection_p[i]*exp(1i*x_coordinate*kx_angle[i])
}
E_squared_reflection_s <- abs(Ey_x_reflection_s)^2 + abs(Ey_x_reflection_s)^2
E_squared_reflection_p <- abs(Ey_x_reflection_p)^2 + abs(Ey_x_reflection_p)^2

png("longitudinal_shift.png", height=1300, width=1300, res=216)
par(mar = c(3, 3, 3, 1)) #  余白の広さを行数で指定．下，左，上，右の順．
par(mgp = c(1.5, 0.3, 0)) #  余白の使い方．説明，ラベル，軸の位置を行で指定．
plot(x_coordinate*10^6, E_squared_reflection_p/max(E_squared_reflection_p),type="l",lwd=2, tcl=0.5, col="red",xlab=expression("x"~(mu~m)), ylab="amplitude")
par(new=T)
#axis(side=3, xlim=c(min(scattering_rate),max(scattering_rate)), at=c(221,261,313,384), labels=c("bulk",170,80,60), tck=1, col="blue", lty="dashed")
#mtext("Grain diameter (nm)*", side=3, line=1.5)
#axis(side=4, xlim=c(min(propagation_length),max(propagation_length)), at=c(9.0,12.0,14.7), labels=c("","",""), tck=1, col="green", lty="dashed")
#text(390,9.0, "wavelength: 800 nm")
dev.off()

print(sum(E_squared_incident*x_coordinate))
print(sum(E_squared_reflection_s*x_coordinate))
print(sum(E_squared_reflection_p*x_coordinate))

