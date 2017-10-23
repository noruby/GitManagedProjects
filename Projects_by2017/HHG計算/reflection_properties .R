n1 <-1 
n2 <-1.5

incident_angle <- seq(0, pi/2, length=100)
refraction_angle <- asin(n1/n2*sin(incident_angle)) #snell's law 

refraction_angle

reflection_coefficient_Ppolarized <- tan(incident_angle-refraction_angle)/tan(incident_angle+refraction_angle)
reflection_coefficient_Spolarized <- -sin(incident_angle-refraction_angle)/sin(incident_angle+refraction_angle)

xrange <- c(0,90)
yrange <- c(0,1)
png("reflection.png", width = 700, height = 500) 
#plot(data, xlim=xrange, ylim=yrange, col="black")  
#par(new=T)
plot(incident_angle/pi*180, reflection_coefficient_Ppolarized^2, xlim=xrange, ylim=yrange, col="red", type="l")  
par(new=T)
plot(incident_angle/pi*180, reflection_coefficient_Spolarized^2, xlim=xrange, ylim=yrange, col="blue", type="l")  
dev.off() 