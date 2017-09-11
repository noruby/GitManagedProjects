Reflection_efficiency <- 0.9
Conversion_efficiency <- sqrt(1-R*R)
grating_pitch <- 1
theta <- 20/180*pi

time_division <- 1000
pulse_duration <- 150 * 10^(-15)
#150 fs
groove_number <- 100
pulse_width <- 50
intensity <- 1
wavelength <- 1.2 
beam_center <- 50

wavevector <- 2*pi/wavelength*c(cos(theta),sin(theta))
SPP_field <- numeric(groove_number)
Input_field <- numeric(groove_number)
Output_field <- numeric(groove_number)
for(i in 2:time_division){
  for(j in 1:groove_number)
  delta <- wavevector*6/time_division
  length <- cos(theta)*distance - sin(theta)*
  width <- cos(theta)*distance
  Input_field[j] <- intensity *exp(-(distance_from_center-beam_center)^2/pulse_width^2 - delay/pulse_length^2 *exp( %*% wavevector)
  distance <- i *grating_pitch
  phase <- i *wavelength *sin(theta) 
  Output_field[i] <- Input_field[i]*Reflection_efficiency + SPP_field[i-1]*Conversion_efficiency
  SPP_field[i] <- -Input_field[i]*Conversion_efficiency + SPP_field[i-1]*Reflection_efficiency
  }
}

plot[]