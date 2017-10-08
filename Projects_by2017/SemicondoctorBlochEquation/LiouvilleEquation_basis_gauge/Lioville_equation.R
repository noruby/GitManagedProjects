#I use SI unit system for calculation
# first written in 6 Sep. 2017
# modified in 11 Sep. 2017
# modified in 7 Oct. 2017
# modified for basis-gauge independent calculation

options(scipen=0)

#Physical constants 
Plank_constant_bar <- 6.62607004 * 10^-34 #m^2 kg / s
Boltzmann_constant <-  1.38064852 *10^-23 # m^2 kg s^-2 K^-1
speed_of_light  <- 299792458 # m/s
elementary_charge <- 1.60217662 * 10^-19 # Coulombs
mass_of_electron <-  9.10938356 * 10^-31 #kg

#Parameters of experiment / material properties
Temperature <- 300 #K
Fermi_energy <- # Chemical potential
lattice_constant_a_of_GaSe <- 3.755 * 2/sqrt(3) *10^-10# m (3.755 ångström)

# Time t
range_t <- 100 *10^(-15) # fs
div_num_t <- 10000 # n
delta_t <-  range_t/ div_num_t
t_vector <- seq( -range_t/2, range_t/2, length= (div_num_t+1))
t_vector2 <- seq( -range_t/2, range_t/2, length= (2*div_num_t+1))
t_vector4 <- seq( -range_t/2, range_t/2, length= (4*div_num_t+1))
frequency <- seq(0, 1/delta_t, length= (div_num_t+1))

# External electric field (assuming the gaussian envelope)
E_amplitude_max <- 10* 10^8 # V/m
envelope_width <- 30* 10^(-15) # 1/e^2 width; 50 fs 
carrier_frequency <- 100 *10^12 # 33 THz
carrier_envelope_phase <- 0 # [rad]
Electric_field <- function(t){
	EF <- E_amplitude_max*exp(-t^2/envelope_width^2) *cos(2*pi* carrier_frequency*t +carrier_envelope_phase) 
}
E_vector <-  Electric_field(t_vector)
E_vector2 <-  Electric_field(t_vector2)
E_vector4 <-  Electric_field(t_vector4)

# External electric field (step function)
#E_vector[1:(div_num_t+1)] <-  rep(0,div_num_t+1)
#E_vector[101:(div_num_t+1)] <-  rep(E_amplitude_max, div_num_t-99)
#E_vector2[1:(2*div_num_t+1)] <-  rep(0,2*div_num_t+1)
#E_vector2[201:(2*div_num_t+1)] <-  rep(E_amplitude_max, 2*div_num_t-199)

#Vector potential
A_vector <- rep(0, div_num_t+1)
for(i in 1:div_num_t){
	A_vector[i+1] <- A_vector[i] - (E_vector2[2*i+1]+4*E_vector2[2*i]+E_vector2[2*1-1])*delta_t/6
}
A_vector2 <- rep(0, 2*div_num_t+1)
for(i in 1:(2*div_num_t)){
	A_vector2[i+1] <- A_vector2[i] - (E_vector4[2*i+1]+4*E_vector4[2*i]+E_vector4[2*1-1])*delta_t/12
}
#debug plot
if(0){
print("Exporting debug_ElectricF_VectorP.png ...")
png("./debug_ElectricF_VectorP.png", width = 700, height = 500)  
plot(t_vector*10^15, E_vector, xlab="t (fs)", ylab= "J", col="red", type="l")  
par(new=T)
plot(t_vector*10^15, A_vector, xlab="", ylab="", col="blue", type="l")  
par(new=T)
plot(t_vector2*10^15, A_vector2, xlab="", ylab="", col="green", type="l")  
dev.off() 
}	

# fixed k normalized by inverse lattice constant 1/a in length gauge
div_num_ka <- 300 # m
delta_k <- 2*pi / lattice_constant_a_of_GaSe / div_num_ka 
ka_vector_redundant <- seq( -pi, pi, length= (div_num_ka+1)) #BZの両端で冗長
ka_vector <- ka_vector_redundant[1:div_num_ka]

#band structure of GaSe (two band model assuming the Tight binding model (1D))
band_width_valence <- 1 *1.60218*10^-19 # J (1eV)
band_width_conduction <- 1.5 *1.60218*10^-19 # J (1.5eV) 
band_gap <- 2 *1.60218*10^-19 # J (2eV) 
energy_valence <- function(ka){
	band_width_valence* (cos( ka ) -1) /2 # J
}
energy_valence_vector <- energy_valence(ka_vector)
energy_conduction <- function(ka){
	band_width_conduction* (1 - cos( ka ) )/2  +  band_gap # J 
}
energy_conduction_vector <- energy_conduction(ka_vector)
energy_difference <- function(ka){
	energy_conduction(ka) - energy_valence(ka) # J 
}
energy_difference_vector <- energy_difference(ka_vector)
if(0){
print("Exporting debug_band_structure_velocity.png ...")
png("./debug_band_structure.png", width = 700, height = 500)  
krange <- c(-pi, pi)
eVrange <- c(-2,6)
plot(ka_vector, energy_valence_vector/(1.60218*10^-19), xlim=krange, ylim=eVrange, xlab="ka (rad)", ylab= "energy (eV)", col="red", type="l")  
par(new=T)
plot(ka_vector, energy_conduction_vector/(1.60218*10^-19), xlim=krange, ylim=eVrange, xlab="", ylab="", col="blue", type="l")  
par(new=T)
plot(ka_vector, energy_difference_vector/(1.60218*10^-19), xlim=krange, ylim=eVrange, xlab="", ylab="", col="green", type="l")  
dev.off() 
}

#microscopic current
w_1 <- rep(0, div_num_ka)
w1 <- rep(0, div_num_ka)
grad_k_midpoint <- function(w0){
	w_1[1] <- w0[div_num_ka]
	w_1[2:div_num_ka] <- w0[1:(div_num_ka-1)]
	w1[div_num_ka] <- w0[1]
	w1[1:(div_num_ka-1)] <- w0[2:div_num_ka]
	return ((w1-w_1)/ 2 / delta_k)
}
grad_k <- function(w0){
	w1[div_num_ka] <- w0[1]
	w1[1:(div_num_ka-1)] <- w0[2:div_num_ka]
	return ((w0-w1) / delta_k)
}
grad_enegy_valence_vector <- grad_k_midpoint(energy_valence_vector)
grad_enegy_conduction_vector <- grad_k_midpoint(energy_conduction_vector)
if(0){
print("Exporting debug_grad_energy.png ...")
png("./debug_grad_energy.png", width = 700, height = 500)  
krange <- c(-pi, pi)
eVrange <- 0.03*eVrange*(1.60218*10^-19)/delta_k
plot(ka_vector, grad_enegy_valence_vector, xlim=krange, ylim=eVrange, xlab="BZ", ylab= "energy (eV)", col="red", type="l")  
par(new=T)
plot(ka_vector,grad_enegy_conduction_vector, xlim=krange, ylim=eVrange, xlab="", ylab="", col="blue", type="l")  
dev.off() 
}

#dipole transition matrix
Rabi_energy_max <- 0.6*1.60218*10^-19 #1.2eV; converstion from eV to J
transition_matrix <- Rabi_energy_max/(50*10^8)/elementary_charge
#transition_matrix <- Rabi_energy_max/E_amplitude_max/elementary_charge/10
d_vv <- rep(0, div_num_ka)
d_cc <- rep(0, div_num_ka)
d_vc <- rep(transition_matrix, div_num_ka)
#d_vc <- rep(0, div_num_ka)
if(0){
print("Exporting debug_dipole_transition.png ...")
png("./debug_dipole_transition.png", width = 700, height = 500)  
plot(ka_vector, Re(d_vv) , xlim=krange, xlab="ka (rad)", ylab= "dipole_transition", col="red", type="l")  
par(new=T)
plot(ka_vector,Re(d_cc), xlim=krange, xlab="", ylab="", col="blue", type="l")
par(new=T)
plot(ka_vector,Re(d_vc), xlim=krange, xlab="", ylab="", col="green", type="l")
dev.off()   
}

#parameters of calculation
Courant_number <- delta_t * E_amplitude_max *  elementary_charge / Plank_constant_bar / delta_k 
print("Courant number (length gauge): " )
print( Courant_number )

#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------
#calculate ---------------------------------------------------------------------------------------------------------------
gauge <- "length"
#gauge <- "velocity"
#relaxation
relaxation_constant_length <-  20 *10^(-15) # 7fs 
#dephasing
dephasing_constant <- 50000 *10^(-15) # 1.1fs 
#relxation time
#relaxation_constant_velocity <-7*10^-15 #fs


if(gauge=="length"){


diff_number_valence_length <- function(E, mpl, nvl){
	dnvl <- - E *elementary_charge*( grad_k_midpoint(nvl) + 2*Im( mpl* Conj(d_vc)) ) * delta_t/ Plank_constant_bar
	#dnvl <- dnvl - (nvl-1)* relaxation_constant_length * delta_t
	return (dnvl)
}
diff_number_conduction_length <- function(E, mpl, ncl){
	dncl <- - E *elementary_charge*( grad_k_midpoint(nvl) - 2*Im( mpl* Conj(d_vc)) ) * delta_t/ Plank_constant_bar
	#dncl <- dncl - ncl / relaxation_constant_length * delta_t
	return (dncl)
}
diff_micro_polarization_length <- function(E, mpl, nvl, ncl){
	dmpl <- ( 1i* energy_difference_vector *mpl - elementary_charge* E*( grad_k_midpoint(mpl) + 1i*(d_vv-d_cc)*mpl +1i*d_vc*(ncl- nvl) ) ) * delta_t/ Plank_constant_bar
	#dmpl <- dmpl -mpl / dephasing_constant * delta_t
	return (dmpl)
}

nvl_kt <- array(1, dim=c(div_num_ka, (div_num_t+1))) #number_valence_length
ncl_kt <- array(0, dim=c(div_num_ka, (div_num_t+1))) #number_conduction_length
mpl_kt  <- array(0, dim=c(div_num_ka, (div_num_t+1))) #micro_polarization_length

nvl <- rep(1, div_num_ka)
ncl <- rep(0, div_num_ka) 
mpl <- rep(0, div_num_ka)
#solving Lioville equation with the length gauge and bloch-state representation (SBE) by the 4th Runge-Kutta method
print("starting calculation with length gauge")
for( j in 1:(div_num_t) ){
	E <- E_vector2[2*j-1]
	E_0 <- E_vector2[2*j]
	E_1 <- E_vector2[2*j]	
	E_2 <- E_vector2[2*j+1]	

	dnvl_0 <- diff_number_valence_length(E, mpl, nvl)
	dncl_0 <- diff_number_conduction_length(E, mpl, nvl)
	dmpl_0 <- diff_micro_polarization_length(E, mpl, nvl, ncl)

	dnvl_1 <- diff_number_valence_length(E_0, mpl+dmpl_0/2, nvl+dnvl_0/2)
	dncl_1 <- diff_number_conduction_length(E_0, mpl+dmpl_0/2, ncl+dncl_0/2)
	dmpl_1 <- diff_micro_polarization_length(E_0, mpl+dmpl_0/2, nvl+dnvl_0/2, ncl+dncl_0/2)

	dnvl_2 <- diff_number_valence_length(E_1, mpl+dmpl_1/2, nvl+dnvl_1/2)
	dncl_2 <- diff_number_conduction_length(E_1, mpl+dmpl_1/2, ncl+dncl_1/2)
	dmpl_2 <- diff_micro_polarization_length(E_1, mpl+dmpl_1/2, nvl+dnvl_1/2, ncl+dncl_1/2)

	dnvl_3 <- diff_number_valence_length(E_2, mpl+dmpl_2, nvl+dnvl_2)
	dncl_3 <- diff_number_conduction_length(E_2, mpl+dmpl_2, ncl+dncl_2)
	dmpl_3 <- diff_micro_polarization_length(E_2, mpl+dmpl_2, nvl+dnvl_2, ncl+dncl_2)

	nvl <- nvl + (dnvl_0 +2*dnvl_1 +2*dnvl_2 +dnvl_3)/6
	ncl <- ncl + (dncl_0 +2*dncl_1 +2*dncl_2 +dncl_3)/6
	mpl <- mpl + (dmpl_0 +2*dmpl_1 +2*dmpl_2 +dmpl_3)/6

	nvl_kt[,j+1] <- nvl
	ncl_kt[,j+1] <- ncl
	mpl_kt[,j+1] <- mpl
}

Current_valence_length <- rep(0, (div_num_t+1))
Current_conduction_length <- rep(0, (div_num_t+1))
Polarization_current_length <- rep(0, (div_num_t+1))
for(i in 1:(div_num_ka)){
	Current_valence_length <-  Current_valence_length - elementary_charge/Plank_constant_bar* grad_enegy_valence_vector[i] * nvl_kt[i,] *delta_k
	Current_conduction_length <- Current_conduction_length - elementary_charge/Plank_constant_bar* grad_enegy_conduction_vector[i] * ncl_kt[i,] *delta_k 
	Polarization_current_length <- Polarization_current_length + elementary_charge/Plank_constant_bar* energy_difference_vector[i]* Im(mpl_kt[i,] * d_vc[i]) *delta_k
}

Current_length <- Current_valence_length + Current_conduction_length
E_HHG_length <- Current_length + Polarization_current_length

E_spectrum <- fft(E_vector)
E_HHG_length_spectrum <- fft(E_HHG_length)
Current_length_spectrum <-  fft(Current_length)
Polarization_current_length_spectrum <- fft(Polarization_current_length)

}

#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------

gauge<-"velocity"
if(gauge=="velocity"){

diff_number_valence_velocity <- function(A, mpv){
	dnvv <- -2i*elementary_charge* A *Re(mpv*Conj(d_vc)) *energy_difference_vector* delta_t / Plank_constant_bar^2
	#dnvv <- dnvv - (nvv-1) /relaxation_constant_velocity *delta_t
	return (dnvv)
}
diff_number_conduction_velocity <- function(A, mpv){
	dncv <- 2i*elementary_charge* A *Re(mpv*Conj(d_vc)) *energy_difference_vector* delta_t / Plank_constant_bar^2
	#dncv <- dncv - ncv /relaxation_constant_velocity *delta_t
	return (dncv)
}
diff_micro_polarization_velocity <- function(A, mpv, nvv, ncv){
	dmpv <- ( 1i *energy_difference_vector *mpv /Plank_constant_bar
	-1i* elementary_charge* A /Plank_constant_bar^2 * ( energy_difference_vector* d_vc * (ncv-nvv) 
	+(grad_enegy_conduction_vector - grad_enegy_valence_vector)*mpv ) ) * delta_t
	#dmpv <- dmpv -mpv / dephasing_constant * delta_t #dephasing
	return (dmpv)
}

nvv_kt <- array(1, dim=c(div_num_ka, (div_num_t+1))) #number_valence_velocity
ncv_kt <- array(0, dim=c(div_num_ka, (div_num_t+1))) #number_conduction_velocity
mpv_kt  <- array(0, dim=c(div_num_ka, (div_num_t+1))) #micro_polarization_velocity

nvv <- rep(1, div_num_ka)
ncv <- rep(0, div_num_ka) 
mpv <- rep(0, div_num_ka)
#solving Lioville equation with the velocity gauge and bloch-state representation (SBE) by the 4th Runge-Kutta method

print("starting calculation with velocity gauge")
for( j in 1:(div_num_t) ){
	A <- A_vector2[2*j-1]
	A_0 <- A_vector2[2*j]
	A_1 <- A_vector2[2*j]	
	A_2 <- A_vector2[2*j+1]	
	
	dnvv_0 <- diff_number_valence_velocity(A, mpv)
	dncv_0 <- diff_number_conduction_velocity(A, mpv)
	dmpv_0 <- diff_micro_polarization_velocity(A, mpv, nvv, ncv)

	dnvv_1 <- diff_number_valence_velocity(A_0, mpv+dmpv_0/2)
	dncv_1 <- diff_number_conduction_velocity(A_0, mpv+dmpv_0/2)
	dmpv_1 <- diff_micro_polarization_velocity(A_0, mpv+dmpv_0/2, nvv+dnvv_0/2, ncv+dncv_0/2)
	
	dnvv_2 <- diff_number_valence_velocity(A_1, mpv+dmpv_1/2)
	dncv_2 <- diff_number_conduction_velocity(A_1, mpv+dmpv_1/2)
	dmpv_2 <- diff_micro_polarization_velocity(A_1, mpv+dmpv_1/2, nvv+dnvv_1/2, ncv+dncv_1/2)

	dnvv_3 <- diff_number_valence_velocity(A_2, mpv+dmpv_2)
	dncv_3 <- diff_number_conduction_velocity(A_2, mpv+dmpv_2)
	dmpv_3 <- diff_micro_polarization_velocity(A_2, mpv+dmpv_2, nvv+dnvv_2, ncv+dncv_2)

	nvv <- nvv + (dnvv_0 +2*dnvv_1 +2*dnvv_2 +dnvv_3)/6
	ncv <- ncv + (dncv_0 +2*dncv_1 +2*dncv_2 +dncv_3)/6
	mpv <- mpv + (dmpv_0 +2*dmpv_1 +2*dmpv_2 +dmpv_3)/6

	nvv_kt[,j+1] <- nvv
	ncv_kt[,j+1] <- ncv
	mpv_kt[,j+1] <- mpv

}


Charge_valence_velocity <- rep(0, (div_num_t+1))
Charge_conduction_velocity <- rep(0, (div_num_t+1))
Current_valence_velocity <- rep(0, (div_num_t+1))
Current_conduction_velocity <- rep(0, (div_num_t+1))
Polarization_current_velocity <- rep(0, (div_num_t+1))
Current_vectorP_velocity <- rep(0, (div_num_t+1))
for(i in 1:(div_num_ka)){
	Charge_valence_velocity <- Charge_valence_velocity + nvv_kt[i,] *delta_k
	Charge_conduction_velocity <- Charge_conduction_velocity + ncv_kt[i,]	* delta_k
	Current_valence_velocity <-  Current_valence_velocity - elementary_charge/Plank_constant_bar* grad_enegy_valence_vector[i] * nvv_kt[i,] *delta_k
	Current_conduction_velocity <- Current_conduction_velocity - elementary_charge/Plank_constant_bar* grad_enegy_conduction_vector[i] * ncv_kt[i,] *delta_k 
	Polarization_current_velocity <- Polarization_current_velocity + elementary_charge/Plank_constant_bar* energy_difference_vector[i]* Im(mpv_kt[i,] * d_vc[i]) *delta_k
}
Current_vectorP_velocity <- elementary_charge^2 / mass_of_electron* A_vector * (Charge_valence_velocity+ Charge_conduction_velocity)

Current_velocity <- Current_valence_velocity + Current_conduction_velocity
E_HHG_velocity <- Current_velocity +  Polarization_current_velocity #+ Current_vectorP_velocity

E_spectrum <- fft(E_vector)
E_HHG_velocity_spectrum <- fft(E_HHG_velocity)
Current_velocity_spectrum <-  fft(Current_velocity)
Polarization_current_velocity_spectrum <- fft(Polarization_current_velocity)
Current_vectorP_velocity_spectrum <- fft(Current_vectorP_velocity)
}
#---------------------------------------------------------------------------------------------------------------

E_HHGrange <- c( -max(abs(E_HHG_length)), max(abs(E_HHG_length)) ) 
Erange <- c( -max(E_vector), max(E_vector) ) 
trange <- c(-range_t/2, range_t/2)*10^15
png("./Bloch_lengthG/debug_ElectricF_length.png", width = 700, height = 500)  
plot(t_vector*10^15, Re(E_vector), xlab="", ylab= "Electric field (arb. unit)", xlim=trange, axes=FALSE, yaxt="n",col="black", type="l")  
par(new=T)
plot(t_vector*10^15, Re(E_HHG_length), xlab="t (fs)", ylab="", xlim=trange,  ylim=E_HHGrange, col="red", type="l")  
par(new=T)
plot(t_vector*10^15, Re(Current_length), xlab="", ylab="", xlim=trange, ylim=E_HHGrange, axes=FALSE, yaxt="n", col="blue", type="l")  
par(new=T)
plot(t_vector*10^15, Re(Polarization_current_length), xlab="", ylab= "", xlim=trange, ylim=E_HHGrange, axes=FALSE, yaxt="n", col="green", type="l")  
dev.off()    

E_HHGrange <- c( min(abs(E_HHG_length_spectrum)^2), max(abs(E_HHG_length_spectrum)^2) ) 
frange <- c(0, 20)
png("./Bloch_lengthG/debug_spectrum_length.png", width = 700, height = 500) 
plot(frequency/carrier_frequency, abs(E_spectrum)^2, xlim=frange, log="y", xlab="f/carrier frequency", ylab="power (arb. unit)", axes=FALSE, yaxt="n", col="black", type="l", lwd = 1)  
par(new=T)
plot(frequency/carrier_frequency, abs(E_HHG_length_spectrum)^2, xlim=frange, ylim=E_HHGrange, log="y", xlab="", ylab= "", col="red", type="l", lwd = 1)  
par(new=T)
plot(frequency/carrier_frequency, abs(Current_length_spectrum)^2, xlim=frange, ylim=E_HHGrange, log="y", axes=FALSE, yaxt="n",xlab="", ylab= "", col="blue", type="l", lwd = 1)  
par(new=T)
plot(frequency/carrier_frequency, abs(Polarization_current_length_spectrum)^2, xlim=frange, ylim=E_HHGrange, log="y", axes=FALSE, yaxt="n",xlab="", ylab= "", col="blue", type="l", lwd = 1)  
abline(v=0:20, col='black', lty="dotted")
dev.off()    

E_HHGrange <- c( -max(abs(E_HHG_velocity)), max(abs(E_HHG_velocity)) ) 
Erange <- c( -max(E_vector), max(E_vector) ) 
trange <- c(-range_t/2, range_t/2)*10^15
png("./Bloch_velocityG/debug_ElectricF_velocity.png", width = 700, height = 500)  
plot(t_vector*10^15, Re(E_vector), xlab="", ylab= "Electric field (arb. unit)", xlim=trange, ylim=Erange, axes=FALSE, yaxt="n", col="black", type="l")  
par(new=T)
plot(t_vector*10^15, Re(E_HHG_velocity), xlab="t (fs)", ylab="", xlim=trange,  ylim=E_HHGrange, col="red", type="l")  
par(new=T)
plot(t_vector*10^15, Re(Current_velocity), xlab="", ylab="", xlim=trange, ylim=E_HHGrange, axes=FALSE, yaxt="n", col="blue", type="l")
par(new=T)
plot(t_vector*10^15, Re(Polarization_current_velocity), xlab="", ylab= "", xlim=trange, ylim=E_HHGrange, axes=FALSE, yaxt="n", col="green", type="l")  
par(new=T)
plot(t_vector*10^15, Re(Current_vectorP_velocity), xlab="", ylab="", xlim=trange, ylim=E_HHGrange, axes=FALSE, yaxt="n", col="orange", type="l")  
dev.off()    

E_HHGrange <- c( min(abs(E_HHG_velocity_spectrum)^2), max(abs(E_HHG_velocity_spectrum)^2) ) 
frange <- c(0, 20)
png("./Bloch_velocityG/debug_spectrum_velocity.png", width = 700, height = 500) 
plot(frequency/carrier_frequency, abs(E_spectrum)^2, xlim=frange, log="y", xlab="f/carrier frequency", ylab="power (arb. unit)", axes=FALSE, yaxt="n", col="black", type="l")  
par(new=T)
plot(frequency/carrier_frequency, abs(E_HHG_velocity_spectrum)^2, xlim=frange, ylim=E_HHGrange, log="y", xlab="", ylab= "", col="red", type="l")  
par(new=T)
plot(frequency/carrier_frequency, abs(Current_velocity_spectrum)^2, xlim=frange, ylim=E_HHGrange, log="y", axes=FALSE, yaxt="n",xlab="", ylab= "", col="blue", type="l")  
par(new=T)
plot(frequency/carrier_frequency, abs(Polarization_current_velocity_spectrum)^2, xlim=frange, ylim=E_HHGrange, log="y", axes=FALSE, yaxt="n",xlab="", ylab= "", col="green", type="l")  
par(new=T)
plot(frequency/carrier_frequency, abs(Current_vectorP_velocity_spectrum)^2, xlim=frange, ylim=E_HHGrange, log="y", axes=FALSE, yaxt="n",xlab="", ylab= "", col="orange", type="l") 
abline(v=0:20, col='black', lty="dotted")
dev.off()     


if(0){
print("Exporting debug_np_length.gif ...")
numt <- round( 2*10^-15 /delta_t )
nplim=c(-1.2, 1.2)
library(animation)
saveGIF({
  ani.options(loop = TRUE)
  for (i in seq(1,(div_num_t+1), by=numt)) {
    .main = paste(t_vector[i]*10^15,"(fs)", seq="")
    plot(ka_vector, nvl_kt[,i], xlim=c(-pi,pi), ylim=nplim, xlab="BZ", ylab= "numbers", main=.main,col="green", type="l")
    par(new=T)
    plot(ka_vector, ncl_kt[,i], xlim=c(-pi,pi), ylim=nplim, xlab="", ylab= "", main=.main, col="blue", type="l")
    par(new=T)    
    plot(ka_vector, Re(mpl_kt[,i]), xlim=c(-pi,pi), ylim=nplim, xlab="BZ", ylab= "polarization", main=.main,col="red" , type="l")
  }
}, movie.name="debug_np_length.gif", interval=0.2)

print("Exporting debug_np_velocity.gif ...")
numt <- round( 2*10^-15 /delta_t ) #一コマ 2fs
print("numt" )
print(numt)
nplim=c(-1.2, 1.2)
library(animation)
saveGIF({
  ani.options(loop = TRUE)
  for (i in seq(1,(div_num_t+1), by=numt)) {
    .main = paste(t_vector[i]*10^15,"(fs)", seq="")
    plot(ka_vector, nvv_kt[,i], xlim=c(-pi,pi), ylim=nplim, xlab="BZ", ylab= "number/polarization", main=.main,col="green", type="l")
    par(new=T)
    plot(ka_vector, ncv_kt[,i], xlim=c(-pi,pi), ylim=nplim, xlab="", ylab= "", main=.main, col="blue", type="l")
    par(new=T)    
    plot(ka_vector, Re(mpv_kt[,i]), xlim=c(-pi,pi), ylim=nplim, xlab="BZ", ylab= "", main=.main,col="red" , type="l")
  }
}, movie.name="debug_np_velocity.gif", interval=0.2)
}



#======================================================================
#stacks 
if(0){
#debug comment out beginning
#write.table(number_electron_kt, "./number_electron_kt.txt")
#png("number.png", width = 700, height = 500)  
#plot(number_electron_kt[], xlab="", ylab="")  
#dev.off() 
#comment out end
}
#comment out
if(0){
w_2 <- rep(0, div_num_ka)
w_1 <- rep(0, div_num_ka)
w1 <- rep(0, div_num_ka)
w2 <- rep(0, div_num_ka)
ww <- rep(0, div_num_ka)
grad_k <- function(w0){
	w1[div_num_ka] <- w0[1]
	w1[1:(div_num_ka-1)] <- w0[2:div_num_ka]
	return ((w0-w1) / delta_k)
}
grad_k2 <- function(w0){
	w_2[1:2] <- w0[(div_num_ka-1):div_num_ka]
	w_2[3:div_num_ka] <- w0[1:(div_num_ka-2)]		
	w_1[1] <- w0[div_num_ka]
	w_1[2:div_num_ka] <- w0[1:(div_num_ka-1)]
	w1[div_num_ka] <- w0[1]
	w1[1:(div_num_ka-1)] <- w0[2:div_num_ka]
	w2[(div_num_ka-1):div_num_ka] <- w0[1:2]
	w2[1:(div_num_ka-2)] <- w0[3:div_num_ka]

	w_left <- 2*w1 + 3*w0  - 6*w_1 + w_2		
	w_right <- - w2 + 6*w1 - 3*w0  - 2*w_1
	ww <- w_left
	for(i in 1:div_num_ka){
		if(w_left[i] < 0)
		ww[i] <- w_left[i]
	}
	return (ww / 6 / delta_k)
}
#	ww <- replace( w_right, which(w_right > 0), w_left )  
grad_k3 <- function(E, w0){
	w_2[1:2] <- w0[(div_num_ka-1):div_num_ka]
	w_2[3:div_num_ka] <- w0[1:(div_num_ka-2)]		
	w_1[1] <- w0[div_num_ka]
	w_1[2:div_num_ka] <- w0[1:(div_num_ka-1)]
	w1[div_num_ka] <- w0[1]
	w1[1:(div_num_ka-1)] <- w0[2:div_num_ka]
	w2[(div_num_ka-1):div_num_ka] <- w0[1:2]
	w2[1:(div_num_ka-2)] <- w0[3:div_num_ka]

	ww <- (- w2 + 8*w1 - 8*w_1 + w_2	)/ 12 / delta_k
	return (ww)
}
}

#debug plot
if(0){
print("Exporting debug_band_velocity.gif ...")
numt <- round( 2*10^-15 /delta_t )
print("numt" )
print(numt)
erange <- c(-2,5) 
jj <- (elementary_charge /Plank_constant_bar * band_width_conduction* lattice_constant_a_of_GaSe)
jrange <- c(-jj,jj)
library(animation)
saveGIF({
  ani.options(loop = TRUE)
  for (i in seq(1,(div_num_t+1), by=numt)) {
    .main = paste(t_vector[i]*10^15,"(fs)", seq="")
    ka_vector <- ka_to_ka(ka_vector, A_vector[i])
    plot(ka_vector, energy_valence(ka_vector)/(1.60218*10^-19), xlim=c(-pi,pi), ylim=erange, xlab="BZ", ylab= "energy (eV)", main=.main, col="green", type="l")
    par(new=T)
    plot(ka_vector, energy_conduction(ka_vector)/(1.60218*10^-19), xlim=c(-pi,pi), ylim=erange, xlab="", ylab= "", main=.main, col="blue", type="l")  
    par(new=T)
    plot(ka_vector, micro_current_valence(ka_vector), xlim=c(-pi,pi), ylim=jrange, xlab="BZ", ylab= "energy (eV)", main=.main, col="green", type="p")
    par(new=T)
    plot(ka_vector, micro_current_conduction(ka_vector), xlim=c(-pi,pi), ylim=jrange, xlab="", ylab= "", main=.main, col="blue", type="p")  
    }
}, movie.name="debug_band_velocity.gif", interval=0.2)
}

