#I use SI unit system for calculation
# first written in 6 Sep. 2017
# modified in 11 Sep. 2017
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
range_t <- 250 *10^(-15) # fs
div_num_t <- 5000 # n
delta_t <-  range_t/ div_num_t
t_vector <- seq( -range_t/2, range_t/2, length= (div_num_t+1))
t_vector2 <- seq( -range_t/2, range_t/2, length= (2*div_num_t+1))
t_vector4 <- seq( -range_t/2, range_t/2, length= (4*div_num_t+1))
frequency <- seq(0, 1/delta_t, length= (div_num_t+1))

# External electric field (assuming the gaussian envelope)
E_amplitude_max <- 10* 10^8 # V/m
envelope_width <- 50* 10^(-15) # 1/e^2 width; 50 fs 
carrier_frequency <- 33 *10^12 # 33 THz
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
	A_vector[i+1] <- A_vector[i] + (E_vector2[2*i+1]+4*E_vector2[2*i]+E_vector2[2*1-1])*delta_t/6
}
A_vector2 <- rep(0, 2*div_num_t+1)
for(i in 1:(2*div_num_t)){
	A_vector2[i+1] <- A_vector2[i] + (E_vector4[2*i+1]+4*E_vector4[2*i]+E_vector4[2*1-1])*delta_t/6
}
#debug plot
if(0){
print("Exporting debug_ElectricF_VectorP.png ...")
png("./debug_ElectricF_VectorP.png", width = 700, height = 500)  
plot(t_vector*10^15, E_vector, xlab="t (fs)", ylab= "J", col="red", type="l")  
par(new=T)
plot(t_vector*10^15, A_vector, xlab="", ylab="", col="blue", type="l")  
dev.off() 
}	

# fixed k normalized by inverse lattice constant 1/a in length gauge
div_num_fka <- 150 # m
delta_fk <- 2*pi / lattice_constant_a_of_GaSe / div_num_fka 
fka_vector_redundant <- seq( -pi, pi, length= (div_num_fka+1)) #BZの両端で冗長
fka_vector <- fka_vector_redundant[1:div_num_fka]
fka_to_ka <- function(fka, A){
	ka <-  (fka+ A* elementary_charge* lattice_constant_a_of_GaSe/Plank_constant_bar ) %% (2*pi) #modulo 2pi 
} 

#band structure of GaSe (two band model assuming the Tight binding model (1D))
band_width_valence <- 1 *1.60218*10^-19 # J (1eV)
band_width_conduction <- 1.5 *1.60218*10^-19 # J (1.5eV) 
band_gap <- 2 *1.60218*10^-19 # J (2eV) 
energy_valence <- function(ka){
	band_width_valence* (cos( ka ) -1) /2 # J
}
energy_conduction <- function(ka){
	band_width_conduction* (1 - cos( ka ) )/2  +  band_gap # J 
}
energy_difference <- function(ka){
	energy_conduction(ka) - energy_valence(ka) # J 
}


#microscopic current
micro_current_valence <- function(ka){
	-elementary_charge / Plank_constant_bar * band_width_valence *lattice_constant_a_of_GaSe *(-sin(ka)) /2
}
micro_current_conduction <- function(ka){
	 -elementary_charge /Plank_constant_bar * band_width_conduction *lattice_constant_a_of_GaSe *sin(ka)/2
}

#dipole transition matrix
Rabi_energy_max <- 12 *1.60218*10^-19 #1.2eV; converstion from eV to J
dipole_transition_max <- Rabi_energy_max / E_amplitude_max 
ed_length <- energy_difference(fka_vector)

dipole_transition_length <- rep(dipole_transition_max, div_num_fka)  #/ ed_length
dipole_transition_velocity <- rep(dipole_transition_max, div_num_fka)
if(0){
print("Exporting debug_dipole_transition.png ...")
png("./debug_dipole_transition.png", width = 700, height = 500)  
plot(fka_vector, dipole_transition_length , xlab="ka (rad)", ylab= "dipole_transition", type="l")  
dev.off()   
}

#parameters of calculation
Courant_number <- delta_t * E_amplitude_max *  elementary_charge / Plank_constant_bar / delta_fk 
print("Courant number (length gauge): " )
print( Courant_number )

w_2 <- rep(0, div_num_fka)
w_1 <- rep(0, div_num_fka)
w1 <- rep(0, div_num_fka)
w2 <- rep(0, div_num_fka)
ww <- rep(0, div_num_fka)

grad_midpoint <- function(w0){
	w_1[1] <- w0[div_num_fka]
	w_1[2:div_num_fka] <- w0[1:(div_num_fka-1)]
	w1[div_num_fka] <- w0[1]
	w1[1:(div_num_fka-1)] <- w0[2:div_num_fka]
	return ((w1-w_1)/ 2 / delta_fk)
}

#comment out
if(0){
grad_k1 <- function(w0){
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

#relaxation
relaxation_constant_length <-  50 *10^(-15) # 7fs 
#dephasing 
dephasing_constant_length <- 7 *10^(-15) # 1.1fs 
diff_number_valence_length <- function(E, mpl, nvl){
	dnvl <- ( 2*Im(E*dipole_transition_length*mpl) + elementary_charge* E* grad_midpoint(nvl) )  * delta_t/ Plank_constant_bar
	dnvl <- dnvl - (nvl-1)* relaxation_constant_length
	return (dnvl)
}
diff_number_conduction_length <- function(E, mpl, ncl){
	dncl <- ( 2*Im(E*dipole_transition_length*mpl) + elementary_charge* E* grad_midpoint(ncl) )* delta_t/ Plank_constant_bar
	dncl <- dncl - ncl / relaxation_constant_length 
	return (dncl)
}
diff_micro_polarization_length <- function(E, mpl, nvl, ncl){
	dmpl <- ( -1i *ed_length *mpl + elementary_charge* E* grad_midpoint(mpl) + 1i* E* dipole_transition_length * (nvl- ncl) ) * delta_t/ Plank_constant_bar
	#dmpl <- dmpl -mpl / dephasing_constant_length * delta_t
	return (dmpl)
}

nvl_kt <- array(1, dim=c(div_num_fka, (div_num_t+1))) #number_valence_length
ncl_kt <- array(0, dim=c(div_num_fka, (div_num_t+1))) #number_conduction_length
mpl_kt  <- array(0, dim=c(div_num_fka, (div_num_t+1))) #micro_polarization_length

nvl <- rep(1, div_num_fka)
ncl <- rep(0, div_num_fka) 
mpl <- rep(0, div_num_fka)

Current_valence_length <- rep(0, (div_num_t+1))
Current_conduction_length <- rep(0, (div_num_t+1))
Polarization_length <- rep(0, (div_num_t+1))
Displacement_current_length <- rep(0, (div_num_t))

#solving Lioville equation with the length gauge and bloch-state representation (SBE) by the 4th Runge-Kutta method
if(0){
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

mcv_length <- micro_current_valence(fka_vector)
mcc_length <- micro_current_conduction(fka_vector)
for(i in 1:(div_num_fka)){
	Current_valence_length <-  Current_valence_length + mcv_length[i] * nvl_kt[i,]
	Current_conduction_length <- Current_conduction_length + mcc_length[i] * ncl_kt[i,]
	Polarization_length <- Polarization_length + mpl_kt[i,] * dipole_transition_length[i]
}

Current_length <- Current_valence_length + Current_conduction_length
Displacement_current_length_short <- (Polarization_length[2:(div_num_t+1)] - Polarization_length[1:div_num_t]) /delta_t
Displacement_current_length <- append( Displacement_current_length_short, c(0.0))
E_HHG_length <- Current_length + Displacement_current_length_short

E_spectrum <- fft(E_vector)
E_HHG_length_spectrum <- fft(E_HHG_length)
Current_length_spectrum <-  fft(Current_length)
Displacement_current_length_spectrum <- fft(Displacement_current_length)

if(0){
print("Exporting debug_np_length.gif ...")
numt <- round( 2*10^-15 /delta_t )
nplim=c(-1.2, 1.2)
library(animation)
saveGIF({
  ani.options(loop = TRUE)
  for (i in seq(1,(div_num_t+1), by=numt)) {
    .main = paste(t_vector[i]*10^15,"(fs)", seq="")
    plot(fka_vector, nvl_kt[,i], xlim=c(-pi,pi), ylim=nplim, xlab="BZ", ylab= "numbers", main=.main,col="green", type="l")
    par(new=T)
    plot(fka_vector, ncl_kt[,i], xlim=c(-pi,pi), ylim=nplim, xlab="", ylab= "", main=.main, col="blue", type="l")
    par(new=T)    
    plot(fka_vector, Re(mpl_kt[,i]), xlim=c(-pi,pi), ylim=nplim, xlab="BZ", ylab= "polarization", main=.main,col="red" , type="l")
  }
}, movie.name="debug_np_length.gif", interval=0.2)
}

E_HHGrange <- c( -max(abs(E_HHG_length)), max(abs(E_HHG_length)) ) 
Erange <- c( -max(E_vector), max(E_vector) ) 
trange <- c(-range_t/2, range_t/2)*10^15
png("./Bloch_lengthG/debug_ElectricF_length.png", width = 700, height = 500)  
plot(t_vector*10^15, Re(E_vector), xlab="", ylab= "Electric field (arb. unit)", xlim=trange, ylim=Erange, col="black", type="l")  
par(new=T)
plot(t_vector*10^15, Re(E_HHG_length), xlab="t (fs)", ylab="", xlim=trange,  ylim=E_HHGrange, col="red", type="l")  
par(new=T)
plot(t_vector*10^15, Re(Displacement_current_length), xlab="", ylab= "", xlim=trange, ylim=E_HHGrange, col="green", type="l")  
par(new=T)
plot(t_vector*10^15, Re(Current_length), xlab="", ylab="", xlim=trange, ylim=E_HHGrange, col="blue", type="l")  
dev.off()    

E_HHGrange <- c( min(abs(E_HHG_length_spectrum)^2), max(abs(E_HHG_length_spectrum)^2) ) 
frange <- c(0, 20)
png("./Bloch_lengthG/debug_spectrum_length.png", width = 700, height = 500) 
plot(frequency/carrier_frequency, abs(E_spectrum)^2, xlim=frange, log="y", xlab="f/carrier frequency", ylab="power (arb. unit)",col="black", type="l")  
par(new=T)
plot(frequency/carrier_frequency, abs(E_HHG_length_spectrum)^2, xlim=frange, ylim=E_HHGrange, log="y", axes=FALSE, yaxt="n",xlab="", ylab= "", col="red", type="l")  
par(new=T)
plot(frequency/carrier_frequency, abs(Current_length_spectrum)^2, xlim=frange, ylim=E_HHGrange, log="y", axes=FALSE, yaxt="n",xlab="", ylab= "", col="green", type="l")  
par(new=T)
plot(frequency/carrier_frequency, abs(Displacement_current_length_spectrum)^2, xlim=frange, ylim=E_HHGrange, log="y", axes=FALSE, yaxt="n",xlab="", ylab= "", col="blue", type="l")  
abline(v=0:20, col='black', lty="dotted")
dev.off()     
}

#---------------------------------------------------------------------------------------------------------------

#smoothing time
smoothing_time_velocity <- 0.5*10^-15 #fs
#relxation time
relaxation_constant_velocity <-7*10^-15 #fs
xv1 <- rep(0, div_num_fka)
xv_1 <- rep(0, div_num_fka)
curvature <- function(xv){
	xv1[1] <- xv[div_num_fka]
	xv1[2:div_num_fka] <- xv[1:(div_num_fka-1)]
	xv_1[div_num_fka] <- xv[1]
	xv_1[1:(div_num_fka-1)] <- xv[2:div_num_fka]
  	(xv1- 2*xv + xv_1)
}
diff_number_valence_velocity <- function(A, ka, mpv){
	dnvv <- 2*A *Re(dipole_transition_velocity*mpv) *energy_difference(ka) *delta_t / Plank_constant_bar^2
	#dnvv <- dnvv - (nvv-1) /relaxation_constant_velocity *delta_t
	dnvv <- dnvv + curvature(nvv) / smoothing_time_velocity *delta_t
	return (dnvv)
}
diff_number_conduction_velocity <- function(A, ka, mpv){
	dncv <- -2*A *Re(dipole_transition_velocity*mpv) *energy_difference(ka) *delta_t / Plank_constant_bar^2 
	#dncv <- dncv - ncv /relaxation_constant_velocity *delta_t
	dncv <- dncv + curvature(ncv) / smoothing_time_velocity *delta_t
	return (dncv)
}
diff_micro_polarization_velocity <- function(A, ka, mpv, nvv, ncv){
	dmpv <- ( -1i /Plank_constant_bar *energy_difference(ka) *mpv 
	- A* dipole_transition_velocity * (nvv-ncv) *energy_difference(ka) /Plank_constant_bar^2 ) * delta_t
	#dmpv <- dmpv -mpv / dephasing_constant * delta_t #dephasing
	#dmpv <- dmpv + curvature(mpv) * delta_t / smoothing_time
	return (dmpv)
}

nvv_kt <- array(1, dim=c(div_num_fka, (div_num_t+1))) #number_valence_velocity
ncv_kt <- array(0, dim=c(div_num_fka, (div_num_t+1))) #number_conduction_velocity
mpv_kt  <- array(0, dim=c(div_num_fka, (div_num_t+1))) #micro_polarization_velocity

nvv <- rep(0.55, div_num_fka)
ncv <- rep(0.45, div_num_fka) 
mpv <- rep(0, div_num_fka)

Current_valence_velocity <- rep(0, (div_num_t+1))
Current_conduction_velocity <- rep(0, (div_num_t+1))
Polarization_velocity <- rep(0, (div_num_t+1))
Displacement_current_velocity <- rep(0, (div_num_t))

#solving Lioville equation with the velocity gauge and bloch-state representation (SBE) by the 4th Runge-Kutta method
if(1){
print("starting calculation with velocity gauge")
for( j in 1:(div_num_t) ){
	A <- A_vector2[2*j-1]
	A_0 <- A_vector2[2*j]
	A_1 <- A_vector2[2*j]	
	A_2 <- A_vector2[2*j+1]	
	
	ka <- fka_to_ka(fka_vector, A)
	ka_0 <- fka_to_ka(fka_vector, A_0)
	ka_1 <- fka_to_ka(fka_vector, A_1)
	ka_2 <- fka_to_ka(fka_vector, A_2)
	
	dnvv_0 <- diff_number_valence_velocity(A, ka, mpv)
	dncv_0 <- diff_number_conduction_velocity(A, ka, mpv)
	dmpv_0 <- diff_micro_polarization_velocity(A, ka, mpv, nvv, ncv)

	dnvv_1 <- diff_number_valence_velocity(A_0, ka_0, mpv+dmpv_0/2)
	dncv_1 <- diff_number_conduction_velocity(A_0, ka_0, mpv+dmpv_0/2)
	dmpv_1 <- diff_micro_polarization_velocity(A_0, ka_0, mpv+dmpv_0/2, nvv+dnvv_0/2, ncv+dncv_0/2)
	
	dnvv_2 <- diff_number_valence_velocity(A_1, ka_1, mpv+dmpv_1/2)
	dncv_2 <- diff_number_conduction_velocity(A_1, ka_1, mpv+dmpv_1/2)
	dmpv_2 <- diff_micro_polarization_velocity(A_1, ka_1, mpv+dmpv_1/2, nvv+dnvv_1/2, ncv+dncv_1/2)

	dnvv_3 <- diff_number_valence_velocity(A_2, ka_2, mpv+dmpv_2)
	dncv_3 <- diff_number_conduction_velocity(A_2, ka_2, mpv+dmpv_2)
	dmpv_3 <- diff_micro_polarization_velocity(A_2, ka_2, mpv+dmpv_2, nvv+dnvv_2, ncv+dncv_2)

	nvv <- nvv + (dnvv_0 +2*dnvv_1 +2*dnvv_2 +dnvv_3)/6
	ncv <- ncv + (dncv_0 +2*dncv_1 +2*dncv_2 +dncv_3)/6
	mpv <- mpv + (dmpv_0 +2*dmpv_1 +2*dmpv_2 +dmpv_3)/6

	nvv_kt[,j+1] <- nvv
	ncv_kt[,j+1] <- ncv
	mpv_kt[,j+1] <- mpv
	
	mcv_velocity <- micro_current_valence(ka)
	mcc_velocity <- micro_current_conduction(ka)
	k <- j+1
	for(i in 1:(div_num_fka)){
		
		Current_valence_velocity[k] <-  Current_valence_velocity[k] + mcv_velocity[i] * nvv_kt[i,k]
		Current_conduction_velocity[k] <- Current_conduction_velocity[k] + mcc_velocity[i] * ncv_kt[i,k]
		Polarization_velocity[k] <- Polarization_velocity[k] + mpv_kt[i,k] * dipole_transition_velocity[i]
	}	
}

Current_velocity <- Current_valence_velocity + Current_conduction_velocity
Displacement_current_velocity_short <- (Polarization_velocity[2:(div_num_t+1)] - Polarization_velocity[1:div_num_t]) /delta_t
Displacement_current_velocity <- append( Displacement_current_velocity_short, c(0.0))
E_HHG_velocity <- Current_velocity + Displacement_current_velocity_short

E_spectrum <- fft(E_vector)
E_HHG_velocity_spectrum <- fft(E_HHG_velocity)
Current_velocity_spectrum <-  fft(Current_velocity)
Displacement_current_velocity_spectrum <- fft(Displacement_current_velocity)

if(1){
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
    plot(fka_vector, nvv_kt[,i], xlim=c(-pi,pi), ylim=nplim, xlab="BZ", ylab= "number/polarization", main=.main,col="green", type="l")
    par(new=T)
    plot(fka_vector, ncv_kt[,i], xlim=c(-pi,pi), ylim=nplim, xlab="", ylab= "", main=.main, col="blue", type="l")
    par(new=T)    
    plot(fka_vector, Re(mpv_kt[,i]), xlim=c(-pi,pi), ylim=nplim, xlab="BZ", ylab= "", main=.main,col="red" , type="l")
  }
}, movie.name="debug_np_velocity.gif", interval=0.2)
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
    ka_vector <- fka_to_ka(fka_vector, A_vector[i])
    plot(fka_vector, energy_valence(ka_vector)/(1.60218*10^-19), xlim=c(-pi,pi), ylim=erange, xlab="BZ", ylab= "energy (eV)", main=.main, col="green", type="l")
    par(new=T)
    plot(fka_vector, energy_conduction(ka_vector)/(1.60218*10^-19), xlim=c(-pi,pi), ylim=erange, xlab="", ylab= "", main=.main, col="blue", type="l")  
    par(new=T)
    plot(fka_vector, micro_current_valence(ka_vector), xlim=c(-pi,pi), ylim=jrange, xlab="BZ", ylab= "energy (eV)", main=.main, col="green", type="p")
    par(new=T)
    plot(fka_vector, micro_current_conduction(ka_vector), xlim=c(-pi,pi), ylim=jrange, xlab="", ylab= "", main=.main, col="blue", type="p")  
    }
}, movie.name="debug_band_velocity.gif", interval=0.2)
}

E_HHGrange <- c( -max(abs(E_HHG_velocity)), max(abs(E_HHG_velocity)) ) 
Erange <- c( -max(E_vector), max(E_vector) ) 
trange <- c(-range_t/2, range_t/2)*10^15
png("./Bloch_velocityG/debug_ElectricF_velocity.png", width = 700, height = 500)  
plot(t_vector*10^15, Re(E_vector), xlab="", ylab= "Electric field (arb. unit)", xlim=trange, ylim=Erange, col="black", type="l")  
par(new=T)
plot(t_vector*10^15, Re(E_HHG_velocity), xlab="t (fs)", ylab="", xlim=trange,  ylim=E_HHGrange, col="red", type="l")  
par(new=T)
plot(t_vector*10^15, Re(Displacement_current_velocity), xlab="", ylab= "", xlim=trange, ylim=E_HHGrange, col="green", type="l")  
par(new=T)
plot(t_vector*10^15, Re(Current_velocity), xlab="", ylab="", xlim=trange, ylim=E_HHGrange, col="blue", type="l")  
dev.off()    

E_HHGrange <- c( min(abs(E_HHG_velocity_spectrum)^2), max(abs(E_HHG_velocity_spectrum)^2) ) 
frange <- c(0, 20)
png("./Bloch_velocityG/debug_spectrum_velocity.png", width = 700, height = 500) 
plot(frequency/carrier_frequency, abs(E_spectrum)^2, xlim=frange, log="y", xlab="f/carrier frequency", ylab="power (arb. unit)",col="black", type="l")  
par(new=T)
plot(frequency/carrier_frequency, abs(E_HHG_velocity_spectrum)^2, xlim=frange, ylim=E_HHGrange, log="y", axes=FALSE, yaxt="n",xlab="", ylab= "", col="red", type="l")  
par(new=T)
plot(frequency/carrier_frequency, abs(Current_velocity_spectrum)^2, xlim=frange, ylim=E_HHGrange, log="y", axes=FALSE, yaxt="n",xlab="", ylab= "", col="green", type="l")  
par(new=T)
plot(frequency/carrier_frequency, abs(Displacement_current_velocity_spectrum)^2, xlim=frange, ylim=E_HHGrange, log="y", axes=FALSE, yaxt="n",xlab="", ylab= "", col="blue", type="l")  
abline(v=0:20, col='black', lty="dotted")
dev.off()     

}

if(0){
#debug comment out beginning
#write.table(number_electron_kt, "./number_electron_kt.txt")
#png("number.png", width = 700, height = 500)  
#plot(number_electron_kt[], xlab="", ylab="")  
#dev.off() 
#comment out end
}

