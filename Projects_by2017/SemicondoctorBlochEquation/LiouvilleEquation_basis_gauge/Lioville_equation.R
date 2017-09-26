#I use SI unit system for calculation
# first written in 6 Sep. 2017
# modified in 11 Sep. 2017
# modified for basis-gauge independent calculation
print("hoge")
options(scipen=0)

#Physical constants 
Plank_constant_bar <- 6.62607004 * 10^-34 #m^2 kg / s
Boltzmann_constant <-  1.38064852 *10^-23 # m^2 kg s^-2 K^-1
speed_of_light  <- 299792458 # m/s
elementary_charge <- 1.60217662 * 10^-19 # Coulombs

#Parameters of experiment / material properties
Temperature <- 300 #K
Fermi_energy <- # Chemical potential
lattice_constant_a_of_GaSe <- 3.755 * 2/sqrt(3) *10^-10# m (3.755 ångström)
print("hoge")
# Time t
range_t <- 300 *10^(-15) # fs
div_num_t <- 10000 # n
delta_t <-  range_t/ div_num_t
t_vector <- seq( -range_t/2, range_t/2, length= (div_num_t+1))
t_vector2 <- seq( -range_t/2, range_t/2, length= (2*div_num_t+1))
t_vector4 <- seq( -range_t/2, range_t/2, length= (4*div_num_t+1))
frequency <- seq(0, 1/delta_t, length= (div_num_t+1))

# External electric field (assuming the gaussian envelope)
E_amplitude_max <- 30* 10^8 # V/m
envelope_width <- 50* 10^(-15) # 1/e^2 width; 50 fs 
carrier_frequency <- 33 *10^12 # 33 THz
carrier_envelope_phase <- 0 # [rad]
Electric_field <- function(t){
	EF <- E_amplitude_max*exp(-t^2/envelope_width^2) *cos(2*pi* carrier_frequency*t +carrier_envelope_phase) 
}
E_vector <-  Electric_field(t_vector)
E_vector2 <-  Electric_field(t_vector2)
E_vector4 <-  Electric_field(t_vector4)
print("hoge")
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
png("debug_E_A.png", width = 700, height = 500)  
plot(t_vector*10^15, E_vector, xlab="t (fs)", ylab= "J", col="red", type="l")  
par(new=T)
plot(t_vector*10^15, A_vector, xlab="", ylab="", col="blue", type="l")  
dev.off() 
}	
print("hoge")
# fixed k normalized by inverse lattice constant 1/a in length gauge
div_num_fka <- 150 # m
delta_fk <- 2*pi / lattice_constant_a_of_GaSe / div_num_fka 
fka_vector_redundant <- seq( -pi, pi, length= (div_num_fka+1)) #BZの両端で冗長
fka_vector <- fka_vector_redundant[1:div_num_fka]
fka_to_ka <- function(fka, A){
	ka <-  (fka+ A* elementary_charge* lattice_constant_a_of_GaSe/Plank_constant_bar ) %% (2*pi) #modulo 2pi 
} 
print("hoge")
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


print("hoge")

#debug plot
if(0){
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
}, movie.name="debug_band.gif", interval=0.2)
}

print("hoge")

#dipole transition matrix
Rabi_energy_max <- 0.12 *1.60218*10^-19 #1.2eV; converstion from eV to J
dipole_transition_max <- Rabi_energy_max / E_amplitude_max 
ed_length <- energy_difference(fka_vector)
dipole_transition_length <- dipole_transition_max * min(ed_length) / ed_length
if(0){
png("debug_dipole_transition.png", width = 700, height = 500)  
plot(fka_vector, dipole_transition_length , xlab="ka (rad)", ylab= "dipole_transition", type="l")  
dev.off()   
}

#relaxation
#relaxation_constant <-  50 *10^(-15) # 7fs 
#dephasing 
dephasing_constant <- 5 *10^(-15) # 1.1fs 

#parameters of calculation
Courant_number <- delta_t * E_amplitude_max *  elementary_charge / Plank_constant_bar / delta_fk 
print("Courant number (length gauge): " )
print( Courant_number )

w_2 <- rep(0, div_num_fka)
w_1 <- rep(0, div_num_fka)
w1 <- rep(0, div_num_fka)
w2 <- rep(0, div_num_fka)
ww <- rep(0, div_num_fka)

print("hoge")

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

diff_number_valence_length <- function(E, mpl, nvl){
	dnvl <- ( 2*Im(E*dipole_transition_length*mpl) + elementary_charge* E* grad_midpoint(nvl) )  * delta_t/ Plank_constant_bar
	#- nvl* relaxation_constant ) 
	return (dnvl)
}
diff_number_conduction_length <- function(E, mpl, ncl){
	dncl <- ( 2*Im(E*dipole_transition_length*mpl) + elementary_charge* E* grad_midpoint(ncl) )* delta_t/ Plank_constant_bar
	#- ncl / relaxation_constant ) 
	return (dncl)
}
diff_micro_polarization_length <- function(E, mpl, nvl, ncl){
	dmpl <- ( -1i *ed_length *mpl + elementary_charge* E* grad_midpoint(mpl) + 1i* E* dipole_transition_length * (nvl- ncl) ) * delta_t/ Plank_constant_bar
	dmpl <- dmpl -mpl / dephasing_constant * delta_t
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

	dnvl_3 <- diff_number_valence_length(E_2, mpl+dmpl_2/2, nvl+dnvl_2/2)
	dncl_3 <- diff_number_conduction_length(E_2, mpl+dmpl_2/2, ncl+dncl_2/2)
	dmpl_3 <- diff_micro_polarization_length(E_2, mpl+dmpl_2/2, nvl+dnvl_2/2, ncl+dncl_2/2)

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
E_HHG <- Current_length + Displacement_current_length_short

E_spectrum <- fft(E_vector)
E_HHG_spectrum <- fft(E_HHG)
Current_length_spectrum <-  fft(Current_length)
Displacement_current_length_spectrum <- fft(Displacement_current_length)

if(0){
numt <- round( 2*10^-15 /delta_t )
print("numt" )
print(numt)
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
}, movie.name="number_polarization.gif", interval=0.2)
}

E_HHGrange <- c( -max(abs(E_HHG)), max(abs(E_HHG)) ) 
Erange <- c( -max(E_vector), max(E_vector) ) 
trange <- c(-range_t/2, range_t/2)*10^15
png("E.png", width = 700, height = 500)  
plot(t_vector*10^15, Re(E_vector), xlab="", ylab= "Electric field (arb. unit)", xlim=trange, ylim=Erange, col="black", type="l")  
par(new=T)
plot(t_vector*10^15, Re(E_HHG), xlab="t (fs)", ylab="", xlim=trange,  ylim=E_HHGrange, col="red", type="l")  
par(new=T)
plot(t_vector*10^15, Re(Displacement_current_length), xlab="", ylab= "", xlim=trange, ylim=E_HHGrange, col="green", type="l")  
par(new=T)
plot(t_vector*10^15, Re(Current_length), xlab="", ylab="", xlim=trange, ylim=E_HHGrange, col="blue", type="l")  
dev.off()    

E_HHGrange <- c( min(abs(E_HHG_spectrum)^2), max(abs(E_HHG_spectrum)^2) ) 
frange <- c(0, 20)
png("debug_spectrum_length.png", width = 700, height = 500) 
plot(frequency/carrier_frequency, abs(E_spectrum)^2, xlim=frange, log="y", xlab="f/carrier frequency", ylab="power (arb. unit)",col="black", type="l")  
par(new=T)
plot(frequency/carrier_frequency, abs(E_HHG_spectrum)^2, xlim=frange, ylim=E_HHGrange, log="y", axes=FALSE, yaxt="n",xlab="", ylab= "", col="red", type="l")  
par(new=T)
plot(frequency/carrier_frequency, abs(Current_length_spectrum)^2, xlim=frange, ylim=E_HHGrange, log="y", axes=FALSE, yaxt="n",xlab="", ylab= "", col="green", type="l")  
par(new=T)
plot(frequency/carrier_frequency, abs(Displacement_current_length_spectrum)^2, xlim=frange, ylim=E_HHGrange, log="y", axes=FALSE, yaxt="n",xlab="", ylab= "", col="blue", type="l")  
abline(v=0:20, col='black', lty="dotted")
dev.off()     

---------------------------------------------------------------------------------------------------------------

diff_number_valence_length <- function(E, mpl, nvl){
	dnvl <- ( 2*Im(E*dipole_transition_length*mpl) + elementary_charge* E* grad_midpoint(nvl) )  * delta_t/ Plank_constant_bar
	#- nvl* relaxation_constant ) 
	return (dnvl)
}
diff_number_conduction_length <- function(E, mpl, ncl){
	dncl <- ( 2*Im(E*dipole_transition_length*mpl) + elementary_charge* E* grad_midpoint(ncl) )* delta_t/ Plank_constant_bar
	#- ncl / relaxation_constant ) 
	return (dncl)
}
diff_micro_polarization_length <- function(E, mpl, nvl, ncl){
	dmpl <- ( -1i *ed_length *mpl + elementary_charge* E* grad_midpoint(mpl) + 1i* E* dipole_transition_length * (nvl- ncl) ) * delta_t/ Plank_constant_bar
	dmpl <- dmpl -mpl / dephasing_constant * delta_t
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

	dnvl_3 <- diff_number_valence_length(E_2, mpl+dmpl_2/2, nvl+dnvl_2/2)
	dncl_3 <- diff_number_conduction_length(E_2, mpl+dmpl_2/2, ncl+dncl_2/2)
	dmpl_3 <- diff_micro_polarization_length(E_2, mpl+dmpl_2/2, nvl+dnvl_2/2, ncl+dncl_2/2)

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
E_HHG <- Current_length + Displacement_current_length_short

E_spectrum <- fft(E_vector)
E_HHG_spectrum <- fft(E_HHG)
Current_length_spectrum <-  fft(Current_length)
Displacement_current_length_spectrum <- fft(Displacement_current_length)


if(0){
png("HHG_para.png", width = 500, height = 500)  
filled.contour(CEP/pi, frequency/ carrier_frequency, -log10(HHG_CEP), ylim=c(0,10), zlim=c(-15,-5),col=rainbow(256),nlevels=(256))
dev.off() 
}

warnings()

if(1){
numt <- round( 2*10^-15 /delta_t )
print("numt" )
print(numt)
nplim=c(-0.7, 0.7)
library(animation)
saveGIF({
  ani.options(loop = TRUE)
  for (i in seq(1,(div_num_t+1), by=numt)) {
    .main = paste(t_vector[i]*10^15,"(fs)", seq="")
    plot(ka_vector, number_electron_kt[,i], xlim=c(-pi,pi), ylim=nplim, xlab="BZ", ylab= "numbers", main=.main,col="green", type="l")
    par(new=T)
    plot(ka_vector, number_hole_kt[,i], xlim=c(-pi,pi), ylim=nplim, xlab="", ylab= "", main=.main, col="blue", type="l")
    par(new=T)    
    plot(ka_vector, Re(polarization_kt[,i]), xlim=c(-pi,pi), ylim=nplim, xlab="BZ", ylab= "polarization", main=.main,col="red" , type="l")
  }
}, movie.name="number_polarization.gif", interval=0.2)
}

if(0){
#debug comment out beginning


Jrange <- 3*c( min(Re(J)), max(Re(J)))
png("J.png", width = 700, height = 500)  
plot(t_vector*10^15, Re(J), xlab="t (fs)", ylab= "J", ylim=Jrange, col="red", type="l")  
par(new=T)
plot(t_vector*10^15, Re(J_electron), xlab="", ylab= "", ylim=Jrange, col="green", type="l")  
par(new=T)
plot(t_vector*10^15, Re(J_hole), xlab="", ylab="", ylim=Jrange, col="blue", type="l")  
dev.off()    

Prange <- c( -max(abs(P)), max(abs(P)))
png("P.png", width = 700, height = 500)  
plot(t_vector*10^15, Re(P), xlab="t (fs)", ylab= "P", ylim=Prange, col="red", type="l")  
par(new=T)
plot(t_vector*10^15, Im(P), xlab="", ylab="", ylim=Prange, col="blue", type="l")  
dev.off()       

E_HHGrange <- c( -max(abs(E_HHG)), max(abs(E_HHG)) ) 
Erange <- c( -max(E_vector), max(E_vector) ) 
trange <- c(-100, 100)
png("E.png", width = 700, height = 500)  
plot(t_vector*10^15, Re(E_HHG), xlab="t (fs)", ylab="Electric field (arb. unit)", xlim=trange,  ylim=E_HHGrange, axes=FALSE, yaxt="n", col="red", type="l")  
par(new=T)
plot(t_vector*10^15, Re(E_vector), xlab="", ylab= "", xlim=trange, ylim=Erange, col="black", type="l")  
#par(new=T)
#plot(t_vector*10^15, Re(dP_dt), xlab="", ylab= "", ylim=Erange, col="green", type="l")  
#par(new=T)
#plot(t_vector*10^15, Re(J), xlab="", ylab="", ylim=Erange, col="blue", type="l")  
dev.off()       


#write.table(number_electron_kt, "./number_electron_kt.txt")


#png("number.png", width = 700, height = 500)  
#plot(number_electron_kt[], xlab="", ylab="")  
#dev.off() 

frange <- c(0, 20)
Erange <- 10^c(max(2*log10(abs(E_spectrum)))-15, max(2*log10(abs(E_spectrum)))+3)
E_HHGrange <- 10^c(max(2*log10(abs(E_HHG_spectrum)))-13, max(2*log10(abs(E_HHG_spectrum)))+5)
#Jrange <- 10^c(max(2*log10(abs(J_spectrum)))-13, max(2*log10(abs(J_spectrum)))+5)
#dP_dtrange <- 10^c(max(2*log10(abs(dP_dt_spectrum)))-13, max(2*log10(abs(dP_dt_spectrum)))+5)

png("E_spectrum.png", width = 700, height = 500) 
plot(frequency/carrier_frequency, abs(E_spectrum)^2, xlim=frange, ylim=Erange, log="y", xlab="f/carrier frequency", ylab="power (arb. unit)",col="black", type="l")  
par(new=T)
plot(frequency/carrier_frequency, abs(E_HHG_spectrum)^2, xlim=frange, ylim=E_HHGrange, log="y", axes=FALSE, yaxt="n",xlab="", ylab= "", col="red", type="l")  
par(new=T)
plot(frequency/carrier_frequency, abs(J_spectrum)^2, xlim=frange, ylim=E_HHGrange, log="y", axes=FALSE, yaxt="n",xlab="", ylab= "", col="green", type="l")  
par(new=T)
plot(frequency/carrier_frequency, abs(dP_dt_spectrum)^2, xlim=frange, ylim=E_HHGrange, log="y", axes=FALSE, yaxt="n",xlab="", ylab= "", col="blue", type="l")  
abline(v=0:20, col='black', lty="dotted")
dev.off()     

max(abs(E_HHG_spectrum)^2)

#comment out end
}

