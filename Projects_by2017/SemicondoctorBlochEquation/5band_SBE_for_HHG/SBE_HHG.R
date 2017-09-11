#I use SI unit system for calculation
# first written in 6 Sep. 2017
#modified in 11 Sep. 2017

options(scipen=300)

#Physical constants 
Plank_constant <- 6.62607004 * 10^-34 #m^2 kg / s
Boltzmann_constant <-  1.38064852 *10^-23 # m^2 kg s^-2 K^-1
speed_of_light  <- 299792458 # m/s
elementary_charge <- 1.60217662 * 10^-19 # Coulombs

#Parameters of experiment / material properties
Temperature <- 300 #K
Fermi_energy <- # Chemical potential
lattice_constant_a_of_GaSe <- 3.755 *10^-10# m (3.755 ångström)

# Time t
range_t <- 200 *10^(-15) # fs
div_num_t <- 3000 # n
delta_t <-  range_t/ div_num_t
t_vector <- seq( -range_t/2, range_t/2, length= (div_num_t+1))
t_vector2 <- seq( -range_t/2, range_t/2, length= (2*div_num_t+1))
frequency <- seq(0, 1/delta_t, length= (div_num_t+1))

# External electric field (assuming the gaussian envelope)
E_amplitude_max <- 0.1* 10^8 # V/m
envelope_width <- 50* 10^(-15) # 1/e^2 width; 50 fs 
carrier_frequency <- 33 *10^12 # 33 THz
carrier_envelope_phase <- 0 # [rad]
E_vector <-  E_amplitude_max * exp(-t_vector^2/envelope_width^2) *cos(2*pi* carrier_frequency* t_vector+carrier_envelope_phase) 
E_vector2 <-  E_amplitude_max * exp(-t_vector2^2/envelope_width^2) *cos(2*pi* carrier_frequency* t_vector2+carrier_envelope_phase) 

# External electric field (step function)
#E_amplitude_max<- 0.01* 10^8 # V/m
#E_vector[1:(div_num_t+1)] <-  rep(0,div_num_t+1)
#E_vector[101:(div_num_t+1)] <-  rep(E_amplitude_max, div_num_t-99)

# k normalized by inverse lattice constant 1/a 
div_num_ka <- 3000 # m
delta_k <- 2*pi / lattice_constant_a_of_GaSe / div_num_ka 
ka_vector_redundant <- seq( -pi, pi, length= (div_num_ka+1))
ka_vector <- ka_vector_redundant[1:div_num_ka]

#band structure of GaSe (two band model assuming the Tight binding model (1D))
band_width_valence <- 1  # eV 
band_width_conduction <- 1.5  # eV 
band_gap <- 2 # eV 
band_energy_valence <- band_width_valence* (cos( ka_vector ) -1) /2 # eV 
band_energy_conduction <- band_width_conduction* (1 - cos( ka_vector ) )/2  +  band_gap # eV 
png("band_structure.png", width = 700, height = 500) 
plot(ka_vector, band_energy_valence, xlab="ka (rad)", ylab="energy (eV)", ylim=c(-2,3),col="red", type="l")  
par(new=T)
plot(ka_vector, band_energy_conduction, xlab="", ylab="", ylim=c(-2,3),col="blue", type="l")  
dev.off()            
band_energy_valence <- band_energy_valence *1.60218*10^-19  # conversion from eV to J
band_energy_conduction <- band_energy_conduction *1.60218*10^-19  # conversion from eV to J     
energy_diff <- band_energy_conduction -band_energy_valence
J_hole_k <- elementary_charge / Plank_constant * band_width_valence *lattice_constant_a_of_GaSe *(-sin(ka_vector)) /2
J_electron_k <- elementary_charge /Plank_constant * band_width_conduction *lattice_constant_a_of_GaSe *sin(ka_vector)/2
png("J_k.png", width = 700, height = 500)  # 描画デバイスを開く 
plot(ka_vector, J_hole_k, xlab="ka (rad)", ylab="j",col="red", type="l")  
par(new=T)
plot(ka_vector, J_electron_k, xlab="", ylab="", col="blue", type="l")  
dev.off()            

#dipole transition matrix
Rabi_energy_max <- 2.0 *1.60218*10^-19 # converstion from eV to J ; 1.2eV
dipole_transition <- Rabi_energy_max / E_amplitude_max

#relaxation
relaxation_constant <- 20 *10^(-15) # 7fs 

#dephasing 
dephasing_constant <- 5 *10^(-15) # 1.1fs 

number_hole_kt <- array(0, dim=c(div_num_ka, (div_num_t+1)))
number_electron_kt <- array(0, dim=c(div_num_ka, (div_num_t+1)))
polarization_kt  <- array(0, dim=c(div_num_ka, (div_num_t+1)))

number_hole <- rep(0, div_num_ka)
number_electron <- rep(0, div_num_ka) 
polarization <- rep(0, div_num_ka)
dn <- rep(0, div_num_ka)
dp <- rep(0, div_num_ka)

J_hole <- rep(0, (div_num_t+1))
J_electron <- rep(0, (div_num_t+1))
P <- rep(0, (div_num_t+1))
P_diff <- rep(0, (div_num_t+1))

#solving SBE by the 4th Runge-Kutta method
fn <- function(E, p, n){
	dn[1] <- n[div_num_ka]
	dn[2:div_num_ka] <- n[1:(div_num_ka-1)]
	nn <- ( 2*Im(E*dipole_transition*p) / Plank_constant
	+ elementary_charge* E* (n-dn) / delta_k / Plank_constant
	- n / relaxation_constant ) * delta_t
	return (nn)
}

fp <- function(E, p, nh, ne){
	dp[1] <- p[div_num_ka]
	dp[2:div_num_ka] <- p[1:(div_num_ka-1)]	
	pp <- ( -1i *energy_diff *p / Plank_constant
	+ elementary_charge* E* (p-dp) / delta_k  / Plank_constant
	+ 1i* E* dipole_transition * (1- nh- ne) / Plank_constant
	- p / dephasing_constant )* delta_t
	return (pp)
}

for( j in 1:(div_num_t) ){
	E <- E_vector2[2*j-1]
	E_0 <- E_vector2[2*j]
	E_1 <- E_vector2[2*j]	
	E_2 <- E_vector2[2*j+1]	

	ne_0 <- fn(E, polarization, number_electron)
	nh_0 <- fn(E, polarization, number_hole )
	p_0 <- fp(E, polarization, number_hole, number_electron)
	
	ne_1 <- fn(E_0, polarization +p_0/2, number_electron +ne_0/2 )
	nh_1 <- fn(E_0, polarization +p_0/2, number_hole + nh_0/2 )
	p_1 <- fp(E_0, polarization+p_0/2, number_hole + nh_0/2, number_electron +ne_0/2)

	ne_2 <- fn(E_1, polarization +p_1/2, number_electron +ne_1/2 )
	nh_2 <- fn(E_1, polarization +p_1/2, number_hole + nh_1/2 )
	p_2 <- fp(E_1, polarization+p_1/2, number_hole + nh_1/2, number_electron +ne_1/2)
	
	ne_3 <- fn(E_2, polarization +p_2, number_electron +ne_2 )
	nh_3 <- fn(E_2, polarization +p_2, number_hole + nh_2 )
	p_3 <- fp(E_2, polarization+p_2, number_hole + nh_2, number_electron +ne_2)	

	number_hole <- number_electron + (nh_0 +2*nh_1 +2*nh_2 +nh_3)/6
	number_electron <- number_hole + (ne_0 +2*ne_1 +2*ne_2 +ne_3)/6
	polarization <- polarization + (p_0 +2*p_1 +2*p_2 +p_3)/6

	number_hole_kt[,j] <- number_electron
	number_electron_kt[,j] <- number_hole
	polarization_kt[,j] <- polarization
}


	for(i in 1:(div_num_ka)){
		J_hole <-  J_hole + J_hole_k[i] * number_hole_kt[i,]
		J_electron <- J_electron + J_electron_k[i] * number_electron_kt[i,]
		P <- P + polarization_kt[i,] * dipole_transition	
	}

P_diff[2:(div_num_t+1)] <- P[1:div_num_t]
dP_dt <- (P-P_diff)/delta_t
J <- J_hole + J_electron
E_HHG <- J +dP_dt

warnings()

if(1){
library(animation)
saveGIF({
  ani.options(loop = TRUE)
  for (i in seq(1,(div_num_t+1), by=150)) {
    .main = paste(t_vector[i]*10^15,"(fs)", seq="")
    plot(ka_vector, number_electron_kt[,i], xlim=c(-pi,pi), ylim=c(-1, 1), xlab="BZ", ylab= "numbers", main=.main,col="green", type="l")
    par(new=T)
    plot(ka_vector, number_hole_kt[,i], xlim=c(-pi,pi), ylim=c(-1, 1), xlab="", ylab= "", main=.main, col="blue", type="l")
    par(new=T)    
    plot(ka_vector, Re(polarization_kt[,i]), xlim=c(-pi,pi), ylim=c(-1, 1), xlab="BZ", ylab= "polarization", main=.main,col="red" , type="l")
  }
}, movie.name="number_polarization.gif", interval=0.3)
}

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

Erange <- c( min(Re(E_HHG)), max(Re(E_HHG)) ) 
png("E.png", width = 700, height = 500)  
plot(t_vector*10^15, E_vector, xlab="t (fs)", ylab="Electric field (V/m)",col="black", type="l")  
par(new=T)
plot(t_vector*10^15, Re(E_HHG), xlab="", ylab= "", ylim=Erange, col="red", type="l")  
par(new=T)
plot(t_vector*10^15, Re(dP_dt), xlab="", ylab= "", ylim=Erange, col="green", type="l")  
par(new=T)
plot(t_vector*10^15, Re(J), xlab="", ylab="", ylim=Erange, col="blue", type="l")  
dev.off()       


#write.table(number_electron_kt, "./number_electron_kt.txt")


#png("number.png", width = 700, height = 500)  
#plot(number_electron_kt[], xlab="", ylab="")  
#dev.off() 

E_spectrum <- fft(E_vector)
E_HHG_spectrum <- fft(E_HHG)


frange <- c(0, 20)
png("E_spectrum.png", width = 700, height = 500)  
plot(frequency/carrier_frequency, abs(E_spectrum)^2, xlim=frange, log="y", xlab="order", ylab="power spectrum",col="blue", type="l")  
par(new=T)
plot(frequency/carrier_frequency, abs(E_HHG_spectrum)^2, xlim=frange, log="y", xlab="", ylab="", col="red", type="l")  
abline(v=0:20, col='black', lty="dotted")
dev.off()     


