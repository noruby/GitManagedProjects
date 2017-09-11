#I use SI unit system for calculation
# first written in 6 Sep. 2017
#modified in 11 Sep. 2017

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
div_num_t <- 2000 # n
delta_t <-  range_t/ div_num_t
delta_t_half <- delta_t
t_vector <- seq( -range_t/2, range_t/2, length= (div_num_t+1))

# External electric field (assuming the gaussian envelope)
E_amplitude_max <- 2* 10^8 # V/m
envelope_width <- 70* 10^(-15) # 1/e^2 width; 50 fs 
carrier_frequency <- 33 *10^12 # 33 THz
carrier_envelope_phase <- 0 # [rad]
E_vector <-  E_amplitude_max * exp(-t_vector^2/envelope_width^2) *cos(2*pi* carrier_frequency* t_vector+carrier_envelope_phase) 

# External electric field (step function)
#E_amplitude_max<- 0.01* 10^8 # V/m
#E_vector[1:(div_num_t+1)] <-  rep(0,div_num_t+1)
#E_vector[101:(div_num_t+1)] <-  rep(E_amplitude_max, div_num_t-99)

# k normalized by inverse lattice constant 1/a 
div_num_ka <- 1000 # m
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
plot(ka_vector, band_energy_valence, xlab="ka (rad)", ylab="energy (eV)", ylim=c(-2,3),col="red")  
par(new=T)
plot(ka_vector, band_energy_conduction, xlab="", ylab="", ylim=c(-2,3),col="blue")  
dev.off()            
band_energy_valence <- band_energy_valence *1.60218*10^-19  # conversion from eV to J
band_energy_conduction <- band_energy_conduction *1.60218*10^-19  # conversion from eV to J     
J_hole_k <- elementary_charge / Plank_constant * band_width_valence *lattice_constant_a_of_GaSe *(-sin(ka_vector)) /2
J_electron_k <- elementary_charge /Plank_constant * band_width_conduction *lattice_constant_a_of_GaSe *sin(ka_vector)/2
png("J_k.png", width = 700, height = 500)  # 描画デバイスを開く 
plot(ka_vector, J_hole_k, xlab="ka (rad)", ylab="j",col="red")  
par(new=T)
plot(ka_vector, J_electron_k, xlab="", ylab="", col="blue")  
dev.off()            

#dipole transition matrix
Rabi_energy_max <- 0.012 *1.60218*10^-19 # converstion from eV to J 
dipole_transition <- Rabi_energy_max / E_amplitude_max

#relaxation
relaxation_constant <- 7 *10^(-15) # 7fs 

#dephasing 
dephasing_constant <- 1.1 *10^(-15) # 1.1fs 

number_hole_kt <- array(0, dim=c(div_num_ka, (div_num_t+1)))
number_electron_kt <- array(0, dim=c(div_num_ka, (div_num_t+1)))
polarization_kt  <- array(0, dim=c(div_num_ka, (div_num_t+1)))

number_hole <- rep(0, div_num_ka)
number_hole_diff <- rep(0, div_num_ka)
number_electron <- rep(0, div_num_ka) 
number_electron_diff <- rep(0, div_num_ka) 
polarization <- rep(0, div_num_ka)
polarization_diff <- rep(0, div_num_ka)

J_hole <- rep(0, (div_num_t+1))
J_electron <- rep(0, (div_num_t+1))
P <- rep(0, (div_num_t+1))
P_diff <- rep(0, (div_num_t+1))

#solving SBE by the 4th Runge-Kutta method
fn <- function(E, p, n){
	+2*Im(E*dipole_transition*p) / Plank_constant
	+elementary_charge* E* (number_hole-number_hole_diff) / delta_k / Plank_constant
	- number_hole / relaxation_constant )
	* delta_t
}
fp <- function(E, p, nh, ne){
	-1i *energy_diff
	 elementary_charge* E* (polarization-polarization_diff) / delta_k ) *polarization
	1i* E* dipole_transition * (1- number_hole- number_electron) )/Plank_constant
	- polarization / dephasing_constant )
	* delta_t
}

for( j in 1:(div_num_t+1) ){
	E <- E_vector[j]
	
	number_hole_diff[1] <- number_hole[div_num_ka]
	number_hole_diff[2:div_num_ka] <- number_hole[1:(div_num_ka-1)]
	number_hole_new <- number_hole + ( ( 2* Im( E* dipole_transition * polarization ) + elementary_charge* E* (number_hole-number_hole_diff) / delta_k ) / Plank_constant  - number_hole / relaxation_constant ) * delta_t
	number_hole <- number_hole_new
	number_hole_kt[,j] <- number_hole_new

	number_electron_diff[1] <- number_electron[div_num_ka]
	number_electron_diff[2:div_num_ka] <- number_electron[1:(div_num_ka-1)]
	number_electron_new <- number_electron + ( ( 2* Im( E* dipole_transition * polarization ) + elementary_charge* E* (number_electron-number_electron_diff) / delta_k ) /Plank_constant  - number_electron/ relaxation_constant ) * delta_t
	number_electron <- number_electron_new
	number_electron_kt[,j] <- number_electron_new
		
	energy_diff <- band_energy_conduction - band_energy_valence
	polarization_diff[1] <- polarization[div_num_ka]
	polarization_diff[2:div_num_ka] <- polarization[1:(div_num_ka-1)]		
	polarization_new  <- polarization + ( ( ( -1i *energy_diff + elementary_charge* E* (polarization-polarization_diff) / delta_k ) *polarization  +1i* E* dipole_transition * (1- number_hole- number_electron) )/Plank_constant - polarization / dephasing_constant ) * delta_t
	polarization <- polarization_new
	polarization_kt[,j] <- polarization_new
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
  for (i in seq(1,(div_num_t+1), by=50)) {
    .main = paste(t_vector[i]*10^15,"(fs)", seq="")
    plot(ka_vector, number_electron_kt[,i], xlim=c(-pi,pi), ylim=c(0, 1), xlab="BZ", ylab= "numbers", main=.main,col="green")
    par(new=T)
    plot(ka_vector, number_hole_kt[,i], xlim=c(-pi,pi), ylim=c(0, 1), xlab="", ylab= "", main=.main, col="blue")
    par(new=T)    
    plot(ka_vector, Re(polarization_kt[,i]), xlim=c(-pi,pi), ylim=c(0, 1), xlab="BZ", ylab= "polarization", main=.main,col="red")
  }
}, movie.name="number_polarization.gif", interval=0.3)
}

Jrange <- 3*c( min(Re(J)), max(Re(J)))
png("J.png", width = 700, height = 500)  
plot(t_vector*10^15, Re(J), xlab="t (fs)", ylab= "J", ylim=Jrange, col="red")  
par(new=T)
plot(t_vector*10^15, Re(J_electron), xlab="", ylab= "", ylim=Jrange, col="green")  
par(new=T)
plot(t_vector*10^15, Re(J_hole), xlab="", ylab="", ylim=Jrange, col="blue")  
dev.off()    

Prange <- c( min(abs(P)), max(abs(P)))
png("P.png", width = 700, height = 500)  
plot(t_vector*10^15, Re(P), xlab="t (fs)", ylab= "P", ylim=Prange, col="red")  
par(new=T)
plot(t_vector*10^15, Im(P), xlab="", ylab="", ylim=Prange, col="blue")  
dev.off()       

Erange <- c( min(Re(E_HHG)), max(Re(E_HHG)) ) 
png("E.png", width = 700, height = 500)  
plot(t_vector*10^15, E_vector, xlab="t (fs)", ylab="Electric field (V/m)",col="black")  
par(new=T)
plot(t_vector*10^15, Re(E_HHG), xlab="", ylab= "", ylim=Erange, col="red")  
par(new=T)
plot(t_vector*10^15, Re(dP_dt), xlab="", ylab= "", ylim=Erange, col="green")  
par(new=T)
plot(t_vector*10^15, Re(J), xlab="", ylab="", ylim=Erange, col="blue")  
dev.off()       


#write.table(number_electron_kt, "./number_electron_kt.txt")


#png("number.png", width = 700, height = 500)  
#plot(number_electron_kt[], xlab="", ylab="")  
#dev.off() 

#E_HHG_w <- fft(E_HHG_t)

