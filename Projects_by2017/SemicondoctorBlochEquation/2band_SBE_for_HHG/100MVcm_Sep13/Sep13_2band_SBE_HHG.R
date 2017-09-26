#I use SI unit system for calculation
# first written in 6 Sep. 2017
#modified in 11 Sep. 2017

options(scipen=0)

#Physical constants 
Plank_constant <- 6.62607004 * 10^-34 #m^2 kg / s
Boltzmann_constant <-  1.38064852 *10^-23 # m^2 kg s^-2 K^-1
speed_of_light  <- 299792458 # m/s
elementary_charge <- 1.60217662 * 10^-19 # Coulombs

#Parameters of experiment / material properties
Temperature <- 300 #K
Fermi_energy <- # Chemical potential
lattice_constant_a_of_GaSe <- 3.755 * 2/sqrt(3) *10^-10# m (3.755 ångström)

# Time t
range_t <- 200 *10^(-15) # fs
div_num_t <- 1000000 # n
delta_t <-  range_t/ div_num_t
t_vector <- seq( -range_t/2, range_t/2, length= (div_num_t+1))
t_vector2 <- seq( -range_t/2, range_t/2, length= (2*div_num_t+1))
frequency <- seq(0, 1/delta_t, length= (div_num_t+1))

# External electric field (assuming the gaussian envelope)
E_amplitude_max <- 100* 10^8 # V/m
envelope_width <- 50* 10^(-15) # 1/e^2 width; 50 fs 
carrier_frequency <- 33 *10^12 # 33 THz
carrier_envelope_phase <- 0 # [rad]
E_vector <-  E_amplitude_max * exp(-t_vector^2/envelope_width^2) *cos(2*pi* carrier_frequency* t_vector+carrier_envelope_phase) 
E_vector2 <-  E_amplitude_max * exp(-t_vector2^2/envelope_width^2) *cos(2*pi* carrier_frequency* t_vector2+carrier_envelope_phase) 
E_spectrum <- fft(E_vector)

# External electric field (step function)
#E_vector[1:(div_num_t+1)] <-  rep(0,div_num_t+1)
#E_vector[101:(div_num_t+1)] <-  rep(E_amplitude_max, div_num_t-99)
#E_vector2[1:(2*div_num_t+1)] <-  rep(0,2*div_num_t+1)
#E_vector2[201:(2*div_num_t+1)] <-  rep(E_amplitude_max, 2*div_num_t-199)


# k normalized by inverse lattice constant 1/a 
div_num_ka <- 150 # m
delta_k <- 2*pi / lattice_constant_a_of_GaSe / div_num_ka 
ka_vector_redundant <- seq( -pi, pi, length= (div_num_ka+1))
ka_vector <- ka_vector_redundant[1:div_num_ka]

#band structure of GaSe (two band model assuming the Tight binding model (1D))
band_width_valence <- 1  # eV 
band_width_conduction <- 1.5  # eV 
band_gap <- 2 # eV 
band_energy_valence <- band_width_valence* (cos( ka_vector ) -1) /2 # eV 
band_energy_conduction <- band_width_conduction* (1 - cos( ka_vector ) )/2  +  band_gap # eV 
erange <- c(-2,5)
png("band_structure.png", width = 700, height = 500) 
plot(ka_vector, band_energy_valence, xlab="ka (rad)", ylab="energy (eV)", ylim=erange,col="red", type="l")  
par(new=T)
plot(ka_vector, band_energy_conduction, xlab="", ylab="", ylim=erange,col="blue", type="l")  
dev.off()            
band_energy_valence <- band_energy_valence *1.60218*10^-19  # conversion from eV to J
band_energy_conduction <- band_energy_conduction *1.60218*10^-19  # conversion from eV to J     
energy_diff <- band_energy_conduction -band_energy_valence
J_hole_k <- elementary_charge / Plank_constant * band_width_valence *lattice_constant_a_of_GaSe *(-sin(ka_vector)) /2
J_electron_k <- -elementary_charge /Plank_constant * band_width_conduction *lattice_constant_a_of_GaSe *sin(ka_vector)/2
png("J_k.png", width = 700, height = 500)  # 描画デバイスを開く 
plot(ka_vector, J_hole_k, xlab="ka (rad)", ylab="j",col="red", type="l")  
par(new=T)
plot(ka_vector, J_electron_k, xlab="", ylab="", col="blue", type="l")  
dev.off()            


#dipole transition matrix
Rabi_energy_max <- 0.12 *1.60218*10^-19 #1.2eV; converstion from eV to J
dipole_transition_max <- Rabi_energy_max / E_amplitude_max 
dipole_transition <- dipole_transition_max * min(energy_diff) / energy_diff 

#relaxation
relaxation_constant <-  50 *10^(-15) # 7fs 
#dephasing 
dephasing_constant <- 5 *10^(-15) # 1.1fs 

#parameters of calculation
Courant_number <- delta_t * E_amplitude_max *  elementary_charge / Plank_constant / delta_k 
print("Courant number: " )
print( Courant_number )

w_2 <- rep(0, div_num_ka)
w_1 <- rep(0, div_num_ka)
w1 <- rep(0, div_num_ka)
w2 <- rep(0, div_num_ka)
ww <- rep(0, div_num_ka)

grad_k0 <- function(E, w0){
	w1[div_num_ka] <- w0[1]
	w1[1:(div_num_ka-1)] <- w0[2:div_num_ka]
	return ((w0-w1) / delta_k)
}

grad_k1 <- function(E, w0){
	w_1[1] <- w0[div_num_ka]
	w_1[2:div_num_ka] <- w0[1:(div_num_ka-1)]
	w1[div_num_ka] <- w0[1]
	w1[1:(div_num_ka-1)] <- w0[2:div_num_ka]
	return ((w1-w_1)/ 2 / delta_k)
}

grad_k2 <- function(E, w0){
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

fn <- function(E, p, n){
	nn <- ( 2*Im(E*dipole_transition*p) / Plank_constant
	+ elementary_charge* E* grad_k1(E, n) / Plank_constant
	- n / relaxation_constant ) * delta_t
	return (nn)
}

fp <- function(E, p, nh, ne){
	pp <- ( -1i *energy_diff *p / Plank_constant
	+ elementary_charge* E* grad_k1(E, p) / Plank_constant
	+ 1i* E* dipole_transition * (1- nh- ne) / Plank_constant
	- p / dephasing_constant )* delta_t
	return (pp)
}

number_hole_kt <- array(0, dim=c(div_num_ka, (div_num_t+1)))
number_electron_kt <- array(0, dim=c(div_num_ka, (div_num_t+1)))
polarization_kt  <- array(0, dim=c(div_num_ka, (div_num_t+1)))

number_hole <- rep(0, div_num_ka)
number_electron <- rep(0, div_num_ka) 
polarization <- rep(0, div_num_ka)

J_hole <- rep(0, (div_num_t+1))
J_electron <- rep(0, (div_num_t+1))
P <- rep(0, (div_num_t+1))
P_diff <- rep(0, (div_num_t+1))

#solving SBE by the 4th Runge-Kutta method

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
		P <- P + polarization_kt[i,] * dipole_transition	[i]
	}



P_diff[2:(div_num_t+1)] <- P[1:div_num_t]
dP_dt <- (P-P_diff)/delta_t
J <- J_hole + J_electron
E_HHG <- J +dP_dt

E_HHG_spectrum <- fft(E_HHG)
J_spectrum <-  fft(J)
dP_dt_spectrum <- fft(dP_dt)

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


