#I use SI unit system for calculation
# first written in 6 Sep. 2017
# modified in 8 Jan. 2018

#assumptions
#1 dimention
#the dipole transition amplitude is independent on the crystal momentum k
#
#
#

#Physical constants 
Planck_constant_bar <- 6.62607004 * 10^-34 #m^2 kg / s
Boltzmann_constant <-  1.38064852 *10^-23 # m^2 kg s^-2 K^-1
speed_of_light  <- 299792458 # m/s
elementary_charge <- 1.60217662 * 10^-19 # Coulombs
mass_of_electron <- 9.10938356 * 10^-31 #kg
effective_mass_valence <- mass_of_electron
effective_mass_conduction <- mass_of_electron

#Parameters of experiment / material properties
lattice_constant <- 3.755 * 2/sqrt(3) *10^-10# m (3.755 ångström)

# Time t
range_t <- 300 *10^(-15) # fs
div_num_t <- 500 # n
delta_t <-  range_t/ div_num_t
t_vector <- seq( -range_t/2, range_t/2, length= (div_num_t+1))
t_vector2 <- seq( -range_t/2, range_t/2, length= (2*div_num_t+1))
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

# fixed k normalized by inverse lattice constant 1/a in length gauge
div_num_ka <- 50 # m
delta_k <- 2*pi / lattice_constant / div_num_ka 
ka_vector0_redundant <- seq( -pi, pi, length= (div_num_ka+1)) #BZの両端で冗長
ka_vector0 <- ka_vector0_redundant[1:div_num_ka]
k_vector0 <- ka_vector0 / lattice_constant

#band structure of GaSe (two band model assuming the Tight binding model (1D))
band_gap <- 2 *1.60218*10^-19 # J (2eV) 
band_width_valence <- 1 *1.60218*10^-19 # J (1eV)
band_width_conduction <- -1.5 *1.60218*10^-19 # J (1.5eV) 
mean_energy_valence <- 0 
mean_energy_conduction <- mean_energy_valence +band_gap +(band_width_valence-band_width_conduction)/2
energy_valence <- function(ka){
	band_width_valence/2*cos(ka) +mean_energy_valence # J
}
energy_conduction <- function(ka){
	band_width_conduction/2*cos(ka) +mean_energy_conduction # J 
}
energy_valence_vector0 <- energy_valence(ka_vector0)
energy_conduction_vector0 <- energy_conduction(ka_vector0)
energy_difference_vector0 <- energy_conduction_vector0 - energy_valence_vector0
w_1 <- rep(0, div_num_ka)
w1 <- rep(0, div_num_ka)
grad_k_midpoint <- function(w0){
	w_1[1] <- w0[div_num_ka]
	w_1[2:div_num_ka] <- w0[1:(div_num_ka-1)]
	w1[div_num_ka] <- w0[1]
	w1[1:(div_num_ka-1)] <- w0[2:div_num_ka]
	return ((w1-w_1)/ 2 / delta_k)
}

#microscopic current
grad_energy_valence_vector0 <- grad_k_midpoint(energy_valence_vector0)
grad_energy_conduction_vector0 <- grad_k_midpoint(energy_conduction_vector0)
if(0){
print("Exporting debug_grad_energy.png ...")
png("./debug_grad_energy.png", width = 700, height = 500)  
krange <- c(-pi, pi)
eVrange <- 0.03*eVrange*(1.60218*10^-19)/delta_k
plot(ka_vector0, grad_energy_valence_vector0, xlim=krange, ylim=eVrange, xlab="BZ", ylab= "energy (eV)", col="red", type="l")  
par(new=T)
plot(ka_vector0,grad_energy_conduction_vector0, xlim=krange, ylim=eVrange, xlab="", ylab="", col="blue", type="l")  
dev.off() 
}

#dipole transition matrix
Rabi_energy_max <- 1.2*1.60218*10^-19 #1.2eV; converstion from eV to J
dipole_transition_amplitude <- Rabi_energy_max/(50*10^8)
d_vv <- rep(0, div_num_ka)
d_cc <- rep(0, div_num_ka)
d <- d_vc <- rep(dipole_transition_amplitude, div_num_ka)

#parameters of calculation
Courant_number <- delta_t * E_amplitude_max *  elementary_charge / Planck_constant_bar / delta_k 
print("Courant number (length gauge): " )
print( Courant_number )

#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------
#calculate ---------------------------------------------------------------------------------------------------------------
calc_bloch <- "no"
#dephasing
#dephasing_constant <- 500 *10^(-15) # 1.1fs 

if(calc_bloch=="yes"){
diff_fc <- function(E, p, fc){
	dfc <- E/Planck_constant_bar*(-2*Im(p*Conj(d)) +elementary_charge*grad_k_midpoint(fc))*delta_t
	return (dfc)
}
diff_fv <- function(E, p, fv){
	dfv <- E/Planck_constant_bar*(2*Im(p*Conj(d)) +elementary_charge*grad_k_midpoint(fv))*delta_t
	return (dfv)
}
diff_p <- function(E, p, fc, fv){
	dp <- (1i*energy_difference_vector0*p +elementary_charge*E*grad_k_midpoint(p)-1i*E*d*(fv-fc))*delta_t/Planck_constant_bar
	#dp <- dp -p/dephasing_constant*delta_t
	return (dp)
}

fv_kt <- array(1, dim=c(div_num_ka, (div_num_t+1))) #number_valence_length
fc_kt <- array(0, dim=c(div_num_ka, (div_num_t+1))) #number_conduction_length
p_kt  <- array(0, dim=c(div_num_ka, (div_num_t+1))) #micro_polarization_length

fv <- rep(1, div_num_ka)
fc <- rep(0, div_num_ka) 
p <- rep(0, div_num_ka)

print("starting calculation with length gauge")
for( j in 1:(div_num_t) ){
	E <- E_vector2[2*j-1]
	E_0 <- E_vector2[2*j]
	E_1 <- E_vector2[2*j]	
	E_2 <- E_vector2[2*j+1]	

	dfc_0 <- diff_fc(E, p, fc)
	dfv_0 <- diff_fv(E, p, fv)
	dp_0 <- diff_p(E, p, fc, fv)

	dfc_1 <- diff_fc(E_0, p+dp_0/2, fc+dfc_0/2)
	dfv_1 <- diff_fv(E_0, p+dp_0/2, fv+dfv_0/2)
	dp_1 <- diff_p(E_0, p+dp_0/2, fc+dfc_0/2, fv+dfv_0/2)

	dfc_2 <- diff_fc(E_1, p+dp_1/2, fc+dfc_1/2)
	dfv_2 <- diff_fv(E_1, p+dp_1/2, fv+dfv_1/2)
	dp_2 <- diff_p(E_1, p+dp_1/2, fc+dfc_1/2, fv+dfv_1/2)

	dfc_3 <- diff_fc(E_2, p+dp_2, fc+dfc_2)
	dfv_3 <- diff_fv(E_2, p+dp_2, fv+dfv_2)
	dp_3 <- diff_p(E_2, p+dp_2, fc+dfc_2, fv+dfv_2)

	fc <- fc + (dfc_0 +2*dfc_1 +2*dfc_2 +dfc_3)/6
	fv <- fv + (dfv_0 +2*dfv_1 +2*dfv_2 +dfv_3)/6
	p <- p + (dp_0 +2*dp_1 +2*dp_2 +dp_3)/6

	fc_kt[,j+1] <- fc
	fv_kt[,j+1] <- fv
	p_kt[,j+1] <- p
}

if(0){
Current_valence_length <- rep(0, (div_num_t+1))
Current_conduction_length <- rep(0, (div_num_t+1))
Polarization_current_length <- rep(0, (div_num_t+1))
Energy_free_length <- rep(0, (div_num_t+1))
Energy_interaction_length <- rep(0, (div_num_t+1))
for(i in 1:(div_num_ka)){
	Current_valence_length <-  Current_valence_length - elementary_charge/Planck_constant_bar* grad_energy_valence_vector0[i] * nvl_kt[i,] *delta_k
	Current_conduction_length <- Current_conduction_length - elementary_charge/Planck_constant_bar* grad_energy_conduction_vector0[i] * ncl_kt[i,] *delta_k 
	Polarization_current_length <- Polarization_current_length + elementary_charge/Planck_constant_bar* energy_difference_vector0[i]* 2* Im(mpl_kt[i,] * Conj(d_vc[i])) *delta_k
	Energy_free_length <- Energy_free_length + (energy_valence_vector0[i] * nvl_kt[i,] + energy_conduction_vector0[i] * ncl_kt[i,] ) *delta_k
	Energy_interaction_length <- Energy_interaction_length -elementary_charge* E_vector *(d_vv[i] * nvl_kt[i,] +d_cc[i] * ncl_kt[i,] +2*Re(Conj(d_vc[i])*mpl_kt[i,])) *delta_k
}

Current_length <- Current_valence_length + Current_conduction_length
E_HHG_length <- Current_length + Polarization_current_length
Energy_length <- Energy_free_length + Energy_interaction_length

E_spectrum <- fft(E_vector)
E_HHG_length_spectrum <- fft(E_HHG_length)
Current_length_spectrum <-  fft(Current_length)
Polarization_current_length_spectrum <- fft(Polarization_current_length)

print("warnings: ")
warnings()
}

#---------------------------------------------------------------------------------------------------------------
if(0){
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
par(new=T)
plot(t_vector*10^15, Re(Current_valence_length), xlab="", ylab= "", xlim=trange, ylim=E_HHGrange, axes=FALSE, yaxt="n", col="yellow", type="l")  
par(new=T)
plot(t_vector*10^15, Re(Current_valence_length), xlab="", ylab= "", xlim=trange, ylim=E_HHGrange, axes=FALSE, yaxt="n", col="pink", type="l")  
dev.off()    
}

if(0){
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
}

#
if(1){
dir.exists <- function(d) {
    de <- file.info(d)$isdir
    ifelse(is.na(de), FALSE, de)
}
print("Exporting fp.gif ...")
numt <- round( 2*10^-15 /delta_t )
nplim=c(-1.2, 1.2)
library(animation)
saveGIF({
  ani.options(loop = TRUE)
  for (i in seq(1,(div_num_t+1), by=numt)) {
    .main = paste(t_vector[i]*10^15,"(fs)", seq="")
    plot(ka_vector0, fc_kt[,i], xlim=c(-pi,pi), ylim=nplim, xlab="BZ", ylab= "numbers", main=.main,col="green", type="l")
    par(new=T)
    plot(ka_vector0, fv_kt[,i], xlim=c(-pi,pi), ylim=nplim, xlab="", ylab= "", main=.main, col="blue", type="l")
    par(new=T)    
    plot(ka_vector0, Re(p_kt[,i]), xlim=c(-pi,pi), ylim=nplim, xlab="BZ", ylab= "polarization", main=.main,col="red" , type="l")
  }
}, movie.name="fp.gif", interval=0.2)
}

}

#==============================================
basis <- "Wannier"
if(basis=="Wannier"){
CV<-C<- matrix(0, nrow=div_num_ka,ncol=div_num_ka)
V <- diag(div_num_ka)
dCV<-dV<-dC<- matrix(0, nrow=div_num_ka,ncol=div_num_ka)
CV_kt<-V_kt<-C_kt<- array(0, dim=c(div_num_ka,div_num_ka,(div_num_t+1)))
pW_kt<-fvW_kt<-fcW_kt<- array(0, dim=c(div_num_ka,(div_num_t+1)))

Blavais_lattice <- seq(-lattice_constant*div_num_ka/2, by=lattice_constant, length=div_num_ka)
RsubR <- rep(1,length=div_num_ka)%*%t(Blavais_lattice)-Blavais_lattice%*%t(rep(1,length=div_num_ka)) 
eRsubR <- elementary_charge* RsubR
expF<- array(0, dim=c(div_num_ka,div_num_ka,div_num_ka))
for(i in 1:div_num_ka){
  expF[,,i] <- exp(1i*k_vector0[i]*RsubR)
}
EV <- EC <- matrix(0,ncol=div_num_ka, nrow=div_num_ka)
for(i in 1:div_num_ka){
for(j in 1:div_num_ka){
　　if(i==j){
    EC[i,j] <- mean_energy_conduction
    EV[i,j] <- mean_energy_valence
  }
  if(i==(j+1)||(i+1)==j){
    EC[i,j] <- band_width_conduction/2
    EV[i,j] <- band_width_valence/2
  }
}
}
EC[1,div_num_ka] <- EC[div_num_ka,1]<- band_width_conduction/2
EV[1,div_num_ka] <- EV[div_num_ka,1]<- band_width_valence/2

diff_C <- function(E,C,CV){
 dC <- C%*%EC-EC%*%C -E*2i*Im(Conj(d_vc)*CV) +E*eRsubR*C 
		  return (dC/1i/Planck_constant_bar*delta_t)
}
diff_V <- function(E,V,CV){
 dC <- V%*%EV-EV%*%V +E*2i*Im(Conj(d_vc)*CV) +E*eRsubR*V 
		  return (dV/1i/Planck_constant_bar*delta_t)
}
diff_CV <- function(E,C,V,CV){
 dCV <- CV%*%EV-EC%*%CV -E*Conj(d_vc)*(C-V) +E*eRsubR*CV 
		  return (dV/1i/Planck_constant_bar*delta_t)
}

print("starting calculation with the Wannier basis")
nn <- 0.1
for( j in 1:(div_num_t) ){
	E <- E_vector2[2*j-1]
	E_0 <- E_vector2[2*j]
	E_1 <- E_vector2[2*j]	
	E_2 <- E_vector2[2*j+1]	

   dC_0 <- diff_C(E,C,CV)
   dV_0 <- diff_C(E,V,CV)
   dCV_0 <- diff_CV(E,C,V,CV)

   dC_1 <- diff_C(E_0,C+dC_0/2,CV+dCV_0/2)
   dV_1 <- diff_C(E_0,V+dV_0/2,CV+dCV_0/2)
   dCV_1 <- diff_CV(E_0,C+dC_0/2,V+dV_0/2,CV+dCV_0/2)

   dC_2 <- diff_C(E_1,C+dC_1/2,CV+dCV_1/2)
   dV_2 <- diff_C(E_1,V+dV_1/2,CV+dCV_1/2)
   dCV_2 <- diff_CV(E_1,C+dC_1/2,V+dV_1/2,CV+dCV_1/2)

   dC_3 <- diff_C(E_2,C+dC_2,CV+dCV_2)
   dV_3 <- diff_C(E_2,V+dV_2,CV+dCV_2)
   dCV_3 <- diff_CV(E_2,C+dC_2,V+dV_2,CV+dCV_2)

   C <- C + (dC_0 +2*dC_1+ 2*dC_2+ dC_3)/6
   V <- V + (dV_0 +2*dV_1+ 2*dV_2+ dV_3)/6
   CV <- CV + (dCV_0 +2*dCV_1+ 2*dCV_2+ dCV_3)/6
   
   for(i in 1:div_num_ka){
     fcW_kt[i,j] <- sum(C*expF[,,i])/div_num_ka
     fvW_kt[i,j] <- sum(V*expF[,,i])/div_num_ka
     pW_kt[i,j] <- sum(CV*expF[,,i])/div_num_ka
   }

   V_kt[,,j] <- V

   if((j/div_num_t)>nn){
     print(nn*100)
     nn <- nn + 0.1 
   }
}


if(0){
Current_valence_length <- rep(0, (div_num_t+1))
Current_conduction_length <- rep(0, (div_num_t+1))
Polarization_current_length <- rep(0, (div_num_t+1))
Energy_free_length <- rep(0, (div_num_t+1))
Energy_interaction_length <- rep(0, (div_num_t+1))
for(i in 1:(div_num_ka)){
	Current_valence_length <-  Current_valence_length - elementary_charge/Planck_constant_bar* grad_energy_valence_vector0[i] * nvl_kt[i,] *delta_k
	Current_conduction_length <- Current_conduction_length - elementary_charge/Planck_constant_bar* grad_energy_conduction_vector0[i] * ncl_kt[i,] *delta_k 
	Polarization_current_length <- Polarization_current_length + elementary_charge/Planck_constant_bar* energy_difference_vector0[i]* 2* Im(mpl_kt[i,] * Conj(d_vc[i])) *delta_k
	Energy_free_length <- Energy_free_length + (energy_valence_vector0[i] * nvl_kt[i,] + energy_conduction_vector0[i] * ncl_kt[i,] ) *delta_k
	Energy_interaction_length <- Energy_interaction_length -elementary_charge* E_vector *(d_vv[i] * nvl_kt[i,] +d_cc[i] * ncl_kt[i,] +2*Re(Conj(d_vc[i])*mpl_kt[i,])) *delta_k
}

Current_length <- Current_valence_length + Current_conduction_length
E_HHG_length <- Current_length + Polarization_current_length
Energy_length <- Energy_free_length + Energy_interaction_length

E_spectrum <- fft(E_vector)
E_HHG_length_spectrum <- fft(E_HHG_length)
Current_length_spectrum <-  fft(Current_length)
Polarization_current_length_spectrum <- fft(Polarization_current_length)
}

if(1){
dir.exists <- function(d) {
    de <- file.info(d)$isdir
    ifelse(is.na(de), FALSE, de)
}
print("Exporting fpW.gif ...")
numt <- round( 2*10^-15 /delta_t )
nplim=c(-1.2, 1.2)
library(animation)
saveGIF({
  ani.options(loop = TRUE)
  for (i in seq(1,(div_num_t+1), by=numt)) {
    .main = paste(t_vector[i]*10^15,"(fs)", seq="")
    plot(ka_vector0, fcW_kt[,i], xlim=c(-pi,pi), ylim=nplim, xlab="BZ", ylab= "numbers", main=.main,col="green", type="l")
    par(new=T)
    plot(ka_vector0, fvW_kt[,i], xlim=c(-pi,pi), ylim=nplim, xlab="", ylab= "", main=.main, col="blue", type="l")
    par(new=T)    
    plot(ka_vector0, Re(pW_kt[,i]), xlim=c(-pi,pi), ylim=nplim, xlab="BZ", ylab= "polarization", main=.main,col="red" , type="l")
  }
}, movie.name="fp.gif", interval=0.2)
}

print("warnings: ")
warnings()
}

if(1){
dir.exists <- function(d) {
    de <- file.info(d)$isdir
    ifelse(is.na(de), FALSE, de)
}
print("Exporting contour.gif ...")
numt <- round( 2*10^-15 /delta_t )
nplim=c(-1.2, 1.2)
library(animation)
saveGIF({
  ani.options(loop = TRUE)
  for (i in seq(1,(div_num_t+1), by=numt)) {
    .main = paste(t_vector[i]*10^15,"(fs)", seq="")
    contour(Re(V_kt[,,i]), main=.main)
  }
}, movie.name="contour.gif", interval=0.2)
}

print("warnings: ")
warnings()
}

