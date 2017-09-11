data <- read.csv("Lorenz.csv")
f_Hz_data <- as.matrix(data[3])
power_data <- as.matrix(data[4])*as.matrix(data[4])

min <- 100000000000000
for (f_Hz in seq(0.9041,0.9042,length=50)){
for (FWHM in seq(0.02775,0.028,length=200)) {
for (power in seq(7350,7400,length=100)){

local_deviation <- power_data - power*FWHM^2/((f_Hz_data - f_Hz)^2 + FWHM^2)

tri <- t(local_deviation) %*% local_deviation
if (  tri[1][1] < min ){
 min <- tri[1][1]
 f_Hz_opt <- f_Hz
 FWHM_opt <- FWHM
power_opt <- power
}

}
}
}

min
f_Hz_opt
FWHM_opt
2*FWHM_opt
power_opt