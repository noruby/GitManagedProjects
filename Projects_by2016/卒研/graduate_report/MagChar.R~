
data <- read.table("MCF_15-01-26-18_32_00.VSM",header=T,skip=40,sep=",")
H <- data["H"]
H <- unlist(H)
M <- data["M"]
M <- unlist(M)

par(mar=c(5,5,2,2))
plot(H,M,type="l",xlab=expression(H(Oe)),ylab=expression(M(emu)))
abline(h=0)
abline(v=0)
dev.copy2eps(file="MagCharMCF_cgs.eps")

H <- H/(4*pi)*1000
M <- M/1000
par(mar=c(5,5,2,2))
plot(H,M,type="l",xlab=expression(H(A/m)),ylab=expression(M(A%.%m^2)))
abline(h=0)
abline(v=0)
dev.copy2eps(file="MagCharMCF_SI.eps")

H <- H*(4*pi)*10^(-7)
M <- M
par(mar=c(5,5,2,2))
plot(H,M,type="l",xlab=expression(B[0](T)),ylab=expression(M(A%.%m^2)))
abline(h=0)
abline(v=0)
dev.copy2eps(file="aaa.eps")

data <- read.table("Nisample2_15-01-26-17_41_10.VSM",header=T,skip=40,sep=",")
H <- data["H"]
H <- unlist(H)
M <- data["M"]
M <- unlist(M)

par(mar=c(5,5,2,2))
plot(H,M,type="l",xlab=expression(H(Oe)),ylab=expression(M(emu)))
abline(h=0)
abline(v=0)
dev.copy2eps(file="MagCharNi_cgs.eps")

H <- H/(4*pi)*1000
M <- M/1000
par(mar=c(5,5,2,2))
plot(H,M,type="l",xlab=expression(H(A/m)),ylab=expression(M(A%.%m^2)))
abline(h=0)
abline(v=0)
dev.copy2eps(file="MagCharNi_SI.eps")

#axes=FALSE,
#axis(side=1,pos=0)
#axis(side=2,at=seq(-2,2,0.5),pos=0)
#mtext("a",side=2,line=0)
