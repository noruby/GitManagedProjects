e_xx <- 1.5
e_yy <- 1.6
e_zz <- 1.4
e_yz <- e_zy <- 0.1
e_zx <- e_xz <- 0.2
e_xy <- e_yx <- 0.3
 
u <- e_zz
v <- (e_xx+e_yy)*e_zz -e_xz^2 -e_yz^2
w <- e_xx*e_yy*e_zz -e_xx*e_yz^2 -e_yy*e_zx^2 -e_zz*e_xy^2 +2*e_xy*e_yz*e_zx

p1 <- (v+sqrt(v^2-4*u*w))/(2*u)
p2 <- (v-sqrt(v^2-4*u*w))/(2*u)

a1 <- e_xy*e_yz -e_xz*(e_yy-p1)
b1 <- e_yx*e_xz -e_yz*(e_xx-p1)
a2 <- e_xy*e_yz -e_xz*(e_yy-p2)
b2 <- e_yx*e_xz -e_yz*(e_xx-p2)

k10 <- sqrt(p1)
k20 <- sqrt(p2)

r_xx <- 2/(a1*b2-a2*b1)*(a1*b2/(1+k10)-a2*b1/(1+k20))-1
r_xy <- 2*b1*b2/(a1*b2-a2*b1)*(1/(1+k10)-1/(1+k20))
r_yy <- 2/(b1*a2-b2*a1)*(b1*a2/(1+k10)-b2*a1/(1+k20))-1
r_yx <- 2*a1*a2/(b1*a2-b2*a1)*(1/(1+k10)-1/(1+k20))
r_xx

#fresnel coefficient
n <- sqrt(e_xx) 
r_fresnel <- (n-1)/(n+1)
r_fresnel 
