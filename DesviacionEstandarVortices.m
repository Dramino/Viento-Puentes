clear all
close all
CapturaDeDatos;
puente.integrarFi;
puente.integrarFiV;
%Valor del modo al extremo del puente en direccion Z
fiI = puente.fiZ(1) ;
wZ  = puente.w(2);
mZ  = puente.m2;

%Integral
fiX    = [datos(1:13,8); datos(33:43,8)];
fiZ    = [datos(1:13,7); datos(33:43,7)];
intVor = sum(( fiX.^2 + fiZ.^2) .* puente.Lv);

%Velocidad resonante
Vr   = (puente.dV * wZ ) /(2* pi* puente.St );
V    = puente.dV*puente.w(2)/(2*pi*puente.St);
z    = ((4* puente.m2 * puente.c ) /( viento.rho *puente.B ^2* viento.Kaz )) *( puente.intFi2 / intVor );
b    = fiI /(2^(5/2) * pi ^(7/4) ) * (( viento.rho *puente.dV ^3* viento.lambdaZ ) /( puente.m2 * viento.bZ * viento.Kaz * puente.intFi2) )^(1/2) * ( viento.sigmaQz /( puente.St ^2* viento.az ) ) * (V / Vr ) ^(3/2) * exp ( -(1/2) *((1 -( Vr /V)) / viento.bZ )^2) ;
sAux = ((1 - z) /2 + (((1 - z ) /2) ^2 + b ^2) ^(1/2) ) ^(1/2) ;
s    = sAux * viento.az*puente.dV