clear all
close all
CapturaDeDatos;
puente.integrarFi;
puente.integrarFiV;
%Valor del modo al extremo del puente en direccion Z
fiI = puente.fiZ(1) ;
wZ=puente.w(2);
mZ=puente.m2;

sRz=0.130;

%Velocidad resonante
Vr = ( puente.dV * wZ ) /(2* pi* puente.St );

%Coeficiente de amortiguamiento aerodinámico
%Constantes
kAz=0.2;
aZ=0.4;
zAe=((viento.rho*puente.B^2*kAz)/(4*mZ)).*((1-(sRz./(aZ*puente.dV)).^2)*(puente.intFi2/puente.intFi2V));

%función de respuesta H de frecuencias
H = abs (1 - ( viento.w ./ wZ ) .^2 + 2*1i *( puente.c - zAe ) *(( viento.w ./ wZ ))).^( -1) ;

%Espectro de carga
lambdaZ=2;
sQz=1;
bZ=0.15;
wS = (2* pi* puente.St *Vr)/puente.dV ;
sQ = 2* lambdaZ * puente.dV * (((0.5* viento.rho *Vr ^2* puente.B* sQz) ^2) /( sqrt (pi)* wS * bZ )) * exp( -(((1 - viento.w ./ wS ) ./ bZ ) .^2)) .* puente.intFi2V ;

%Espectro de respuesta
espRes = (( fiI^2* H .^2) ./(( wZ ^2* mZ * puente.intFi2) .^2) ) .* sQ ;

%Respuesta al extremo del puente
t = linspace (0 ,600 ,600) ;
n=length([viento.w]);
dw = ( viento.w (n)- viento.w(1)) / n ;
rz = zeros (1 ,length(t));
for i = 1: n
	rz = rz + sqrt (2* espRes(i)* dw ) .* cos ( viento.w(i)*t - rand (1) *2* pi);
end

%Valor máximo
maxZ = max (abs (rz));
kpz = maxZ / sRz

%Gráfica de la respuesta en el tiempo
plot (t , rz)
grid
xlabel ( 'T [s] ')
ylabel ( ' r_z [m] ')