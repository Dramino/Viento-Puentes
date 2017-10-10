clear all
close all
tiempo = cputime;
CapturaDeDatos;
puente.integrarFi;
%Amortiguamiento aerodinámico
%Definición de integrales
intY   = 0;
intFiY = 0;
intZ   = 0;
intFiZ = 0;

intY   = integrarS(puente.fiY.^2.*puente.D.* puente.cD,puente.L);
intFiY = integrarS(puente.fiY.^2.,puente.L);
intZ   = integrarS(puente.fiZ.^2.*puente.B.*puente.cdl + puente.D.*puente.cD,puente.L);
intFiZ = integrarS(puente.fiZ.^2.,puente.L);


%Definición de amortiguamiento aerodinámico
zetaAeY = -(viento.rho*viento.v) /(2* puente.w(1)* puente.m1) * intY / intFiY ;
zetaAeZ = -(viento.rho*viento.v) /(4* puente.w(2)* puente.m2) * intZ / intFiZ ;

%función de respuesta H de frecuencias
Hy = abs (1 - ( viento.w ./ puente.w(1) ) .^2 + 2*1i *(puente.c - zetaAeY) *(( viento.w ./ puente.w(1) ))) .^( -1) ;
Hz = abs (1 - ( viento.w ./ puente.w(2) ) .^2 + 2*1i *(puente.c - zetaAeZ) *(( viento.w ./ puente.w(2) ))) .^( -1) ;


%Funcion de densidad espectral de Kaimal
[sNu,sNw] = viento.densidadEspectralKaimal;

%Co espectro - Constantes
coU = zeros(1 ,length(viento.w));
coW = zeros(1 ,length(viento.w));
cuY = 9;
cwY = 6;

%Función de aceptancia normalizada
jY = 0;
jZ = 0;
%Expansión y prolongación de datos
div = 4;
puente.D   = puente.prolong ( puente.D   , div );
puente.cD  = puente.prolong ( puente.cD  , div );
puente.cdl = puente.prolong ( puente.cdl , div );
puente.fiY = puente.prolong ( puente.fiY , div );
puente.fiZ = puente.prolong ( puente.fiZ , div );
puente.L   = puente.expand  ( puente.L   , div );

n = length(puente.L);
xi=0;
cuPrueba=0
	% Ingtegral
	for i    = 1: n
 		Di   = puente.D(i);
		Cdi  = puente.cD(i);
		Cdli = puente.cdl(i);
		L    = puente.L(i);
		xi   = xi  + L ;
		fiYi = puente.fiY(i);
        fiZi = puente.fiZ(i);
		xj   = 0;
		for j = 1: n
			Dj   = puente.D(j);
			Cdj  = puente.cD(j);
			Cdlj = puente.cdl(j);
			L    = puente.L(j);
			xj   = xj + L ;
			fiYj = puente.fiY(j);
            fiZj = puente.fiZ(j);
			Dx   = abs (xi - xj );
            %Co-espectro
            coU = exp (- cuY* ( viento.w*Dx ) /(2 * pi * viento.v));
            coW = exp (- cwY* ( viento.w*Dx ) /(2 * pi * viento.v));
			cuPrueba=[cuPrueba coU];
			%j
			%jy2= puente.fiY .* puente.fiY .* ( (((2* puente.Iu )/ puente.B) ^2 .* puente.D .* puente.cD .* puente.D .* puente.cD) .* coU.* sNu + ( puente.cL * puente.Iv ) ^2 .* coW.* sNw ) .* puente.L * puente.L ;
            jY = jY + fiYi * fiYj .* ( (((2* puente.Iu )/ puente.B) ^2* Di * Cdi * Dj * Cdj ) .* coU.* sNu + ( puente.cL * puente.Iv ) ^2 .* coW.* sNw ) * L * L ;
			jZ = jZ + fiZi * fiZj .* ( (2* puente.cL* puente.Iu) ^2 .* coU.* sNu + (( puente.Iv^2*( Cdli +( Di * Cdi )/puente.B) *(Cdlj +( Dj * Cdj )/ puente.B)) .*coW.* sNw )) * L * L ;
		end
		x=cumsum(puente.L)-xi;
		cuPrueba2= @(w)exp (- cuY.* ( w.*x ) ./(2* pi * viento.v));
        %jY =@(w)@(z)integrarS fiYi * fiYj .* ( (((2* puente.Iu )/ puente.B) ^2* Di * Cdi * Dj * Cdj ) .* coU.* sNu + ( puente.cL * puente.Iv ) ^2 .* coW.* sNw ) * L * L ;
    end
e = cputime-tiempo
% Normalización de la función de aceptancia
jNormY = jY /( [puente.intFi1] ) ^2;
jNormZ = jZ /( [puente.intFi2] ) ^2;

%Espectro de respuesta
espResDivFiY = (( viento.rho *viento.v ^2* puente.B) /(2* puente.m1 * puente.w(1)^2) ) ^2 * Hy.^2.* jNormY;
espResDivFiZ = (( viento.rho *viento.v ^2* puente.B) /(2* puente.m2 * puente.w(2)^2) ) ^2 * Hz.^2.* jNormZ;
%Respuesta al extremo del puente
espResY = espResDivFiY*puente.fiY(1)^2;
espResZ = espResDivFiZ*puente.fiZ(1)^2;

%Respuesta en el dominio del tiempo
%Divisions de la frecuencia del viento
Nw = length ( viento.w );
dw = ( viento.w( Nw ) - viento.w(1)) / Nw ;
%Simulación de 10 minutos
t   = linspace (0 ,600 ,600) ;
ry  = zeros ( 1 , length(t) );
rz  = zeros ( 1 , length(t) );
sRy = 0.139;
sRz = 0.076;
%Matriz de valores aleatorios
%fiAl= linspace (1 ,1 ,length(viento.w)) ;
%fiAl=map(fiAl,sin);
%función anonima de valores al azar de 0 a d pi
%crearAleatorio=@(x)x*rand(1)*2*pi;
%aleatorio=arrayfun(crearAleatorio,linspace (1 ,1 ,600) );
for i = 1: Nw
	ry = ry+ sqrt(2* espResY(i)* dw ) .* cos ( viento.w(i)*t - rand (1) *2* pi);
	rz = rz+ sqrt(2* espResZ(i)* dw ) .* cos ( viento.w(i)*t - rand (1) *2* pi);
end

%Factor pico
maxY = max ( abs (ry ) );
maxZ = max(abs(rz));
kpx  = maxY / sRy
kpZ  = maxZ / sRz

%Gráfica de la respuesta
subplot (2 ,1 ,1)
plot (t , ry)
grid
xlabel ( 'T [s] ')
ylabel ( ' r_y [m] ')
hold all
subplot (2 ,1 ,2)
plot (t , rz )
grid
xlabel ( 'T [s] ')
ylabel ( ' r_z [m] ')
e = cputime-tiempo