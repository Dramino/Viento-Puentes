clear all
close all
CapturaDeDatos;
t = cputime;
puente.integrarFi;
%Amortiguamiento aerodinámico
%Definición de integrales
intY   = 0;
intFiY = 0;
intZ   = 0;
intFiZ = 0;

%Dominio de la frecuencia del viento
w=linspace (1 ,2000 ,2000);
intY   = sum(puente.fiY.^2.*puente.D.* puente.cD.*puente.L);
intFiY = sum(puente.fiY.^2.*puente.L);
intZ   = sum(puente.fiZ.^2.*puente.B.*puente.cdl + puente.D.*puente.cD.*puente.L);
intFiZ = sum(puente.fiZ.^2.*puente.L);

%Definición de amortiguamiento aerodinámico
zetaAeY = -(viento.rho*viento.v) /(2* puente.w(1)* puente.m1) * intY / intFiY ;
zetaAeZ = -(viento.rho*viento.v) /(4* puente.w(2)* puente.m2) * intZ / intFiZ ;

%función de respuesta H de frecuencias
Hy = @(w) abs (1 - ( w ./ puente.w(1) ) .^2 + 2*1i *(puente.c - zetaAeY) *(( w ./ puente.w(1) ))) .^( -1) ;
Hz = @(w) abs (1 - ( w ./ puente.w(2) ) .^2 + 2*1i *(puente.c - zetaAeZ) *(( w ./ puente.w(2) ))) .^( -1) ;


%Funcion de densidad espectral de Kaimal
[sNu,sNw] = viento.densidadEspectralKaimal;

%Co espectro - Constantes
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
e = cputime-t;
dMatriz=zeros(1,169);
% xi=puente.L;
% xa=xi;
% xa(:,length(xi))=[];
% xj=[0,xa];
% dx=abs(xi-xj);
	% Ingtegral
	for i    = 1: n
        t = cputime;
		Di   = puente.D(i);
		Cdi  = puente.cD(i);
		Cdli = puente.cdl(i);
		L    = puente.L(i);
		xi   = xi + L ;
		fiYi = puente.fiY(i);
        fiZi = puente.fiZ(i);
		xj   = 0;
		x0   = puente.L(i)
		for j = 1: n
            t = cputime;
			Dj   = puente.D(j);
			Cdj  = puente.cD(j);
			Cdlj = puente.cdl(j);
			L    = puente.L(j);
			xj   = xj + L ;
			fiYj = puente.fiY(j);
            fiZj = puente.fiZ(j);
						
            %Co-espectro
            CoU = @(w,Dx) exp (- cuY .* ( w.*Dx ) /(2 * pi * viento.v));
			CoW = @(w,Dx) exp (- cwY .* ( w.*Dx ) /(2 * pi * viento.v));
			Dx  = arrayfun(@(x abs(x-x0),puente.L);
            
			jY = @(w) puente.fiY .* puente.fiY .* ( (((2* puente.Iu )/ puente.B) ^2* Di * Cdi * Dj * Cdj ) .* CoU(w,Dx).* sNu(w) + ( puente.cL * puente.Iv ) ^2 .* CoW(w,Dx).* sNw(w) ) * L * L ;
			jZ = @(w) jZ + puente.fiZ .* puente.fiZ .* ( (2* puente.cL* puente.Iu) ^2 .* CoU(w).* sNu(w) + (( puente.Iv^2*( Cdli +( Di * Cdi )/puente.B) *(Cdlj +( Dj * Cdj )/ puente.B)) .*CoW(w).* sNw(w) )) * L * L ;
            e = cputime-t;
        end
        e = cputime-t;
        dMatriz(i)=Dx;
    end

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