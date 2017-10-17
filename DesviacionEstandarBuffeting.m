%obtenci贸n de la desviaci贸n est谩ndar

clear all
close all
CapturaDeDatos;
puente.integrarFi;
%Amortiguamiento aerodin谩mico
%Definici贸n de integrales
intSy   = 0;
intFiY = 0;
intSz   = 0;
intFiZ = 0;

intSy   = integrarS(puente.fiY.^2.*puente.D.* puente.cD,puente.L);
intFiY = integrarS(puente.fiY.^2.,puente.L);
intSz   = integrarS(puente.fiZ.^2.*puente.B.*puente.cdl + puente.D.*puente.cD,puente.L);
intFiZ = integrarS(puente.fiZ.^2.,puente.L);

%Definici贸n de amortiguamiento aerodin谩mico
zetaAeY = -(viento.rho*viento.v) /(2* puente.w(1)* puente.m1) * intSy / intFiY ;
zetaAeZ = -(viento.rho*viento.v) /(4* puente.w(2)* puente.m2) * intSz / intFiZ ;

%funci贸n de respuesta H de frecuencias
Hy = abs (1 - ( viento.w ./ puente.w(1) ) .^2 + 2*1i *(puente.c - zetaAeY) *(( viento.w ./ puente.w(1) ))) .^( -1) ;
Hz = abs (1 - ( viento.w ./ puente.w(2) ) .^2 + 2*1i *(puente.c - zetaAeZ) *(( viento.w ./ puente.w(2) ))) .^( -1) ;


%Funcion de densidad espectral de Kaimal
[sNu,sNw] = viento.densidadEspectralKaimal;

%Co espectro - Constantes
coU = zeros(1 ,length(viento.w));
coW = zeros(1 ,length(viento.w));
cuY = 9;
cwY = 6;

%Funci贸n de aceptancia normalizada
jY = 0;
jZ = 0;
%Expansi贸n y prolongaci贸n de datos
div = 4;
puente.D   = puente.prolong ( puente.D   , div );
puente.cD  = puente.prolong ( puente.cD  , div );
puente.cdl = puente.prolong ( puente.cdl , div );
puente.fiY = puente.prolong ( puente.fiY , div );
puente.fiZ = puente.prolong ( puente.fiZ , div );
puente.L   = puente.expand  ( puente.L   , div );

n = length(puente.L);
x=cumsum(puente.L);
jY=0;
jZ=0;

% funci髇 de acepantcia conjunta
for i=1:n
	for j=1:n
        Dx=abs(x(i)-x(j));
        coU = exp (- cuY* ( viento.w*Dx ) /(2 * pi * viento.v));
        coW = exp (- cwY* ( viento.w*Dx ) /(2 * pi * viento.v));
		jY  = jY + puente.fiY(i) * puente.fiY(j) .* ( (((2* puente.Iu )/ puente.B) ^2* puente.D(i) * puente.cD(i) * puente.D(j) * puente.cD(j) ) .* coU.* sNu + ( puente.cL * puente.Iv ) ^2 .* coW.* sNw ) * puente.L(i) * puente.L(j) ;
        jZ  = jZ + puente.fiZ(i) * puente.fiZ(j) .* ( (2* puente.cL* puente.Iu) ^2 .* coU.* sNu + (( puente.Iv^2*( puente.cdl(i) +( puente.D(i) * puente.cD(i) )/puente.B) *(puente.cdl(j) +( puente.D(j) * puente.cD(j) )/ puente.B)) .*coW.* sNw )) * puente.L(i) * puente.L(j) ;
	end
end

% Normalizaci贸n de la funci贸n de aceptancia
jNormY = jY /( [puente.intFi1] ) ^2;
jNormZ = jZ /( [puente.intFi2] ) ^2;
%D omega
Nomega = length ( viento.w);
Domega = ( viento.w ( Nomega ) - viento.w(1) ) /( Nomega -1) ;
%Integral para la desviaci贸n
%a=sum(Hy.^2.*jNormY)
%b=sum(Hz.^2.*jNormZ)
%Calibraci髇 con integral de trapecio de matlab
%a2=trapz(Hy.^2.*jNormY,viento.w);
%b2=trapz(Hz.^2.*jNormZ,viento.w);
intSy = integrar(Hy.^2.*jNormY,Domega);
intSz = integrar(Hz.^2.*jNormZ,Domega);
%desviaci贸 esd谩ndar
yr = 1;
sigmaY=abs(puente.fiY(yr))*viento.rho *viento.v^2*puente.B/(2*puente.m1*puente.w(1)^2)*sqrt(intSy)

zr=1;
sigmaZ=abs(puente.fiZ(zr))*viento.rho *viento.v^2*puente.B/(2*puente.m2*puente.w(2)^2)*sqrt(intSz)