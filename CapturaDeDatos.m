%Guarda los datos en una matriz de datos
datos=dlmread('Inputfile.txt','',1,0);

%crea los objetos
puente=Puente;
viento=Viento;
%Asigna los datos a cada objeto
puente.L=datos(:,2);
puente.x=datos(:,3);
puente.Lv=[datos(1:13,2); datos(33:43,2)];
puente.m=datos(:,5);
puente.fiY=datos(:,6);
puente.fiZ=datos(:,7);
puente.fiX=datos(:,8);
puente.w=[1.401 1.687];
puente.c=0.008;
%Obtencion de la masa modal efectiva
puente.mModal;
%Pasar de kg a kN
puente.m1=puente.m1*1000/9.82;
puente.m2=puente.m2*1000/9.82;

%Geometría
puente.B=11.1;
puente.D=datos(:,4);
puente.dV=3.4;

%Coeficiente aoerdinmámicos
puente.cD=datos(:,9);
puente.cL=0.5;
puente.cdl= zeros(length(puente.L),1);
puente.St=0.11;
%Datos del viento
viento.rho=1.25;
viento.v=38.4;
viento.h=30.08;

%Turbulencia
puente.Iu=0.14;
puente.Iv=0.07;