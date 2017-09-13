classdef Viento

    properties
        rho=0;   %Densidad del viento
        v=0;     %Velocidad media del viento
        h=0;     %Altura a revisar
        w=linspace (0 ,20 ,2000) ;  %frecuencias del viento
    end
    
    methods
        function [sNu,sNw]=densidadEspectralKaimal(obj)
            %Obtiene la función de densidad espectral en todo el dominio w
            %para las dos direeciones
            %Constantes de la funcion
            Au=6.8/(2*pi);
            Aw=9.4/(2*pi);
            H=[obj.h];
            xLu=100*(H/10)^0.3;
            xLw=xLu/12;
            v=[obj.v];
            %Definicion de variables
            w = [obj.w];
            sNu=zeros(1,length(w));
            sNw=zeros(1,length(w));
            %Función de densidad espectral
            sNu= Au*xLu./(v *(1+1.5* Au.*w*(xLu/v )).^(5/3) );
            sNw= Aw*xLw./(v *(1+1.5* Aw.*w*(xLw/v )).^(5/3) );
        end
    end
    
end

