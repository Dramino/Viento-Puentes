classdef Viento

    properties
        rho=0;   %Densidad del viento
        v=0;     %Velocidad media del viento
        h=0;     %Altura a revisar
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
            V=[obj.v];
 
            %Función de densidad espectral
            sNu=@(w) Au.*xLu./(V *(1+1.5* Au.*w*(xLu/V )).^(5/3) );
            sNw=@(w) Aw.*xLw./(V *(1+1.5* Aw.*w*(xLw/V )).^(5/3) );
        end
    end
    
end

