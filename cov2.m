function covx=cov2(fx,,px,Dx)
	%Regresa la autocovarianza de una serie de datos fx para un tiempo de
	%desfase Dx en función a la densidad de probabilidad
	%Dx debe ser entero
	j    = Dx;
	x1   = fx(1:N-j);
	x2   = fx(j+1:N);
    fxi    = integrar(px .* x1 .* x2);
	covx = sum((x1-mean(x1)) .* (x2-mean(x2)))/(N-j);
end