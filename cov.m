function covx=cov(fx,Dx)
	%Regresa la autocovarianza de una serie de datos fx para un tiempo de desfase Dx
	%Dx debe ser entero
	N    = length(fx);
	j    = Dx;
	x1   = fx(1:N-j);
	x2   = fx(j+1:N);
	covx = sum((x1-mean(x1)) .* (x2-mean(x2)))/(N-j);
end