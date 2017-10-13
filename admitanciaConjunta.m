function ji=admitanciaConjuta(x,fi,beta,w)
	n=length(fi);
    ji=0;
    for i=1:n
        for j=1:n
            ji=ji+ fi(i)*fi(j)*exp(-beta*w*abs(x(j)-x(i)));
        end
    end
    ji=ji/n^2;
end