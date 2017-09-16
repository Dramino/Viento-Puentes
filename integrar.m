function int=integrar(fx,Dx)
	%integral de una serie de datos por medio del método de simpson
	n=length(fx);
    int=0;
    for i=1 : n
		if i==1 || i==n 
			int=int+fx(i);
		else
			if mod(i,2)==0
				int=int+4*fx(i);
			else
				int=int+2*fx(i);
			end
		end
    end
    int=int*Dx/3;
end