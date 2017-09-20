function int=integrarS(fx,Dx)
	%integral de una serie de datos por medio del método de simpson con distancias variables
	n=length(fx);
	if n>=3
		x=fx(1:3);
		dx=Dx(1:3);
		fx(1:2)=[];
		Dx(1:2)=[];
		int=integralLocal(x,dx)+integrarS(fx,Dx);
	else
		int=integralLocal(fx,Dx)+0;
    end
end
	
	
function int=integralLocal(fx,Dx)
    int=0;
	if length(fx)==3
		h1=Dx(2)-Dx(1);
		h2=Dx(3)-Dx(2);
		y0=fx(1);
		y1=fx(2);
		y2=fx(3);
		int=(h1+h2)/(6*h1*h2)*(h2^2*(y1-y0)+h1^2*(y1-y2)+2*h1*h2*(y0+y1+y2));
    elseif length(fx)==2
		h=Dx(2)-Dx(1);
		y0=fx(1);
		y1=fx(2);
		int=(y0+y1)*h/2;
    else
        int=0;
	end
end