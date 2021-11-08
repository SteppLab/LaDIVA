function h=diva_hanning(N1,w,factor)
%h=[zeros(1,N1),1];
%h=[zeros(1,N1),.5+.5*cos(pi*(-N2+1:0)/(N2)),.5+.5*cos(pi*(1:N3-1)/(N3))];
h=[zeros(1,N1),1,w*(1-factor)*factor.^(1:100)];
h=h/sum(h);