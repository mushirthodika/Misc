%-----------ionization potential-------
 Ip=load('diffenergies.dat');

%---parameters for pulse--------
     E=0.0360;
 omega=0.08;
  tmin=2.0;
 tmax=10000.0;
  Ton=9000.0;
    n=1000000;
    e=2.718;
%--------------------------------

%*******Pulse*********************************************************
%------- if t < T_on ------------         

 t=linspace(tmin,Ton,floor((Ton/tmax)*n));
 Fe(1:floor((Ton/tmax)*n))=(E).*((sin(0.5*pi*t/Ton)).^2);  %.*sin(omega*t);   

%--------t > T_on ----------------
 T=linspace(Ton+t(2),tmax,floor(((tmax-Ton)/tmax)*n));
 t(floor(((Ton)/tmax)*n)+1:n)=T;
 Fe(floor(((Ton)/tmax)*n)+1:n)=E;   %(E/(omega*omega)).*sin(omega*T) ;
%--------------------------------

 k=(2.*Ip).^(0.5);
 Ne=2;

%--------------------------------

for i1=1:size(Ip) ;

 gamma=(k(i1).*omega)./Fe;

 g=(1.5./gamma).*( (1+(0.5./gamma.^2) ).*asinh(gamma) - ((1+gamma.^2).^(0.5))./(2.*gamma) ) ;

 Left = ((3.*Fe)/(pi.*(k(i1).^3))).*( (2./k(i1)) -1 );

 Middle = (4.*e.*(k(i1).^3))./(((2./k(i1)) - 1).*(Fe)) ;

 Right = exp(-(2.*g.*(k(i1).^3))./(3.*Fe) ) ;

 FCADK = Ne.*((Left).^(0.5)).*(Fe./(8*pi)).*((Middle).^(2/k(i1))).*Right;

 Y(i1)=1-exp(-trapz(t,FCADK));

end
