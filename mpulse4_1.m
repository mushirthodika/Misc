
% program for pulse

 format long

% fs=(1.05457188*10^(-34))/(43.6*10^(-19))
% au=((43.6*10^(-19)/(1.05457188*10^(-34))

 op=fopen('norm.txt','w');
 %op1=fopen('wave3.txt','w');

% for j1=1:10;    % loop for changing E value

%**********Parameter declaration**************************

  omega=0.05 ;  %*41.324*j1*0.50;

% omega=j1*0.0110*41.324; % j1 iteration for omega value 

  tmin=0.00;
  tmax=10000.00;
  Ton=8000.00;
    tp=2000.00; % already in fs
     n=1000000;
%     w=1000.00;  % already in fs
     E=0.0140875;   %*0.5*0.5*j1*j1;

%**********************************************************
%******Eigenvector of field free hamiltonian***************

H=dlmread('HMAT.txt');
V=dlmread('VMAT.txt');
Z=dlmread('ZMAT.txt');

HA=H;   %-1i*V-Z;
[EVEC,EVAL]=eig(H);

%**********************************************************
%**************cosine^2 envelop****************************
 dt=(tmax-tmin)/(n-1); 
%------- if t < T_on ------------         

 t=linspace(tmin,Ton,floor((Ton/tmax)*n));
 wave(1:floor((Ton/tmax)*n))=(E).*((sin(0.5*pi*t/Ton)).^2).*sin(omega*t);   

%--------t > T_on ----------------
 T=linspace(Ton+t(2),tmax,floor(((tmax-Ton)/tmax)*n));
 t(floor(((Ton)/tmax)*n)+1:n)=T;
 wave(floor(((Ton)/tmax)*n)+1:n)=(E).*sin(omega*T) ;

%*********************************************************************

  tic()
 for i1=1:n ;  %(n-1)/4 ;
  i1
 t= tmin+(i1-1)*dt ; 

% wave=E*sin(omega*(t-tp))*exp(-1.38*((t-tp)/w)^2) ;

% if abs(t-tp) < w ;
% wave=E*cos(omega*(t-tp))*(cos((pi*(t-tp))/(2*w)))^2 ;
% else
% wave=0.00;
% end

% fprintf(op1,'%7.5f %10.8f\n',t,wave);
%**********************************************************

 HA=H-0.50.*1i.*V-Z.*wave(i1);

 EVEC1=expm(-1i.*HA.*dt)*EVEC;
 EVEC=EVEC1;
 norm=(conj(EVEC)*EVEC);

 fprintf(op,'%7.5f %10.8f\n',t,norm(1,1));

 end
  toc()

 fprintf(op,'\n')
 fprintf(op,'\n')
 fprintf(op,'\n')

 %fprintf(op1,'\n')
 %fprintf(op1,'\n')
 %fprintf(op1,'\n')

% end

 fclose(op);
 %fclose(op1);
%*********************************************************

