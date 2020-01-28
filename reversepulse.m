 format long

%---parameters for pulse--------
     E=5.635e-04;
 omega=0.01;
  tmin=0.0;
tmax=5000.0;
  Ton=4000.0;
    n=40000;
%--------------------------------

%*******Pulse*********************************************************
%------- if t < T_on ------------         

 t=linspace(tmin,Ton,floor((Ton/tmax)*n));
 alpha(1:floor((Ton/tmax)*n))=(E/(omega*omega)).*((sin(0.5*pi*t/Ton)).^2).*sin(omega*t);   

%--------t > T_on ----------------
 T=linspace(Ton+t(2),tmax,floor(((tmax-Ton)/tmax)*n));
 t(floor(((Ton)/tmax)*n)+1:n)=T;
 alpha(floor(((Ton)/tmax)*n)+1:n)=(E/(omega*omega)).*sin(omega*T) ;

%*********************************************************************

%******Eigenvector of field free hamiltonian*************************
 n1=1001;
 x=linspace(-6.00,6.00,n1);
 dt=t(2)-t(1);
 

 KE=dlmread('KE_MAT.txt');
 [vec,val]=eig(KE);
 EXPT=diag(exp(-1i*(-1i*dt)*diag(val)));

 %L=@(x,alpha) 0.005224852*(x-alpha).^4 - 0.0139753572*(x-alpha).^2;
 
 %L = dlmread('v0kh.dat');

 P1=diag(hcndoublewell(x,0.00));
 H=KE+P1;
 [EVEC,EVAL]=eig(H);

%********************************************************************

 
%********************************************************************

 b=0;
%*****************Propagation Part**************************************

 EVEC1=EVEC(:,1);
 EVEC2(:,1)=x';
 V(:,1)=x';

 for i1=1:n ;
 
  D=hcndoublewell(x,alpha(i1));

 EVEC1=diag(exp(-D*1i*(-1i*dt)*0.50))*EVEC1 ;  
 EVEC1=(vec')*EVEC1;
 EVEC1=EXPT*EVEC1;
 EVEC1=vec*EVEC1;
 EVEC1=diag(exp(-D*1i*(-1i*dt)*0.50))*EVEC1;  
 
 norm = EVEC1'*EVEC1;               % Orthogonalisation
 EVEC1 = (1.0/norm).*EVEC1;         

%------Animation section--------------------------

 EVEC2(:,1+b)=conj(EVEC1).*EVEC1;
     V(:,1+b)=D;
 b=b+1
 subplot (3,1,1)
 plot (t(1:i1),alpha(1:i1))
 axis([tmin tmax -E/(omega*omega)-0.2 E/(omega*omega)+0.2 ])
 
 subplot (3,1,2)
 plot (x,V(:,b));
 axis ([-5.5 5.5 -0.25 0.5]);

 subplot (3,1,3)
 plot (x',EVEC2(:,b));
 axis ([-6.0 6.0 -0.02 0.1]);
 pause(0.01)
%-------------------------------------------------
 end

h=figure(1);
writerObj = VideoWriter('rt2_b1.avi');
open(writerObj);
f(1:100)=struct('cdata',[],'colormap',[]);

for i1=1:20:n
 subplot (3,1,1)
 plot (t(1:i1),alpha(1:i1))
 axis([tmin tmax -E/(omega*omega)-0.2 E/(omega*omega)+0.2 ])

 subplot (3,1,2)
 plot (x,V(:,i1+1));
 axis ([-5.5 5.5 -0.25 0.5]);

 subplot (3,1,3)
 plot (x',EVEC2(:,i1+1));
 axis ([-6.0 6.0 -0.02 0.1]);
 title('omega=0.05 & E=0.0140875')
 f(i1)=getframe(gcf);
  writeVideo(writerObj,f(i1));
end
close(writerObj);
%---------------------------------------------------------
%save('rpulse_data.mat')
%dlmwrite('evec.mat',EVEC2(:,1:1000:n),'delimiter','\t','precision',8)
