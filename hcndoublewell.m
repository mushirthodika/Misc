  function f = hcndoublewell(x,alpha);

 a1 = -0.1422d0;b1 = -1.34d0;c1 = 1.826d0;a2 = 0.0009163d0;

 b2 = -1.46d0;c2 = 0.7081d0;a3 = -83.53d0;b3 = 9.541d0;

 c3 = 15.15d0;a4 = -0.04274d0;b4 = 0.5924d0;c4 = 0.8776d0;

 a5 = -0.02282d0;b5 = -0.4495d0;c5 = 1.005d0;a6 = -72.59d0;

 b6 = -10.69d0;c6 = 12.97d0;a7 = -0.01293d0;b7 = 1.429d0;

 c7 = 0.8177d0;a8 = -0.1139d0;b8 = 1.69d0;c8 = 1.789d0;

fn1 = a1*exp(-(((x-alpha)-b1)/c1).^2) + a2*exp(-(((x-alpha)-b2)/c2).^2) + a3*exp(-(((x-alpha)-b3)/c3).^2) + a4*exp(-(((x-alpha)-b4)/c4).^2);
fn2 = a5*exp(-(((x-alpha)-b5)/c5).^2) + a6*exp(-(((x-alpha)-b6)/c6).^2) + a7*exp(-(((x-alpha)-b7)/c7).^2) + a8*exp(-(((x-alpha)-b8)/c8).^2);

f = fn1 + fn2 + 93.1582106d0;


end
% a=gn(x,alpha);
%plot(x,gn(x,alpha));
