%%    INDUCTION HEATING 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Workpiece in Stainless Steel X5CrNi 18/9 (1.4301)
%Physics Constants
f = 10000;                   %[Hz]
w = 2*pi*f;
Tf=950;                      %[Gradi Celsius]
Tamb=20;                     %[Gradi Celsius]
vk = 77*exp(1i*w*100);       %[V]%100 al posto di t
I = 40*exp(1i*w*100);        %[A]%100 al posto di t
mu = 1.256637e-6;            %[H/m]
stb=5.670374e-8;             %[W*m^-2*K^-4] %Stefan-Boltzmann
deltaH=-stb*(Tf^4-Tamb^4);
      %Sigma Costants
a=4.9659e-7;
b=8.4121e-10;
c=-3.7246e-13;
d=6.1960e-17;
      %Lambda costants
g=0.11215;
q=1.4087e-4;

sigma=1/(a+b*Tf+c*Tf^2+d*Tf^3);   %[?]
lambda=100*(g+q*Tf);              %[W*m^-1*K^-1]
Rw=20e-3;                          %[m] Raggio del Workpiece
Rc1=25e-3;                        %[m] Raggio interno del Workpiece
Rc2=30e-3;                        %[m] Raggio esetrno del Workpiece
d=5e-3;                           %[m] Diametro cavo coil

%% Problema Elettromagnetico

phi=EM_Eq(sigma, w, vk, mu, Rw, Rc1, Rc2);
phi
plot(abs(phi));

