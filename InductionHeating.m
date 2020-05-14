%%    INDUCTION HEATING 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Workpiece in Stainless Steel X5CrNi 18/9 (1.4301)
f = 10000;                    %[Hz]
w = 2*pi*f;
Tf=950;                       %[Gradi Celsius]
Tamb=20;                      %[Gradi Celsius]
vk = 77;                      %[V] presasolo l'ampiezza
I = 40*exp(1i*w*t);           %[A]
mu = 1.256637e-6;             %[H/m]
stb=5.670374e-8;              %[W*m^-2*K^-4] %Stefan-Boltzmann
%deltaH=-stb*(Tf^4-Tamb^4);

      %Sigma Costants
a=4.9659e-7;
b=8.4121e-10;
c=-3.7246e-13;
d=6.1960e-17;

      %Lambda costants
g=0.11215;
q=1.4087e-4;

sigma=1/(a+b*Tf+c*Tf^2+d*Tf^3);   %[]
lambda=100*(g+q*Tf);              %[W*m^-1*K^-1]
R=15e-3;                          %[m]
R0=0;

%Electromagnetic prob
%Studio con il problema elettromagnetico 0<r<10R
%R raggio del workpiece, 10R un punto lontano dal workpiece
%Esprimo il N come numero di Nodi tra 0 e R
N=20;
% Problem Definitions 
%bvpfunE=@(r)(sigma*vk/(2*pi*r));
phi=EM_Eq(sigma, w, vk, N, R, R0, mu);
