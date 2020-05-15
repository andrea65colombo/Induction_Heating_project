%%    INDUCTION HEATING 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Workpiece in Stainless Steel X5CrNi 18/9 (1.4301)
f = 10000;              %[Hz]
w = 2*pi*f;
Tf=950;                 %[Gradi Celsius]
Tamb=20;                %[Gradi Celsius]
vk = 77*exp(1i*w*100);      %[V]%100 al posto di t
I = 40*exp(1i*w*100);       %[A]%100 al posto di t
mu = 1.256637e-6;       %[H/m]
stb=5.670374e-8;        %[W*m^-2*K^-4] %Stefan-Boltzmann
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
R=20e-3;                   %[m]
R0=0;
Rc1=25e-3;
Rc2=30e-3;
d=5e-3;                   %[m] diametro filo coil

%Electromagnetic prob
%Studio con il problema elettromagnetico 0<r<R
%R raggio del workpiece 
N=15;
% Problem Definitions 
%bvpfunE=@(r)(sigma*vk/(2*pi*r));
Vo=EM_out(sigma, w, vk, N, Rc2, Rc1, mu);
V1=EM_in(N, Rc1, R, mu, Vo);
phi=EM_Eq(sigma, w, vk, N, R, mu, V1);
phi
% %Allungo phi per valori negativi di r
% % -R<r>R per usare condizioni al contorno sui 2 bordi
% phi=zeros(2*length(phi0)-1,1)';
% phi0=phi0';
% phi(1:length(phi0))=fliplr(phi0);
% phi(length(phi0):end)=phi0;
% phi=phi';
% % 
% % %Heat
% % %Studio l'equazione del calore per valori -R<r<R
% N=N*2+1;
% [hr]=BVP(R, R1, N);
% for l=1:(N+2)
%     bvpfun(l)=((sigma/2)*(abs(-1i*w*phi(l)+vk/(2*pi*hr(l))).^2));
% end
% [T] = HeatEquation(lambda, N, deltaH,...
% R, R1, bvpfun);
% 
% 
% %Grafico da rivedere
% plot(T,'o');



