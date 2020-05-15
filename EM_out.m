function [Vo] = EM_out(sigma, w, vk, N, Rc2, Rc1, mu)
% Studio tramite Differenze finite il valore di phi
% tra la soperficie esterna ed interna del coil
%Cambio di variabili s=r^2, V=sqrt(s)*phi
S1=Rc1^2;
S2=Rc2^2;

h=(S2-S1)/(N+1);
hr=(linspace(S1,S2,N+2))';
hmu=4/(h^2*mu);

Afd=zeros(N+1,N+1);

for j = 2:N
Afd(j,j-1) = -hmu*hr(j+1); 
Afd(j,j) = 2*hmu*hr(j+1)+1i*sigma*w; 
Afd(j,j+1) = -hmu*hr(j+1);
end
% Condizioni di Neumann
Afd(1,2)=-2*hmu*hr(2);
Afd(1,1)=2*hmu*hr(2)+1i*w*sigma;
%Condizioni di Dirichlet
Afd(N+1,N)=-hmu*hr(N+2);
Afd(N+1,N+1)=2*hmu*hr(N+2)+1i*sigma*w;

f=zeros(N+1,1);
f(1:N+1)=vk*sigma/(2*pi);
V=Afd\f;
Vo=V(1);
end

