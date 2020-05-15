function [phi] = EM_Eq(sigma, w, vk, N, R, mu, V1)
% Studio tramite Differenze finite il valore di phi
% tra 0 e la superficie del workpiece

S=R^2;

h=S/(N+1);
hr=(linspace(0,S,N+2))';
hmu=4/(h^2*mu);

Afd=zeros(N,N);

for j = 2:N-1
Afd(j,j-1) = -hmu*hr(j+1); 
Afd(j,j) = 2*hmu*hr(j+1)+1i*sigma*w; 
Afd(j,j+1) = -hmu*hr(j+1);
end

%Condizioni di Dirichlet
Afd(1,1)=2*hmu*hr(2)+1i*w*sigma;
Afd(1,2)=-1*hmu*hr(2);
Afd(N,N-1)=2*hmu*sqrt(hr(N+1));
Afd(N,N)= 2*hmu*hr(N+1)+1i*sigma*w; 

f=zeros(N,1);
f(1:N-1)=vk*sigma/(2*pi);
f(N)=vk*sigma/(2*pi)+V1*hmu*hr(N+1);
V=Afd\f;

%Torno alla variabile phi
phi=zeros(length(V)+1,1);
for l=2:length(V)+1
    phi(l)=V(l-1)/sqrt(hr(l));
end
end

