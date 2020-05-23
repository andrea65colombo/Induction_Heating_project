function [phi] = EM_Eq(sigma, w, vk, mu, Rw, Rc1, Rc2)
% Studio tramite Differenze finite il valore di phi

S=(Rw*10)^2;
Sw=Rw^2;
Sc1=Rc1^2;
Sc2=Rc2^2;
a=Rw*1e3;
b=Rc1*1e3;
c=Rc2*1e3;
d=Rw*10*1e3;

N=Rw*10*1e3;

h=S/N;
hr=(linspace(0,S,N+1))';
hmu=4/(h^2*mu);

Afd=zeros(N,N);

for j = 2:(a-1)
Afd(j,j-1) = -hmu*hr(j+1); 
Afd(j,j) = 2*hmu*hr(j+1)+1i*sigma*w; 
Afd(j,j+1) = -hmu*hr(j+1);
end

for j = (a+1):(b-1)
Afd(j,j-1)=hmu*sqrt(hr(j+1));
Afd(j,j)=-2*hmu*sqrt(hr(j+1));
Afd(j,j+1)=hmu*sqrt(hr(j+1));
end

for j = (b+1):(c-1)
Afd(j,j-1) = -hmu*hr(j+1); 
Afd(j,j) = 2*hmu*hr(j+1)+1i*sigma*w; 
Afd(j,j+1) = -hmu*hr(j+1);
end

for j = (c+1):(d-1)
Afd(j,j-1)=hmu*sqrt(hr(j+1));
Afd(j,j)=-2*hmu*sqrt(hr(j+1));
Afd(j,j+1)=hmu*sqrt(hr(j+1));
end

%Condizioni Contorno Dirichlet in 0
Afd(1,1)=2*hmu*hr(2)+1i*w*sigma;
Afd(1,2)=-1*hmu*hr(2);
%Condizioni d'interfaccia
Afd(a,a-1)=-1;
Afd(a,a)=2;
Afd(a,a+1)=-1;
Afd(b,b-1)=-1;
Afd(b,b)=2;
Afd(b,b+1)=-1;
Afd(c,c-1)=-1;
Afd(c,c)=2;
Afd(c,c+1)=-1;

%Robin condition
Afd(N,N)=(-2*hmu*sqrt(hr(N+1))-4/(mu*h*sqrt(hr(N+1))));
Afd(N,N-1)=(2*hmu*sqrt(hr(N+1)));

f=zeros(N,1);
f(1:a-1)=vk*sigma/(2*pi);
f(a:b)=0;
f(b+1:c-1)=vk*sigma/(2*pi);
f(c:N)=0;
V=Afd\f;
phi=zeros(length(V)+1,1);
for l=2:length(V)+1
    phi(l)=V(l-1)/sqrt(hr(l));
end

end

