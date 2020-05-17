function [V1] = EM_air(N, Rc1, R, mu, Vo)
% Studio tramite Differenze finite il valore di phi
% tra la superficie interna del coil e la superficie del workpiece

S1=Rc1^2;
S=R^2;
h=(S1-S)/(N+1);
hr=(linspace(S1,S,N+2))';
hmu=4/(h^2*mu);
Afd=zeros(N+1,N+1);

for j = 2:N
Afd(j,j-1)=hmu*sqrt(hr(j+1));
Afd(j,j)=-2*hmu*sqrt(hr(j+1));
Afd(j,j+1)=hmu*sqrt(hr(j+1));
end

%Condizioni di Neumann
Afd(1,2)=2*hmu*sqrt(hr(2));
Afd(1,1)=-2*hmu*sqrt(hr(2));
%Condizioni di Dirichlet
Afd(N+1,N)=hmu*sqrt(hr(N+2));
Afd(N+1,N+1)=-2*hmu*sqrt(hr(N+2));

f=zeros(N+1,1);
f(N+1)=-Vo*hmu*sqrt(hr(N+2));
V=Afd\f;
V1=V(1);
end

