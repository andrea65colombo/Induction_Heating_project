function [phi] = EM_Eq(sigma, w, vk, No, R, R0, mu)
%Risoluzione tramite differenze finite del problema elettromagnetico
%Il problema viene studiato al variare della del raggio
%nel dominio 0<r<10*R 
%R è il raggio del workpiece cilindrico preso in considerazione e 10R un punto
%sufficientemente lontano dal coil da poter applicare le condizioni al contorno


%Cambio di variabili s=r^2
S1=(R*10)^2;
N=((No+1)*10)-1;
S0=R0^2;
S=R^2;
%Discretizzazzione
h=S1/(N+1);
hr=(linspace(0,S1,N+2))';
hmu=4/(h^2*mu);

%Creazione Matrice 

Afd=zeros(N+1,N+1);

for j = 2:No+1 %fino No se utilizzo condizioni interfaccia
Afd(j,j-1) = -hmu*hr(j+1); 
Afd(j,j) = 2*hmu*hr(j+1)+1i*sigma*w; 
Afd(j,j+1) = -hmu*hr(j+1);
end

for j = No+2:N %da No+3 se utilizzo condizioni d'interfaccia
Afd(j,j-1)=hmu*sqrt(hr(j+1));
Afd(j,j)=-2*hmu*sqrt(hr(j+1));
Afd(j,j+1)=hmu*sqrt(hr(j+1));
end

%Condizioni al contorno di Dirichlet

Afd(1,1)=2*hmu*hr(2)+1i*w*sigma;
Afd(1,2)=-1*hmu*hr(2);

%Condizioni di interfaccia sulla superficie del workpiece

% Afd(No+1,No)=-2*hmu*hr(No+2);
% Afd(No+1,No+1)=2*hmu*hr(No+2)+1i*w*sigma;
% Afd(No+2,No+2)=-2*hmu*sqrt(hr(No+3));
% Afd(No+2,No+3)=2*hmu*sqrt(hr(No+3));

%Condizione al contorno di Robin lontano dal workpiece

Afd(N+1,N)=2*hmu*sqrt(hr(N+2));
Afd(N+1,N+1)=-2*hmu*sqrt(hr(N+2))+(-4)/(mu*h*sqrt(hr(N+2)));

f=zeros(N+1,1);
f(1:No+1)=vk*sigma/(2*pi);
V=Afd\f;

%Cambio di variabili da V a Phi

phi=zeros(length(V),1);
for l=1:length(V)
    phi(l)=V(l)/sqrt(hr(l+1));
end

end

