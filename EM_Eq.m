function [phi] = EM_Eq(sigma, w, vk, N, R, mu, V1)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% S1=R*10;
% No=N*10;
% S0=R0^2;
% S=R^2;
%% All'esterno del workpiece
% ho=(S1-S)/(No+1);
% hro=(linspace(S,S1,No+1))';
% hmu=4/(mu*ho^2);
% eo=ones(No+1,1);
% Afdo=spdiags([hmu*eo -2*hmu*eo hmu*eo],...
%     -1:1, No+1, No+1);
% for l=1:No+1
%     Afdo(l,1:No+1)=Afdo(l,1:No+1)*sqrt(hro(l));
% end
% 
% Afdo(1,1)=(4*sqrt(hro(1)))/(ho^2)*(1+1/mu);
% Afdo(1,2)=(-2*hmu*sqrt(hro(1)));
% %Afdo(No+2,No+1)=2*hmu*sqrt(hro(No+2));
% %Afdo(No+2,No+2)=(-2*hmu*sqrt(hro(No+2))-4/(mu*ho^2*sqrt(hro(No+2))));
% g=zeros(No+1,1);
% g(No+1,1)=-hmu;
% Vo=Afdo\g;
% 
% %Interno del Workpiece
% h=(S-S0)/(N+1);
% c=1;
% %hr=(linspace(S0,S,N+1))';
% e=ones(N,1);
% Afd=spdiags([(-4)/(h^2*mu)*e 1*e ...
%     ((-4)/(h^2*mu))*e], -1:1, N, N);
% for l=2:1:N-1
%     Afd(l,c)=Afd(l,c)*l*h;
%     c=c+1;
% end
% for l=1:1:N-1
%     Afd(l,l)=((8*l*h)/(mu*h^2)+(1i*w*sigma));
% end 
% c=2;
% for l=1:1:N-1
%     Afd(l,c)=Afd(l,c)*l*h;
%     c=c+1;
% end
% %Robin Condition
% Afd(N,N)=((8*(N)*h)/(mu*h^2)+(1i*w*sigma));
% Afd(N,N-1)=-4*(N)*h/(mu*h^2);
% %Dirichlet Condition
% Afd(1,1)=((8*h)/(h^2*mu)+1i*w*sigma);
% Afd(1,2)=(-4*h/(mu*h^2));
% e=ones(N,1);
% f=e*vk*sigma/(2*pi);
% f(N)=f(N)+Vo(1)*4*N*h/(h^2*mu);
% V=Afd\f;
% phi=zeros(length(V)+1,1);
% for l=2:N+1
%     phi(l)=V(l-1)/sqrt(h*(l-1));
% end

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

%Condizioni
Afd(1,1)=2*hmu*hr(2)+1i*w*sigma;
Afd(1,2)=-1*hmu*hr(2);

% Afd(No+1,No)=-2*hmu*hr(No+2);
% Afd(No+1,No+1)=2*hmu*hr(No+2)+1i*w*sigma;
% Afd(No+2,No+2)=-2*hmu*sqrt(hr(No+3));
% Afd(No+2,No+3)=2*hmu*sqrt(hr(No+3));

Afd(N,N-1)=2*hmu*sqrt(hr(N+1));
Afd(N,N)= 2*hmu*hr(N+1)+1i*sigma*w; 

% 
% e=[1:N];
% e1=[1:N-1];
% e2=[1:N-1];
% 
% e1(1:No)=((-hmu)*hr(3:No+1));
% e1(No+1:N)=((hmu)*sqrt(hr(No+2:N)));
% 
% 
% e1(1:No)=((-hmu)*hr(3:No+1));
% e1(No+1:N)=((hmu)*sqrt(hr(No+2:N)));
% 
% e(1:No+1)=((2*hmu*hr(2:No+2))/(mu*h^2)+...
%          (1i*w*sigma));
% e(No+2:N)=((-2*hmu*sqrt(hr(No+3:N+1))));
% 
% e2(1:No+1)=((-hmu*hr(2:No+2)));
% e2(No+2:N-1)=((hmu*sqrt(hr(No+3:N))));
% 
% Afd=diag([e1*ones(1,N)],-1)+diag([e*ones(1,N)])+diag([e2*ones(1,N)],1);


% c=No+1;
% for l=No+2:N
%      Afd(l,c)=(Afd(l,c)*(hmu)*sqrt(hr(l)));
%      c=c+1;
% end
%  for l=No+2:N
%      Afd(l,l)=(-2*hmu*sqrt(l));
%  end 
%  c=No+3;
% for l=No+2:N-1
%     Afd(l,c)=(Afd(l,c)*(hmu*sqrt(hr(l))));
%     c=c+1;
% end


% f=ones(N,1);
% f(1:No+1)=f(1:No+1)*vk*sigma/(2*pi);
% f(No+2:N-1)=0;
% f(N)=(-hmu*sqrt(hr(N+1))/sqrt(hr(N+2)));
% V=e2\f;
%  for l=2:N+1
%      phi(l)=V(l-1)/sqrt(h*(l-1));
%  end


% S0=R0^2;
% S=R^2;
% h=(S-S0)/(N+1);
% c=1;
% %hr=(linspace(S0,S,N+1))';
% e=ones(N,1);
% Afd=spdiags([(-4)/(h^2*mu)*e 1*e ...
%     ((-4)/(h^2*mu))*e], -1:1, N, N);
% 
% for l=2:1:N
%     Afd(l,c)=Afd(l,c)*l*h;
%     c=c+1;
% end
% for l=1:1:N
%     Afd(l,l)=((8*l*h)/(mu*h^2)+(1i*w*sigma));
% end 
% c=2;
% for l=1:1:N-1
%     Afd(l,c)=Afd(l,c)*l*h;
%     c=c+1;
% end
% % %Robin Condition
% % Afd(N+1,N+1)=((8*(N+1)*h)/(mu*h^2)+(1i*w*sigma)+4/(mu*h));
% % Afd(N+1,N)=-8*(N+1)*h/(mu*h^2);
% %Dirichlet Condition
% % Afd(1,1)=((8*h)/(h^2*mu)+1i*w*sigma);
% % Afd(1,2)=(-4*h/(mu*h^2));
f=zeros(N,1);
f(1:N-1)=vk*sigma/(2*pi);
f(N)=vk*sigma/(2*pi)+V1*hmu*hr(N+1);
V=Afd\f;

phi=zeros(length(V)+1,1);
for l=2:length(V)+1
    phi(l)=V(l-1)/sqrt(hr(l));
end

end

