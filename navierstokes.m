%This is Navier-Stokes in space!!



%Galerkin Method for Stokes

clear all
close all

nv=6:2:12;
error=zeros(3,size(nv,2));



%exact solution
 u1= @(x,y) (cos(pi*x)+1).*sin(2*pi*y);
 u2= @(x,y) (0.5)*sin(pi*x).*(1-cos(2*pi*y));
 p=@(x,y) sin(pi*x).*cos(pi*y);
 
 f1=@(x,y) (5*(pi^2))*cos(pi*x).*sin(2*pi*y)+4*(pi^2)*sin(2*pi*y)+pi*cos(pi*x).*cos(pi*y)+...
     u1(x,y).*(-pi*sin(pi*x).*sin(2*pi*y))+u2(x,y).*(2*pi*(cos(pi*x)+1).*cos(2*pi*y));
 f2=@(x,y) (0.5)*(-5*(pi^2)*sin(pi*x).*cos(2*pi*y)+(pi^2)*sin(pi*x))-pi*sin(pi*x).*sin(pi*y)+...
     u1(x,y).*((0.5)*pi*cos(pi*x).*(1-cos(2*pi*y)))+u2(x,y).*(pi*sin(pi*x).*sin(2*pi*y));
 
 
     
 for iterate=1:size(nv,2) 
    n=nv(iterate);
    n1=n+1;
    N=n+10;
     
     %Calculation of Stiffness Matrix
    k=(0:n-2)';
    S=spdiags(4*k+6,0,n-1,n-1);
    
    %Calculation of Mass matrix
     v1=-2./(2*k+5);
     M=spdiags([v1,2./(2*k+1)+2./(2*k+5),[0;0;v1(1:n-3)]],[-2,0,2],(n-1),(n-1));
     
     %Laplacian
     A=kron(M,S)+kron(S,M);
     
     %construction of B and C matrices
     P=spdiags(-2*ones(n-1,1),1,n-1,n-1);
     
     j=(0:n+1)';
     lnorm=2./(2*j+1);
     Q=spdiags([lnorm(1:n-1),-lnorm(3:n+1)],[0,2],n-1,N+1);
     
     Qt=Q(:,1:n-1);
     
     B1=-kron(Qt,P);
     B1=B1(:,2:end); %delete column 1
     
     B2=-kron(P,Qt);
     B2=B2(:,2:end); %delete column 1
     
     B=[B1;B2];
     
     
     C1=B1';
     C2=B2';
     
     s1=size(A);
     r1=size(B1,2);
     
     %calculated
     n12=(n-1)^2;
     
     %non-linear matrices W1 and W2
     
     [W1,W2]=nlmatrices(n);
     In12=eye(n12);
     
     
     
     %construction of right hand side
     
     %Legendre-GL points
     lnodes=zeros(N+1,1);
     lnodes(2:N)=gnodesjacobi(N-1,1,1);
     lnodes(1)=-1;
     lnodes(N+1)=1;
     
     [x,y]=meshgrid(lnodes,lnodes);
     x=x';
     x=x(:);
     y=y';
     y=y(:);
     
    
     %forward legendre transform
     f1xy=f1(x,y);
     f2xy=f2(x,y);
     
     fcof1=zeros(N+1); %coef of InF w.r.t Legendre polynomials
     fcof2=fcof1;
     
     w=(2/(N*(N+1)))*(legendreP(N,lnodes).^(-2));
     
     %loop starts
     %gamma
     g=zeros(N+1,1);
     j=(0:N-1)';
     g(1:N)=2./(2*j+1);
     g(N+1)=2/N;
     
     for k=1:N+1
         
         lky=legendreP(k-1,lnodes);
         gk=g(k);
            
         for l=1:N+1
             
             llx=legendreP(l-1,lnodes);
             gl=g(l);
             
            fcof1(l,k)=sum(f1xy.*kron(lky.*w,llx.*w))/(gl*gk);
            fcof2(l,k)=sum(f2xy.*kron(lky.*w,llx.*w))/(gl*gk);
             
         end
     end
     
     fcof1=reshape(fcof1,(N+1)^2,1);
     fcof2=reshape(fcof2,(N+1)^2,1);
     
     F1=kron(Q,Q)*fcof1;
     F2=kron(Q,Q)*fcof2;
     
     F=[F1;F2;zeros((n-1)^2-1,1)];
     
     
     %solving the non-linear system
     sol=zeros(3*(n-1)^2-1,1);
     iter=1;
     while iter<=100
         uhold=sol(1:n12);
         vhold=sol(n12+1:2*n12);
         W=kron(In12,uhold')*W1+kron(In12,vhold')*W2;
         GM=[W+A zeros(s1) B1; zeros(s1) W+A B2; C1 C2 zeros(r1)];
         solnew=GM\F;
         if norm(solnew-sol,inf)<1e-13
             break;
         end
         iter=iter+1;
         sol=solnew;
     end
     
     
     
     
     uh=sol(1:n12); %(n-1)^2
     
     vh=sol(n12+1:2*n12); %(n-1)^2
     
     ph=sol(2*n12+1:end); %(n-1)^2-1
     
     %exact
     
     %Legendre-GL points
     lnodes=zeros(n+1,1);
     lnodes(2:n)=gnodesjacobi(n-1,1,1);
     lnodes(1)=-1;
     lnodes(n+1)=1;
     
     [x,y]=meshgrid(lnodes,lnodes);
     x=x';
     x=x(:);
     y=y';
     y=y(:);
     
     ue=u1(x,y);
     ve=u2(x,y);
     
     %backward legendre transforms 
     
     %u
     [in3,in4]=meshgrid(0:n-2,0:n-2);
     in3=in3';
     in3=in3(:);
     in4=in4';
     in4=in4(:);
     
     uhe=zeros((n+1)^2,1);
     vhe=uhe;
     
     for k=1:(n+1)^2
         
        uhe(k)=sum(uh.*(legendreP(in3,x(k))-legendreP(in3+2,x(k))).*(legendreP(in4,y(k))-legendreP(in4+2,y(k))));
        vhe(k)=sum(vh.*(legendreP(in3,x(k))-legendreP(in3+2,x(k))).*(legendreP(in4,y(k))-legendreP(in4+2,y(k))));
        
     end
     
     %p
     gnodes=gnodesjacobi(n-2,0,0);
     
     [x,y]=meshgrid(gnodes,gnodes);
     x=x';
     x=x(:);
     y=y';
     y=y(:);
     
     pe=p(x,y);
     pe=pe-ones(size(pe))*pe(1,1);
     pe=pe(2:end);
     
     in3=in3(2:end);
     in4=in4(2:end);
     
     n21=(n-2)^2-1;
     phe=zeros(n21,1);
     
     for k=1:n21
         phe(k)=sum(ph.*(legendreP(in3,x(k+1))).*(legendreP(in4,y(k+1))));
     end
     
     c=phe(1)-pe(1);
     phe=phe-c*ones(size(phe));
     
     error(:,iterate)=[norm(uhe-ue);norm(vhe-ve);norm(phe-pe)];
     
 end
figure(1)
semilogy(nv,error(1,:),'+-',nv,error(2,:),'b-*',nv,error(3,:),'r-o');


function [W1,W2]=nlmatrices(n)
     m=(0:n-2)';
     ni=(0:n-2)';
     [m,ni]=meshgrid(m,ni);
     m=m';
     m=m(:);
     ni=ni';
     ni=ni(:);
     n12=(n-1)^2;
     
     W1=zeros(n12^2,n12);
     for k=1:n12
             [Nl,DNl]=nonlinearmatrices(n,m(k),ni(k)); %N-n D-m
             W1((k-1)*n12+1:k*n12,:)=kron(Nl,DNl);
     end
     
     W2=zeros(n12^2,n12);
     for k=1:n12
             [Nl,DNl]=nonlinearmatrices(n,ni(k),m(k)); %N-m D-n
             W2((k-1)*n12+1:k*n12,:)=kron(DNl,Nl);
     end
     
end

function [Nl,DNl]=nonlinearmatrices(N,m,n) %N-n D-m
Nl=zeros(N-1);
for i=0:N-2
    for j=0:N-2
        Nl(i+1,j+1)=tp(i,j,n)-tp(i,j,n+2)-tp(i,j+2,n)+tp(i,j+2,n+2)-...
            (tp(i+2,j,n)-tp(i+2,j,n+2)-tp(i+2,j+2,n)+tp(i+2,j+2,n+2));
    end
end
DNl=zeros(N-1);
for i=0:N-2
    for j=0:N-2
        DNl(i+1,j+1)=tdp(i,j,m)-tdp(i,j,m+2)-tdp(i,j+2,m)+tdp(i,j+2,m+2)-...
            (tdp(i+2,j,m)-tdp(i+2,j,m+2)-tdp(i+2,j+2,m)+tdp(i+2,j+2,m+2));
    end
end
end
function d=tdp(k,l,m)
d=0;
  for j=l-1:-2:0
      d=d+(2*j+1)*tp(k,j,m);
  end
end
function p=tp(k,l,m)
p=0;
if (mod(k+l+m,2)==0 && abs(k-l)<=m ) && m<=k+l
    s=(k+l+m)/2;
    p=(2/factorial(2*s+1))*factorial(2*s-2*k)*factorial(2*s-2*l)*factorial(2*s-2*m)*...
     power(factorial(s)/(factorial(s-k)*factorial(s-l)*factorial(s-m)),2);
end
end