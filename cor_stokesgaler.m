%Galerkin Method for Stokes

clear all
close all
clc
nv=(6:2:16);

nn=size(nv,2);
%nv=(4:10);
 error=zeros(3,size(nv,2));
 u1= @(x,y) (cos(pi*x)+1).*sin(2*pi*y);
 u2= @(x,y) (0.5)*sin(pi*x).*(1-cos(2*pi*y));
 p=@(x,y) sin(pi*x).*cos(pi*y);
 f1=@(x,y) (5*(pi^2))*cos(pi*x).*sin(2*pi*y)+4*(pi^2)*sin(2*pi*y)+pi*cos(pi*x).*cos(pi*y);
 f2=@(x,y) (0.5)*(-5*(pi^2)*sin(pi*x).*cos(2*pi*y)+(pi^2)*sin(pi*x))-pi*sin(pi*x).*sin(pi*y);
 
 for iter=1:nn
     n=nv(iter);
    n1=n+1;
    %calculation of Stiffness matrix
    k=(0:n-2)';
    S=spdiags(4*k+6,0,n-1,n-1);
    %calculation of mass matrix
    v1=-2./(2*k+5);
    M=spdiags([v1,2./(2*k+1)+2./(2*k+5),[0;0;v1(1:n-3)]],[-2,0,2],(n-1),(n-1));
    A=(kron(S,M)+kron(M,S));
    %Legendre-GL points
     lnodes=gnodesjacobi(n-1,1,1);
     lnodes=[-1;lnodes;1];
     [x,y]=meshgrid(lnodes,lnodes);
     x=x';
     x=x(:);
     y=y';
     y=y(:);
     %forward legendre transform
     f1xy=f1(x,y);
     f2xy=f2(x,y);
     fcof1=zeros(n1,n1); %coef of InF w.r.t Legendre polynomials
     fcof2=fcof1;
     w=(2/(n*(n+1)))*(legendreP(n,lnodes).^(-2)); %weights
     for k=1:n+1
         lky=legendreP(k-1,lnodes);
         if k==n+1
                 gk=2/n;
             else
                gk=2/(2*k-1);
         end
         for l=1:n+1
             llx=legendreP(l-1,lnodes);
             if l==n+1
                 gl=2/n;
             else
                gl=2/(2*l-1);
             end
            fcof1(l,k)=sum(f1xy.*kron(lky.*w,llx.*w))/(gl*gk);
            fcof2(l,k)=sum(f2xy.*kron(lky.*w,llx.*w))/(gl*gk);
         end
     end
      %constructor F vector
     j=(0:n)';
     lnorm=2./(2*j+1);
     k1=1:n-1;
     k2=3:n1;
     n11=[(n-1)^2 1];
     F1=reshape(fcof1(k1,k1),n11).*kron(lnorm(k1),lnorm(k1))-reshape(fcof1(k1,k2),n11).*kron(lnorm(k2),lnorm(k1))...
         -reshape(fcof1(k2,k1),n11).*kron(lnorm(k1),lnorm(k2))+reshape(fcof1(k2,k2),n11).*kron(lnorm(k2),lnorm(k2));
     F2=reshape(fcof2(k1,k1),n11).*kron(lnorm(k1),lnorm(k1))-reshape(fcof2(k1,k2),n11).*kron(lnorm(k2),lnorm(k1))...
         -reshape(fcof2(k2,k1),n11).*kron(lnorm(k1),lnorm(k2))+reshape(fcof2(k2,k2),n11).*kron(lnorm(k2),lnorm(k2));
     n12=(n-1)^2;
     F=[F1;F2;zeros(n12-1,1)];
     %construction of B and C matrices
     a=spdiags(-2*ones(n-1,1),1,n-1,n-1);
     b=spdiags([lnorm(1:n-1),-lnorm(1:n-1)],[0,2],n-1,n-1);
     C1=-kron(b',a');
     C2=-kron(a',b');
     B1=C1';
     B2=C2';
     C1=C1(2:n12,:);
     C2=C2(2:n12,:);
     B1=B1(:,2:n12);
     B2=B2(:,2:n12);
     GM=[A zeros(n12) B1; zeros(n12) A B2; C1 C2 zeros(n12-1)];
     sol=GM\F;
     %calculated
     uh=sol(1:n12,1);
     vh=sol(n12+1:2*n12,1);
     ph=sol(2*n12+1:3*n12-1,1);
      %exact
     ue=u1(x,y);
     ve=u2(x,y);
     %backward legendre transforms :(
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
     end
     %v
     for k=1:(n+1)^2
        vhe(k)=sum(vh.*(legendreP(in3,x(k))-legendreP(in3+2,x(k))).*(legendreP(in4,y(k))-legendreP(in4+2,y(k))));
     end
     %p
     lnodes=gnodesjacobi(n-3,1,1);
     lnodes=[-1;lnodes;1];
     [x,y]=meshgrid(lnodes,lnodes);
     x=x';
     x=x(:);
     y=y';
     y=y(:);
     pe=p(x,y);
     pe=pe-ones(size(pe))*pe(1,1);
     pe=pe(2:n12);
     in3=in3(2:end);
     in4=in4(2:end);
     phe=zeros(n12-1,1);
     for k=1:n12-1
         phe(k)=sum(ph.*(legendreP(in3,x(k+1))).*(legendreP(in4,y(k+1))));
     end
     c=phe(1)-pe(1);
     phe=phe-c*ones(size(phe));
     error(:,iter)=[norm(uhe-ue);norm(vhe-ve);norm(phe-pe)];
 end
 figure(1)
semilogy(nv,error(1,:),'+-',nv,error(2,:),'b-*',nv,error(3,:),'r-o');