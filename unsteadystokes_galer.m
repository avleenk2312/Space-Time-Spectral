%Galerkin Method P_N-P_{N-2} for the unsteady-Stokes equations
%chebyshev gauss Lobatto in time!

clearvars
close all
clc
nv=(6:2:16);
%nv=(4:10);
%nv=20;
nn=size(nv,2);

 error=zeros(3,size(nv,2));

u1= @(x,y,t) (cos(pi*x)+1).*sin(2*pi*y).*sin(0.5*pi*t);
u2= @(x,y,t) (0.5)*sin(pi*x).*(1-cos(2*pi*y)).*sin(0.5*pi*t);
p= @(x,y,t) sin(pi*x).*cos(pi*y).*sin(0.5*pi*t);
 
f1= @(x,y,t) ((5*(pi^2))*cos(pi*x).*sin(2*pi*y)+4*(pi^2)*sin(2*pi*y)+pi*cos(pi*x).*cos(pi*y)).*sin(0.5*pi*t)...
     +0.5*pi*(cos(pi*x)+1).*sin(2*pi*y).*cos(0.5*pi*t);
 
f2= @(x,y,t) ((0.5)*(-5*(pi^2)*sin(pi*x).*cos(2*pi*y)+(pi^2)*sin(pi*x))-pi*sin(pi*x).*sin(pi*y)).*sin(0.5*pi*t)...
     +(0.25)*pi*sin(pi*x).*(1-cos(2*pi*y)).*cos(0.5*pi*t);
 
 for iter=1:nn
    n=nv(iter);
    n1=n+1;
    %Legendre-GL points
     lnodes=gnodesjacobi(n-1,1,1);
     lnodes=[-1;lnodes;1];
     
     clnodes=zeros(n+1,1);
     clnodes(2:n)=gnodesjacobi(n-1,0.5,0.5);
     clnodes(1)=-1;
     clnodes(n1)=1;
     
     [x,y,t]=meshgrid(lnodes,lnodes,clnodes(2:n1));
     [x,y]=meshgrid(lnodes,lnodes);
     x=x';
     x=x(:);
     y=y';
     y=y(:);
    L=spdiags([-1*ones(n-1,1),ones(n-1,1)],[-2,0],n+1,n-1);
    L=kron(L,L);
    u10=u1(x,y,-1);
    u20=u2(x,y,-1);
    %forward legendre transform
    [x,y,t]=meshgrid(lnodes,lnodes,clnodes(2:n1));
     x=permute(x,[2 1 3]);
     x=x(:);
     y=permute(y,[2 1 3]);
     y=y(:);
     t=t(:);  
     f1xyt=f1(x,y,t);
     f2xyt=f2(x,y,t);
     fcof1=zeros(n1,n1,n); %coef of InF w.r.t Legendre polynomials
     fcof2=fcof1;
     u1cof=zeros(n1,n1);
     u2cof=u1cof;
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
            fc1=f1xyt.*kron(ones(n,1),kron(lky.*w,llx.*w));
            fcof1(l,k,:)=(sum(reshape(fc1,(n1)^2,n))/(gl*gk))';
            fc2=f2xyt.*kron(ones(n,1),kron(lky.*w,llx.*w));
            fcof2(l,k,:)=(sum(reshape(fc2,(n1)^2,n))/(gl*gk))';
            u1cof(l,k)=sum(u10.*kron(lky.*w,llx.*w))/(gl*gk);
            u2cof(l,k)=sum(u20.*kron(lky.*w,llx.*w))/(gl*gk);
         end
     end
     %finding u1 and u2 coefficients
     u1c=L\(reshape(u1cof,(n+1)^2,1));
     u2c=L\(reshape(u2cof,(n+1)^2,1));
    %constructor F vector
     j=(0:n)';
     lnorm=2./(2*j+1);
     k1=1:n-1;
     k2=3:n1;
     set1=kron(ones(n,1),kron(lnorm(k1),lnorm(k1)));
     set2=kron(ones(n,1),kron(lnorm(k2),lnorm(k1)));
     set3=kron(ones(n,1),kron(lnorm(k1),lnorm(k2)));
     set4=kron(ones(n,1),kron(lnorm(k2),lnorm(k2)));
     n11=[((n-1)^2)*n 1];
     F1=reshape(fcof1(k1,k1,:),n11).*set1-reshape(fcof1(k1,k2,:),n11).*set2...
         -reshape(fcof1(k2,k1,:),n11).*set3+reshape(fcof1(k2,k2,:),n11).*set4;
     F2=reshape(fcof2(k1,k1,:),n11).*set1-reshape(fcof2(k1,k2,:),n11).*set2...
         -reshape(fcof2(k2,k1,:),n11).*set3+reshape(fcof2(k2,k2,:),n11).*set4;
     n12=((n-1)^2)*n;
     %time derivative
    Dl=globattojacbD(n,-0.5,-0.5);
    D=Dl(2:n1,2:n1);
    d=Dl(2:n1,1);
    In=eye(n);
    %calculation of Stiffness matrix
    k=(0:n-2)';
    S=spdiags(4*k+6,0,n-1,n-1);
    %calculation of mass matrix
    v1=-2./(2*k+5);
    M=spdiags([v1,2./(2*k+1)+2./(2*k+5),[0;0;v1(1:n-3)]],[-2,0,2],(n-1),(n-1));
    A=kron(S,M)+kron(M,S);
    M2=kron(M,M);
    At=kron(D,M2)+kron(In,A);

    %F vector
    Dt=kron(d,M2);
    F1=F1-Dt*u1c;
    F2=F2-Dt*u2c;
    F=[F1;F2;zeros(n12-n,1)];
     %construction of B and C matrices
     a=spdiags(-2*ones(n-1,1),1,n-1,n-1);
     b=spdiags([lnorm(1:n-1),-lnorm(1:n-1)],[0,2],n-1,n-1);
     n123=(n-1)^2;
     C1=-kron(b',a');
     C2=-kron(a',b');
     C1=C1(2:n123,:);
     C2=C2(2:n123,:);
     C1=kron(In,C1);
     C2=kron(In,C2);
     B1=C1';
     B2=C2';
 
    GM=zeros(3*n12-n);
     GM=[At zeros(n12) B1; zeros(n12) At B2; C1 C2 zeros(n12-n)];
     sol=GM\F;
    
     %calculated
     k=size(F1);
     u1h=sol(1:k);
     u2h=sol(k+1:2*k);
     ap=sol(2*k+1:size(sol));
      %backward legendre transforms
      %u
       ph=ap(((n-1)^2-1)*(n-1)+1:end);
      [x,y]=meshgrid(lnodes,lnodes);
      x=x';
      x=x(:);
      y=y';
      y=y(:);
      u1e=u1(x,y,1);
      u2e=u2(x,y,1);
      [in3,in4]=meshgrid(0:n-2,0:n-2);
     in3=in3';
     in3=in3(:);
     in4=in4';
     in4=in4(:);
     u1he=zeros((n+1)^2,1);
     u2he=u1he;
     uh=u1h((n-1)^3+1:end);
     vh=u2h((n-1)^3+1:end);
     for k=1:(n+1)^2
        u1he(k)=sum(uh.*(legendreP(in3,x(k))-legendreP(in3+2,x(k))).*(legendreP(in4,y(k))-legendreP(in4+2,y(k))));
        u2he(k)=sum(vh.*(legendreP(in3,x(k))-legendreP(in3+2,x(k))).*(legendreP(in4,y(k))-legendreP(in4+2,y(k))));
     end
     %p
     gnodes=gnodesjacobi(n-2,0,0);
     [x,y]=meshgrid(gnodes,gnodes);
     x=x';
     x=x(:);
     y=y';
     y=y(:);
     pe=p(x,y,1);
     pe=pe-ones(size(pe))*pe(1,1);
     pe=pe(2:end);
     in3=in3(2:end);
     in4=in4(2:end);
     phe=zeros(size(pe));
     for k=1:((n-2)^2-1)
         phe(k)=sum(ph.*(legendreP(in3,x(k+1))).*(legendreP(in4,y(k+1))));
     end
     c=phe(1)-pe(1);
     phe=phe-c*ones(size(phe));
     error(:,iter)=[norm(u1he-u1e);norm(u2he-u2e);norm(phe-pe)];
 end
figure(1)
loglog(nv,error(1,:),'+-',nv,error(2,:),'b-*',nv,error(3,:),'r-o');
