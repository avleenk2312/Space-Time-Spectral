%Unsteady Navier-Stokes Jacobi Collocation


clc
close all
clear all

al=0;
be=0;

nv=(6:2:12);
nn=size(nv,2);
error=zeros(3,size(nv,2));

u1= @(x,y,t) (cos(pi*x)+1).*sin(2*pi*y).*sin(0.5*pi*t);
u2= @(x,y,t) (0.5)*sin(pi*x).*(1-cos(2*pi*y)).*sin(0.5*pi*t);
p= @(x,y,t) sin(pi*x).*cos(pi*y).*sin(0.5*pi*t);
 
f1= @(x,y,t) ((5*(pi^2))*cos(pi*x).*sin(2*pi*y)+4*(pi^2)*sin(2*pi*y)+pi*cos(pi*x).*cos(pi*y)).*sin(0.5*pi*t)...
     +0.5*pi*(cos(pi*x)+1).*sin(2*pi*y).*cos(0.5*pi*t)+...
     u1(x,y,t).*(-pi*sin(pi*x).*sin(2*pi*y)).*sin(0.5*pi*t)+...
     u2(x,y,t).*(2*pi*(cos(pi*x)+1).*cos(2*pi*y)).*sin(0.5*pi*t);
 
f2= @(x,y,t) ((0.5)*(-5*(pi^2)*sin(pi*x).*cos(2*pi*y)+(pi^2)*sin(pi*x))-pi*sin(pi*x).*sin(pi*y)).*sin(0.5*pi*t)...
     +(0.25)*pi*sin(pi*x).*(1-cos(2*pi*y)).*cos(0.5*pi*t)+...
     u1(x,y,t).*((0.5)*pi*cos(pi*x).*(1-cos(2*pi*y))).*sin(0.5*pi*t)+...
     u2(x,y,t).*(pi*sin(pi*x).*sin(2*pi*y)).*sin(0.5*pi*t);
 
u10= @(x,y) -(cos(pi*x)+1).*sin(2*pi*y);
u20= @(x,y) -(0.5)*sin(pi*x).*(1-cos(2*pi*y));

 for iter=1:size(nv,2)
     n=nv(iter);
     n1=n+1;
     Dl=globattojacbD(n,al,be);
     D=Dl(2:n1,2:n1); %first derivative row one deleted
     dl=Dl(2:n,2:n);% first derivative row one and end deleted
     d=Dl(2:n1,1);
     Dl=Dl^2;
     Dl=Dl(2:n,2:n); %Dl==second order lobatto derivative
     Dcg=cgaussjacbD(n,al,be);
     dcg=Dcg(2:n1,2:n1);%first derivative row one and end deleted
     Dcg=Dcg^2;
     Dcg=Dcg(2:n1,2:n1); %Dcg==second order closed gauss derivative
     
     
     
     I=eye(n*(n-1));  %I= I_{N(N-1)}
     In=eye(n);        %In
     A1=kron(D,I)-kron(eye(n^2),Dl)-kron(In,kron(Dcg,eye(n-1)));
     A2=kron(D,I)-kron(I,Dcg)-kron(In,kron(Dl,In));
     
     %block for interpolation of p
     Bp=zeros(n-1,n);
     %lobatto nodes
     lnodes=gnodesjacobi(n-1,al+1,be+1);
     lnodes=[-1;lnodes;1];
     %gauss nodes
      gnodes=gnodesjacobi(n,al,be);
      cgnodes=[-1;gnodes;1];
      %jacobi at lobatto nodes
      jacbn=jacobiP(n,al,be,lnodes(2:n));
      djacbn=(0.5*gamma(al+be+n+2)/gamma(n+al+be+1))*jacobiP(n-1,al+1,be+1,gnodes);
      for i=1:n-1
          Bp(i,:)=-jacbn(i)./(djacbn'.*((ones(1,n)*lnodes(i+1)-gnodes').^2));
      end
      
      B1=kron(In,Bp);
      B1=B1(:,2:n^2);
      B1=kron(In,B1);
      B2=kron(Bp,In);
      B2=B2(:,2:n^2);
      B2=kron(In,B2);
      
      %interpolation of u and v
      Bu=zeros(n,n-1);
      for i=1:n
          Bu(i,:)=(be-al-((al+be)*gnodes(i)))./(ones(1,n-1)*gnodes(i)-lnodes(2:n)');
          Bu(i,:)=Bu(i,:)+(1-gnodes(i)^2)./((ones(1,n-1)*gnodes(i)-lnodes(2:n)').^2);
          Bu(i,:)=Bu(i,:)./(jacbn');
          Bu(i,:)=(djacbn(i)/(n*(n+al+be+1)))*Bu(i,:);
      end
      C1=kron(In,Bu);
      C1=C1(2:n^2,:);
      C1=kron(In,C1);
      C2=kron(Bu,In);
      C2=C2(2:n^2,:);
      C2=kron(In,C2);
      
      %interpolation for non-linear terms
      N1=zeros(n,n-1); %lobatto basis at interior gauss nodes
      for i=1:n
          N1(i,:)=((1-gnodes(i)^2)*djacbn(i))./(ones(1,n-1)*gnodes(i)-lnodes(2:n)');
          N1(i,:)=N1(i,:)./(jacbn');
          N1(i,:)=(-1/(n*(n+al+be+1)))*N1(i,:);
      end
      N2=zeros(n-1,n);%closed gauss basis at lobatto nodes
      for j=1:n
          N2(:,j)=(1-(lnodes(2:n).^2)).*jacbn;
          N2(:,j)=N2(:,j)./(lnodes(2:n)-ones(n-1,1)*gnodes(j));
          N2(:,j)=N2(:,j)/((1-gnodes(j)^2)*djacbn(j));
      end
      
      
      
      %exact stuff
  
     [x,y,t]=meshgrid(lnodes(2:n),gnodes,lnodes(2:n1));
     x=permute(x,[2 1 3]);
     x=x(:);
     y=permute(y,[2 1 3]);
     y=y(:);
     t=t(:);
     u1e=u1(x,y,t);
     ind=1:n*(n-1);
     u10e=u10(x(ind),y(ind));
     F1=f1(x,y,t)-kron(d,u10e);
     [x,y,t]=meshgrid(gnodes,lnodes(2:n),lnodes(2:n1));
     x=permute(x,[2 1 3]);
     x=x(:);
     y=permute(y,[2 1 3]);
     y=y(:);
     t=t(:);
     u2e=u2(x,y,t);
     u20e=u20(x(ind),y(ind));
     F2=f2(x,y,t)-kron(d,u20e);
     
     [x,y,t]=meshgrid(gnodes,gnodes,lnodes(2:n1));
     x=permute(x,[2 1 3]);
     x=x(:);
     y=permute(y,[2 1 3]);
     y=y(:);
     t=t(:);
     pe=p(x,y,t);
     pe=pe-kron(diag(pe(1:n^2:n^3)),eye(n^2))*ones(size(pe));
     index=setdiff(1:n^3,1:n^2:n^3);
     pe=pe(index);
  

    %%%%%%%%%%% Direct method 
    
    
    m=(n^2)*(n-1);
    b=[F1;F2;zeros(n^3-n,1)];
    sol=zeros(size(b));
    k=size(F1);
    iterate=1;
    
    while iterate<=100
        uhk=sol(1:k);
        vhk=sol(k+1:2*k);
        W1=diag(uhk)*kron(eye(n^2),dl)+diag(kron(In,kron(N1,N2))*vhk)*kron(In,kron(dcg,eye(n-1)));
        W2=diag(kron(In,kron(N2,N1))*uhk)*kron(I,dcg)+diag(vhk)*kron(In,kron(dl,In));
        A=[W1+A1 zeros(m) B1; zeros(m) W2+A2 B2; C1 C2 zeros(n^3-n)];
        solnew=A\b;
        if norm(solnew-sol,inf)<1e-13
             break;
        end
        iterate=iterate+1;
        sol=solnew;   
    end
        
    u1h=sol(1:k);
    u2h=sol(k+1:2*k);
    ap=sol(2*k+1:size(sol));
    error(:,iter)=[norm(u1h-u1e,inf);norm(u2h-u2e,inf);norm(ap-pe,inf)];
   
    
 end

figure(1)
semilogy(nv,error(1,:),'+-',nv,error(2,:),'r-*',nv,error(3,:),'b-o');
title(['Plot of error for alpha=' num2str(al) ' and beta=' num2str(be)]);
legend('u1','u2','p');
xlabel('n');
ylabel('error');

