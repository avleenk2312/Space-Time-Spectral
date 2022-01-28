% Stokes Jacobi Collocation


clc
close all
clear all

%Square Domain
%variables

fig=0;
nv=(4:3:10);
%nv=4;
%nv=(4:10);
al=-0.5;
be=-0.5;

err=zeros(3,size(nv,2));
err2=err;

normstwo=zeros(6,size(nv,2));
normsinf=normstwo;
ms=zeros(1,size(nv,2));
c=ms;
r=c;
rA=r;
realsp=r;
exA=zeros(2,size(nv,2));
exS=exA;
rexA=r;
rexS=r;
c1=zeros(6,size(nv,2));
u1= @(x,y,t) (cos(pi*x)+1).*sin(2*pi*y).*sin(0.5*pi*t);
u2= @(x,y,t) (0.5)*sin(pi*x).*(1-cos(2*pi*y)).*sin(0.5*pi*t);
p= @(x,y,t) sin(pi*x).*cos(pi*y).*sin(0.5*pi*t);
 
f1= @(x,y,t) ((5*(pi^2))*cos(pi*x).*sin(2*pi*y)+4*(pi^2)*sin(2*pi*y)+pi*cos(pi*x).*cos(pi*y)).*sin(0.5*pi*t)...
     +0.5*pi*(cos(pi*x)+1).*sin(2*pi*y).*cos(0.5*pi*t);
 
f2= @(x,y,t) ((0.5)*(-5*(pi^2)*sin(pi*x).*cos(2*pi*y)+(pi^2)*sin(pi*x))-pi*sin(pi*x).*sin(pi*y)).*sin(0.5*pi*t)...
     +(0.25)*pi*sin(pi*x).*(1-cos(2*pi*y)).*cos(0.5*pi*t);
 
u10= @(x,y) -(cos(pi*x)+1).*sin(2*pi*y);
u20= @(x,y) -(0.5)*sin(pi*x).*(1-cos(2*pi*y));
counter=1;

 for iter=1:size(nv,2)
     n=nv(iter);
     n1=n+1;
     Dl=globattojacbD(n,al,be);
     D=Dl(2:n1,2:n1);
     d=Dl(2:n1,1);
     Dl=Dl^2;
     Dl=Dl(2:n,2:n);
     Dcg=cgaussjacbD(n,al,be);
     Dcg=Dcg^2;
     Dcg=Dcg(2:n1,2:n1);
     I=eye(n*(n-1));
     In=eye(n);
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
      
      S=C1*(A1\B1)+C2*(A2\B2);
       m=(n^2)*(n-1);
      A=[A1 zeros(m) B1; zeros(m) A2 B2; C1 C2 zeros(n^3-n)];
   %    issymmetric(S)
   %    eig(S)
   %    meshgrid
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
  %%%%%%%%%% method 1 Direct i)
    %X=C1*(A1\B1)+C2*(A2\B2);
    Y=C1*(A1\F1)+C2*(A2\F2);
    ap=S\Y;
    u1h=A1\(F1-B1*ap);
    u2h=A2\(F2-B2*ap);
    err(:,iter)=[norm(u1h-u1e,inf);norm(u2h-u2e,inf);norm(ap-pe,inf)];
%     m=svd(S);
%     ms(iter)=m(n^2-1);
    c(iter)=cond(S);
    c1(:,iter)=[cond(A1);cond(A2);cond(B1);cond(B2);cond(C1);cond(C2)];
    eigS=eig(S);
    exS(:,iter)=[max(abs(eigS));min(abs(eigS))];
    e=real(eigS);
    rexS=min(e);
    r(iter)=sum(e==0);
    
    %%%%%%%%%%% method Direct ii)
    b=[F1;F2;zeros(n^3-n,1)];
    sol=A\b;
    k=size(F1);
    u1h=sol(1:k);
    u2h=sol(k+1:2*k);
    ap=sol(2*k+1:size(sol));
    err2(:,iter)=[norm(u1h-u1e,inf);norm(u2h-u2e,inf);norm(ap-pe,inf)];
    figure(counter)
    eigA=eig(A);
    plot(eigA,'x');
    title(['Spectrum of A for n=' num2str(nv(iter)) ' for alpha=' num2str(al) ' and beta=' num2str(be)]);
    counter=counter+1;
    eA=real(eigA);
    realsp(iter)=isreal(eigA);
    rA(iter)=sum(eA==0);
    rexA(iter)=min(eA);
    exA(:,iter)=[max(abs(eigA));min(abs(eig(A)))];
    
    figure(counter)
    plot(eig(S),'.');
    title(['Spectrum of S for n=' num2str(nv(iter)) ' for alpha=' num2str(al) ' and beta=' num2str(be)]);
    counter=counter+1;
    
    normstwo(:,iter)=[norm(A1);norm(A2);norm(B1);norm(B2);norm(C1);norm(C2)];
    normsinf(:,iter)=[norm(A1,inf);norm(A2,inf);norm(B1,inf);norm(B2,inf);norm(C1,inf);norm(C2,inf)];
    
 end
figure(counter)
semilogy(nv,err(1,:),'+-',nv,err(2,:),'r-*',nv,err(3,:),'b-o');
title(['Plot of error for alpha=' num2str(al) ' and beta=' num2str(be)]);
legend('u1','u2','p');
xlabel('n');
ylabel('error');
% figure(2)
% semilogy(nv,ms,'-o',nv,1./sqrt(nv),'-*');
counter=counter+1;
figure(counter)
loglog(nv,c,'b-o',nv, c1(1,:),'r-o', nv, c1(3,:), 'g-o', nv, c1(5,:), 'm-o',nv, nv.^2,'k-x',nv,nv.^4,'c-x');
legend('S','A','B','C','n^2','n^4');
title(['Plot of condition numbers of S, A, B and C, for alpha=' num2str(al) ' and beta=' num2str(be)]);
xlabel('n');
ylabel('condition number');




