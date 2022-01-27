function D=globattojacbD(n,al,be)
    n1=n+1; 
    D=zeros(n1);
    N=n-1;
    k=0:N-1;
    a=al+1;
    b=be+1;
    d=(b^2-a^2)./((2.*k+a+b).*(2.*k+a+b+2));
    if a==0 && b==0
        d(1)=0;
    end
    k=1:N-1;
    dg=sqrt((4*k.*(k+a).*(k+b).*(k+a+b))./((2*k+a+b-1).*((2*k+a+b).^2).*(2*k+a+b+1)));
    gnodes=eig(diag(dg,1)+diag(dg,-1)+diag(d));
    lnodes=[-1;gnodes;1];
    %J(x)=d(J(N-1,al+1,be+1));
    %j1=sqrt((n-1)*(n+al+be+2))*jacobiP(n-2,al+2,be+2,lnodes);
    j1=0.5*(n+al+be+2)*jacobiP(n-2,al+2,be+2,lnodes);
    %differentiation matrix
    for i=1:n1
        lnodesi=lnodes(i);
        j1i=j1(i);
        for j=1:n1
            if j==1
                if i==1
                    D(i,j)=(al-n*(n+al+be+1))/(2*(be+2));
                elseif i==n1
                    D(i,j)=(((-1)^n)*gamma(be+2)*gamma(n+al+1))/(2*gamma(al+2)*gamma(n+be+1));
                else
                   D(i,j)=(((-1)^(n-1))*gamma(n)*gamma(be+2)*(1-lnodesi)*j1i)/(2*gamma(n+be+1));
                end
            elseif j==n1
                if i==1
                    D(i,j)=(((-1)^n1)*gamma(al+2)*gamma(n+be+1))/(2*gamma(be+2)*gamma(n+al+1));
                elseif i==n1
                    D(i,j)=(n*(n+al+be+1)-be)/(2*(al+2));
                else
                    D(i,j)=(gamma(n)*gamma(al+2)*(1+lnodesi)*j1i)/(2*gamma(n+al+1));
                end
            else
                if i==1
                    D(i,j)=(2*((-1)^n)*gamma(n+be+1))/(gamma(n)*gamma(be+2)*(1-lnodes(j))*((1+lnodes(j))^2)*j1(j));
                elseif i==n1
                    D(i,j)=(-2*gamma(n+al+1))/(gamma(n)*gamma(al+2)*((1-lnodes(j))^2)*(1+lnodes(j))*j1(j));
                elseif i==j
                    D(i,j)=(al-be+(al+be)*lnodesi)/(2*(1-(lnodesi^2)));
                else
                    D(i,j)=((1-lnodesi^2)*j1i)/((1-lnodes(j)^2)*j1(j)*(lnodesi-lnodes(j)));
                end
            end
        end
    end
end