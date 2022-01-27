function gnodes=gnodesjacobi(n,al,be)
if al==be && al==-0.5
         k=0:n-1;
         gnodes=-cos(((2.*k+1)*pi)./(2*n));
         gnodes=gnodes';
else
         k=0:n-1;
         d=(be^2-al^2)./((2.*k+al+be).*(2.*k+al+be+2));
         if al==0 && be==0
             d(1)=0;        
         end
         k=1:n-1;
         dg=sqrt((4*k.*(k+al).*(k+be).*(k+al+be))./((2*k+al+be-1).*((2*k+al+be).^2).*(2*k+al+be+1)));
         gnodes=eig(diag(dg,1)+diag(dg,-1)+diag(d));
end
end