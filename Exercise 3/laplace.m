function d=laplace(d)
    function neighbor(di,dk)
        if d.omega(i+di,k+dk)
            d.L(r,d.rr(i+di,k+dk))=1;
        end
    end % neighbor
        d.L=sparse(d.NU);
        for r=1:d.NU
            d.L(r,r)=-4;
            i=d.ii(r);
            k=d.kk(r);
            neighbor(1,0);
            neighbor(-1,0);
            neighbor(0,1);
            neighbor(0,-1);
        end
    h=d.x(2)-d.x(1);
   % d.L=d.L/h^2;
end % laplace