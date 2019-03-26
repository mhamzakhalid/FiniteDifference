function d=domain(N)
d.x=linspace(-3,3,N);
d.y=linspace(-3,3,N);
[X,Y]=meshgrid(d.x,d.y);
d.omega=(abs(X)<1)&(abs(Y)<1)&((X>0)|(Y>0));
r=0;
d.rr=zeros(N,N);
    for i=1:N
        for k=1:N
            if d.omega(i,k)
                r=r+1;
                d.ii(r)=i;
                d.kk(r)=k;
                d.rr(i,k)=r;
            end
        end
    end
d.NU=r;
end % setup domain



