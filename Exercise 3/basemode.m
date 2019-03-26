function basemode(d)
[evec,eval]=eigs(-d.L,1,'sm');
s=sign(sum(evec));
field=zeros(size(d.omega));
    for r=1:d.NU
        field(d.ii(r),d.kk(r))=s*evec(r);
    end
mesh(field);
%axis off
print -depsc ml_logo_m.eps
figure
contour(field, 32);
axis equal
print -depsc ml_logo_c.eps
figure
imagesc(field);
axis equal
print -depsc ml_logo_i.eps
end