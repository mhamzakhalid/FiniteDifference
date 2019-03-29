function Ax = matvec(x,N,h,k)
    Ax = x; 
    index = 2:N;
    Ax(index,index) = (4-k^2*h^2)*x(index,index)-x(index+1,index)-x(index,index+1)-...
                        x(index-1,index)-x(index,index-1);
end
        