function Ax = matvec(task,x,N,h,k)
    Ax = x; 
    index = 2:N;
    switch task
        case 'Order4'
            Ax(index,index) = (-20 + k^2 - k^4*h^2/12)*x(index,index)+x(index-1,index-1) + x(index+1,index+1)+... 
                                                   x(index+1,index-1) + x(index-1,index+1)+...
                                                    4*x(index-1,index) + 4*x(index+1,index)+...
                                                    4*x(index,index-1) + 4*x(index,index+1);
        case 'Order2'
            Ax(index,index) = ((4-k^2*h^2)*x(index,index)-x(index+1,index)-x(index,index+1)-...
                        x(index-1,index)-x(index,index-1));
    end
end
        