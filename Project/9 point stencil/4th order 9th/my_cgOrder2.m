function  [u,res,iter]=my_cgOrder2(task,u0,rhs,N,tol,maxit,k)
   %Initializations
    h = 1/N;
    b  = rhs;
    x0 = u0;
    x = x0;
    r0 = zeros(N+1);
    id = 2:N;
    %CG starts here
    Ax0 = 1/h^2.*matvec(task,x0,N,h,k);
    r0(id,id) = b(id,id) - Ax0(id,id);
    r = r0;
    p0 = r0;
    p = p0;
    error = 1; 
    iter = 0;
    while error>tol   || iter == maxit
        Ap0   = 1/h^2.*matvec(task,p0,N,h,k);
        alpha = sum(dot(r0,r0))/sum(dot(Ap0,p0));  %alpha=(r0,r0)/(Ap0,p0)
        x(id,id) = x0(id,id) + alpha*p0(id,id);    %x = x0 + alpha*p0
        r(id,id) = r0(id,id) - alpha*Ap0(id,id);   %r = r0 - alpha*Ap0
        beta = sum(dot(r,r))/sum(dot(r0,r0));      %beta = (r,r)/(r0,r0)
        p(id,id) = r(id,id) + beta*p0(id,id);      %p= r + beta*p0
        Ax = 1/h^2.*matvec(task,x,N,h,k);
        error = norm(b(id,id) - Ax(id,id))/norm(b(id,id) - Ax0(id,id));
        x0 = x;
        r0 = r;
        p0 = p;
        iter = iter + 1;
        res(iter) = norm(b(id,id) - Ax(id,id));
    end
    
    u = x;     
end

