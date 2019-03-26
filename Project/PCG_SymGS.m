function [U0, iter,err_symGS] = PCG_SymGS(A,F,U0,tol)

    r0 = F - A*U0;
    err_symGS(1)=norm(r0);
    z0=zeros(length(r0),1);
    %Applying Symmetric Gauss Seidal
    z0 = symGS(A,r0); 
    p0 = z0;
    error=1.;
    iter=0;
while error>tol 
    Ap0 = A*p0;
    alpha = dot(r0,z0)/dot(Ap0,p0);
    U1 = U0 + alpha*p0;
    r1 = r0 - alpha*Ap0;
    %Applying Symmetric Gauss Seidal iterations 

     z1=symGS(A,r1);   
    beta = dot(r1,z1)/dot(r0,z0);
    p1 = z1 + beta*p0;
    p0 = p1;
    r0 = r1;
    z0 = z1;
    U0 = U1;

    error = norm(F - A*U0);
    iter = iter +1;
    err_symGS(iter+1) = error;
   % err_est(iter) = 2*estimate^(1*iter)*sqrt(Upcg'*(M\A)*Upcg);

end
end