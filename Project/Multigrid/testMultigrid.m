clear 
close all

%% TASK 1 
% task = 'CG';
%% TASK 2
% task = 'MGV';
%% TASK 3
 task = 'PCG';
% load('Reference_Sol.mat')
%% Initial Values
L = 6;
N = 2^L;
h = 1/N;
x = 0:h:1;
tol = 1e-12;
maxit = 1000;
nu1 = 3;
nu2 = 3;
level = 1;
e = 1.;
it = 0;
k = sqrt(8/h^2)+1;
%% Programs

f = @(x,y) 20*pi^2*sin(2*pi*x).*cos(4*pi*y); 
g = @(x,y) 0;

% Loading Geometry
[U rhs X Y] = geometryHelmholtz(N,f,g,task);

switch task
    case 'CG'
        [U,res,iter]=my_cg(U,rhs,N,tol,maxit,k);

        figure
        semilogy(res)
        tc = sprintf('Convergence CG for N: %d',N);
        title(tc)
        xlabel('Iterations')
        ylabel('Residual')
        grid on
    case 'MGV' 
        max_level = L;
        Au0 = matvec(U,N,h,k);

        Au0 = 1/h^2.*matvecMG(U,N,h,k);
       while e>tol  
        U  = mgv(U,rhs,N,nu1,nu2,level,max_level,k);             
        Au = 1/h^2.*matvec(U,N,h,k); 
        res(it + 1) = norm(rhs(2:N,2:N) - Au(2:N,2:N));
        e = norm(rhs(2:N,2:N) - Au(2:N,2:N))/norm(rhs(2:N,2:N) - Au0(2:N,2:N));
        it = it +1;
       end
       semilogy(res)
       tc = sprintf('Convergence Multigrid V-cycle for N: %d',N);
       title(tc)
       xlabel('Iterations')
       ylabel('Residual')
       grid on
%        
    case 'PCG'
        max_level = L;
        Au0 = matvec(U,N,h,k);
        title ('Initial guess')
        
       [U res it]  = precondCG(U,rhs,N,tol,maxit,nu1,nu2,X,Y,k);             

       semilogy(res)
       tc = sprintf('Convergence PCG for N: %d',N);
       title(tc)
       xlabel('Iterations')
       ylabel('Residual')
       grid on
end
   
