clear
clc
%% Cases
 %   test = 'Poisson1D 2(f)1';
%   test = 'Poisson1D 2(f)2';
%   test = 'Poisson1D 2(f)3';

%% Initial Values
M_h = 14;
for M = 2:M_h
    h = 1/2^M;
    u = zeros(2^M-1,1);
    e = ones(2^M-1,1);
    A = spdiags([-e./h^2 2*e/h^2 -e/h^2],-1:1,2^M-1,2^M-1);
    x = 0:h:1;

    %% Boundary and exact solution
    switch test
        case 'Poisson1D 2(f)1'
            f = @(x) sin(pi.*x);
            exact = @(x) sin(pi.*x)./pi^2;
            rhs = f(x);
            u_exact = exact(x); 
            boundary = zeros(2^M-1,1);
        case 'Poisson1D 2(f)2'
            f = @(x) 1./x;
            exact = @(x) -x.*log(x)+ x;
            rhs = f(x);
            u_exact = exact(x); 
            alpha = 0.;
            beta  = 1.;
            boundary = zeros(2^M-1,1);
            boundary(1) = 0.;
            boundary(end) = 1./h^2;
        case 'Poisson1D 3'
            f = @(x) sin(pi.*x);
            exact = @(x) -x.*log(x)+ x;
            rhs = f(x);
            u_exact = exact(x); 
            alpha = 0.;
            beta  = 1.;
            boundary = zeros(2^M-1,1);
            boundary(1) = 0.;
            boundary(end) = 1./h^2;
    end

    %% Solving the problem
    switch test
        case 'Poisson1D 2(f)1'
           u = pcg(A,rhs(2:2^M)'+boundary);
           error(M) = norm(u - u_exact(2:2^M)');
        case 'Poisson1D 2(f)2'
            F = rhs(2:2^M)'+ boundary;
            u = A\F;
           %u = pcg(A,rhs(2:2^M)'+boundary);
           error(M) = norm(u - u_exact(2:2^M)'); 
    end
end

%% Visualization
m = 1:M_h;
H = 1./2.^m;
loglog(H,error)
hold on
axis('equal')
grid on
xlabel('Step-size')
ylabel('Error')
set(gca, 'FontName', 'Times New Roman')