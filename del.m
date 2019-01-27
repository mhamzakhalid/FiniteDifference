clear

%% Initial Values
M_h = 14;
for M = 2:M_h
    h=3.0518e-05

            u = zeros(2^M,1);
            e = ones(2^M,1);
            A = spdiags([-e 2*e -e],-1:1,2^M,2^M);
            A(end,end-1) = -h;
            A(end,end) = (h+h^2);
            A = 1/h^2.*A;
            x = 0:h:1;

            %% Boundary and exact solution
            f = @(x) sin(pi.*x);            
            rhs = f(x);
            rhs(end) =  0.;
            exact = @(x) 1./pi^2.*sin(pi.*x) + (pi+1)./(2.*pi).*x;
            u_exact = exact(x);
            boundary = zeros(2^M,1);
            boundary(end) = 1.;

 
    %% Solving the problem
            F  = rhs(2:end)'+ boundary;
            u = A\F;
            error(M) = norm(u - u_exact(2:2^M+1)');

    
end

%% Visualization
m = 1:M_h;
H = 1./2.^m;
loglog(H,error,'-*')

hold on
axis('equal')
grid on
xlabel('Step-size')
ylabel('Error')
set(gca, 'FontName', 'Times New Roman')
legend('case 1','case 2','case 3')