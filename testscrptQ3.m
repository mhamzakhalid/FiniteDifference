clear

%% Cases
   % test = 'Poisson1D 3_1';
    % test = 'Poisson1D 3_2';
    test = 'Poisson1D 3_3';
%% Initial Values
M_h = 10;
for M = 2:M_h
    h0 = .1;
    h = h0/2^M;
    switch test
        case 'Poisson1D 3_1'
            x = 0:h:1;
            u = zeros(length(x)-1,1);
            e = ones(length(x)-1,1);
            A = 1/h^2.*spdiags([-e 2*e -e],-1:1,length(x)-1,length(x)-1);
            A(end,end-1) = -1/h;
            A(end,end) = (1+h)/h;
            
            %% Boundary and exact solution
            f = @(x) sin(pi.*x);            
            rhs = f(x);
            rhs(end) =  1.;
            exact = @(x) 1./pi^2.*sin(pi.*x) + (pi+1)./(2.*pi).*x;
            u_exact = exact(x);
            boundary = zeros(length(x)-1,1);

        case 'Poisson1D 3_2'
            x = 0:h:1;
            u = zeros(length(x)-1,1);
            e = ones(length(x)-1,1);
            A = spdiags([-e./h^2 2*e/h^2 -e/h^2],-1:1,length(x)-1,length(x)-1);
            A(end,end-2) = 1/(2*h);
            A(end,end-1) = -2/h;
            A(end,end) = 3/(2*h)+ 1;
            
            %% Boundary and exact solution
            f = @(x) sin(pi.*x);            
            rhs = f(x);
            rhs(end) = 0.;
            exact = @(x) 1./pi^2.*sin(pi.*x) + (pi+1)./(2.*pi).*x;
            u_exact = exact(x);
            boundary = zeros(length(x)-1,1);
            boundary(end) = 1.;
        case 'Poisson1D 3_3'
            x = 0:h:1;
            u = zeros(length(x)-1,1);
            e = ones(length(x)-1,1);
            A = spdiags([-e./h^2 2*e/h^2 -e/h^2],-1:1,length(x)-1,length(x)-1);
            A(end,end-1)=-2/h^2;
            A(end,end)=(2+2*h)/h^2;
            
            %% Boundary and exact solution
            f = @(x) sin(pi.*x);            
            rhs = f(x);
            rhs(end) =h^2*rhs(end);
            exact = @(x) 1./pi^2.*sin(pi.*x) + (pi+1)./(2.*pi).*x;
            u_exact = exact(x);
            boundary = zeros(length(x)-1,1);
            boundary(end) = 2/h;
    end
    %% Solving the problem
            F  = rhs(2:end)'+boundary;
            u = A\F;
            error(M) = norm(u - u_exact(2:length(x))');
    
end
%% Visualization
m = 1:M_h;
H = 1./2.^m;
switch test
    case 'Poisson1D 3_1'
        loglog(H,error,'-*')
    case 'Poisson1D 3_2'
        loglog(H,error,'-o')
    case 'Poisson1D 3_3'
        loglog(H,error,'-+')
end
hold on
axis('equal')
grid on
xlabel('Step-size')
ylabel('Error')
set(gca, 'FontName', 'Times New Roman')
legend('Case 1','Case 2','Case 3')