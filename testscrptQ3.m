clear

%% Cases
    test = 'Poisson1D 3_1';
   %  test = 'Poisson1D 3_2';
   % test = 'Poisson1D 3_3';
load 'ref.mat'
%% Initial Values
M_h = 14;
% M = 14;
for M = 2:M_h
    h = 1/2^M;
    switch test
        case 'Poisson1D 3_1'
            u = zeros(2^M,1);
            e = ones(2^M,1);
            A = spdiags([-e./h^2 2*e/h^2 -e/h^2],-1:1,2^M,2^M);
            A(end,end-1)=-1/h;
            A(end,end)=1+1/h;
            x = 0:h:1;

            %% Boundary and exact solution
            f = @(x) sin(pi.*x);            
            rhs = f(x);
            rhs(end) = 1.;
            exact = @(x) 1./pi^2.*sin(pi.*x) + (pi+1)./(2.*pi).*x;
            u_exact = exact(x);
            boundary = zeros(2^M,1);

        case 'Poisson1D 3_2'
            u = zeros(2^M,1);
            e = ones(2^M,1);
            A = spdiags([-e./h^2 2*e/h^2 -e/h^2],-1:1,2^M,2^M);
            A(end,end-2)=1/(2*h);
            A(end,end-1)=-2/h;
            A(end,end)=3/(2*h)+ 1;
            x = 0:h:1;

            %% Boundary and exact solution
            f = @(x) sin(pi.*x);            
            rhs = f(x);
            rhs(end) = 0.;
            exact = @(x) 1./pi^2.*sin(pi.*x) + (pi+1)./(2.*pi).*x;
            u_exact = exact(x);
            boundary = zeros(2^M,1);
            boundary(end) = 1.;
        case 'Poisson1D 3_3'
            u = zeros(2^M,1);
            e = ones(2^M,1);
            A = spdiags([-e./h^2 2*e/h^2 -e/h^2],-1:1,2^M,2^M);
            A(end,end-1)=-2/h^2;
            A(end,end)=(2+2*h)/h^2;
            x = 0:h:1;

            %% Boundary and exact solution
            f = @(x) sin(pi.*x);            
            rhs = f(x);
            rhs(end) =h^2*rhs(end);
            exact = @(x) 1./pi^2.*sin(pi.*x) + (pi+1)./(2.*pi).*x;
            u_exact = exact(x);
            boundary = zeros(2^M,1);
            boundary(end) = 2/h;
    end
    %% Solving the problem
            F  = rhs(2:end)'+boundary;
            u = A\F;
            error(M) = norm(u - u_exact(2:2^M+1)')
%             F  = rhs(2:end)'+boundary;
%             %u = pcg(A,F);
%             %F  = rhs(2:end)'+boundary; 
%             u = A\F;
%             error(M) = norm(u - u_exact(2:2^M+1)')
    
end
% save('ref.mat','u_ex')

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
