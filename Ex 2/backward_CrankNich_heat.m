clear
%% Case for plots Uncomment Solvers 
plotcase = 'time';
%plotcase = 'space';

switch plotcase
    case 'time'
        M = 20000;
       for stp=1:8
           N = 2^stp;
           T = .1;   
        % Solve the equation.
        
%%***************** UNCOMMENT TO RUN! *******************************        
         %[x,t,U] = backward_euler(M,N,T);
         [x,t,U] = crank_nich(M,N,T);
%%*******************************************************************
        k = T/2^stp;
        h = 1/M;
        r = k/h^2;
        e = time_conv(x,t,U);
        err(stp) = e(end-1)*sqrt(k);
       end
       stp = 1:8;
       dt = T./2.^stp;
       loglog(dt,err) 
       xlabel('Step-size')
       ylabel('Error')
       grid on
       title('Convergence in Time')
       set(gca, 'FontName', 'Times New Roman')

    case 'space'
        N = 10000;
       for stp=1:8
           M = 2^stp;
           T = .1;
           %Solve the equation
%%***************** UNCOMMENT TO RUN! *******************************        
         %[x,t,U] = backward_euler(M,N,T);
         [x,t,U] = crank_nich(M,N,T);
%%*******************************************************************
            k = T/N;
            h = 1/M;
            r = k/h^2;
            e = space_conv(x,t,U);
            err(stp) = norm(e)*sqrt(h);
       end
       stp = 1:8;
       H = 1./2.^stp;
       loglog(H,err)
       xlabel('Step-size')
       ylabel('Error')
       grid on
       title('Convergence in Space')
       set(gca, 'FontName', 'Times New Roman')
end
    
%% Error Calculation
function error = time_conv(x,t,U)
    exact = @(t,x) sin(pi.*x).*exp(-pi^2.*t);     
     for i=1:length(t)
        error(i) = norm(U(:,i) - exact(t(i),x)');
     end
end
function error = space_conv(x,t,U)
    exact = @(t,x) sin(pi.*x).*exp(-pi^2.*t);    
    error = U(:,end-1) - exact(t(end-1),x)';
end
%% Plots
function plot_solution(x,t,U)
    % Plotting the solution of the heat equation.
    figure
    [T,X] = meshgrid(t,x);
    surf(T,X,U);
    shading interp
    colormap cool;
    xlabel('t')
    ylabel('x')
    title('Solution')
    saveas(gca,'Solution_2.png')
end

%% Numerical Schemes
%------------------ Backward Euler
function [x,t,U] = backward_euler(M,N,T)
    h = 1/M;    % Step size in space.
    k = T/N;    % Step size in time.
    r = k/h^2;
    
    % Initializing arrays.
    U = zeros(M+1,N+1);     % Array to store the solution, boundaries included.
    x = linspace(0,1,M+1);  % Gridpoints on the x-axis.
    t = linspace(0,T,N+1);  % Gridpoints on the t-axis.
    U(:,1) = f(x);          % Initial values.
    e = ones(M-1,1);
    A = spdiags([-r*e 1+2*r*e -r*e], -1:1, M-1, M-1);
    F = f(x)';
    U(:,1) = F;
   for i = 1:N
    U(2:M,i+1) = A\U(2:M,i);
   end
end

%****************** Crank Nicholson
function [x,t,U] = crank_nich(M,N,T)
    h = 1/M;    % Step size in space.
    k = T/N;    % Step size in time.
    r = .5*k/h^2;
    
    % Initializing arrays.
    U = zeros(M+1,N+1);     % Array to store the solution, boundaries included.
    x = linspace(0,1,M+1);  % Gridpoints on the x-axis.
    t = linspace(0,T,N+1);  % Gridpoints on the t-axis.
    U(:,1) = f(x);          % Initial values.
    e = ones(M-1,1);
    A1 = spdiags([-r*e 1+2*r*e -r*e], -1:1, M-1, M-1);
    A2 = spdiags([r*e 1-2*r*e r*e], -1:1, M-1, M-1);
    F = f(x)';
    U(:,1) = F;
   for i = 1:N
      rhs = A2*U(2:M,i);
    U(2:M,i+1) = A1\rhs;
    clear rhs
   end
end
%% Function
function Y = f(x)
    % Initial value.
    Y = sin(pi*x);    
end