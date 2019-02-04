
function heat
    M = 20;
    N = 40;
    T = 0.5;
    %global M N T
    % Solve the equation.
    %[x,t,U] = forward_euler(20,400,0.5);
     %[x,t,U] = backward_euler(20,400,0.5);
     %[x,t,U] = crank_nich(20,400,0.5);
     %[x,t,U] = CN_Periodic(20,40,.5);
     [x,t,U] = diffprob(20,400,.5);
     tspan = linspace(0,T,size(U,2));
     for i=1:length(U)
         plot(x,U(:,i))
         pause(0.05)
     end
     % Plot the solution at some points in time.
    figure
    hold on
    %tplots = [0.0,0.1,0.2,0.3,0.4,0.5];
    %for i = 1:6
    %    tn = tplots(i);
    %    plot(x,U(:,i+1),'-o','DisplayName',strcat('t=',num2str(tn)))
    %end
    %xlabel('x')
    %ylabel('u(x,t)')
    %legend()
    %saveas(gca,'Solution_1.png')
    
    % Creating 3D-plot.
    plot_solution(x,t,U)
end



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
% ************* Periodic Boundary Conditions

function [x,t,U] = CN_Periodic(M,N,T)
    h = 1/M;    % Step size in space.
    k = T/N;    % Step size in time.
    r = .5*k/h^2;
    
    % Initializing arrays.
    U = zeros(M+1,N+1);     % Array to store the solution, boundaries included.
    x = linspace(0,1,M+1);  % Gridpoints on the x-axis.
    t = linspace(0,T,N+1);  % Gridpoints on the t-axis.
    U(:,1) = f(x);          % Initial values.
    e = ones(M+1,1);
    A1 = spdiags([-r*e 1+2*r*e -r*e], -1:1, M+1, M+1);
    A2 = spdiags([r*e 1-2*r*e r*e], -1:1, M+1, M+1);
    % imposing periodic BC
    A1(1,end) = -r;
    A1(end,1) = -r;
    A2(1,end) = r;
    A2(end,1) = r;
    F = f(x)';
    U(:,1) = F;
   for i = 1:N
      rhs = A2*U(:,i);
    U(:,i+1) = A1\rhs;
    clear rhs
   end
end

% ***************** Diffusion Problem
function [x,t,U] = diffprob(M,N,T)
     h = 1/M;    % Step size in space.
    k = T/N;    % Step size in time.
      
    % Initializing arrays.
    U = zeros(M+1,1);     % Array to store the solution, boundaries included.
    x = linspace(0,1,M+1);  % Gridpoints on the x-axis.
    t = linspace(0,T,N+1);  % Gridpoints on the t-axis.
    U(:,1) = f(x);          % Initial valuesS
    F = sin(pi.*(x - 0.25)).^100';
    U(:,1) = F;
    [t,y] = ode23s(@(t,y) ode_sys(t,y),[0 10],U(:,1));
    
    U = y';
    
end
% ************************* Forward Euler
function [x,t,U] = forward_euler(M,N,T)
    % Solving u_t=u_xx on [0,1] with Dirichlet conditions u(0,t)=u(1,t)=0
    % and initial value u(x,0)=f(x) over the time interval [0,T].
    % Input parameters:
    %       M,N: number of grid intervals in the x- and t directions.
    %       T  : end of integration.
    % Input parameters:
    %       x,t: the gridpoints in the x- and t- directions.
    %       U  : An array with the numerical solution.
    
    % Setting step sizes
    h = 1/M;    % Step size in space.
    k = T/N;    % Step size in time.
    r = k/h^2;
    
    % Initializing arrays.
    U = zeros(M+1,N+1);     % Array to store the solution, boundaries included.
    x = linspace(0,1,M+1);  % Gridpoints on the x-axis.
    t = linspace(0,T,N+1);  % Gridpoints on the t-axis.
    U(:,1) = f(x);          % Initial values.
    
    % Main loop.
    for i = 1:N
        U(2:M,i+1) = U(2:M,i)+r*(U(3:M+1,i)-2*U(2:M,i)+U(1:M-1,i));
    end
end

function dydt = ode_sys(t,y)
        h = 1/(length(y)-1);
        e = ones(length(y),1);
        A = spdiags([e./h^2 -2*e./h^2 e./h^2], -1:1, length(y), length(y));
        A(1,1) = 3/(2*h);
        A(1,2) = -2/h;
        A(1,3) = 1/(2*h);
        A(end,end) = 3/(2*h);
        A(end,end-1) = -2/h;
        A(end,end-2) = 1/(2*h);
        dydt = 0.01*A*y + y.*(ones(length(y),1) - y);
end
function Y = f(x)
    % Initial value.
    Y = sin(pi*x);    
end