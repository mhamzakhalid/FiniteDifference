clear
     [x,t,U] = diffprob(20,40000,10);
      plot_solution(x,t,U)    

% ***************** Diffusion Problem
function [x,t,U] = diffprob(M,N,T)
     h = 1/M;    % Step size in space.
     k = T/N;    % Step size in time.
    % Initializing arrays.
    U = zeros(M+1,1);     % Array to store the solution, boundaries included.
    x = linspace(0,1,M+1);  % Gridpoints on the x-axis.
    t = linspace(0,T,N+1);  % Gridpoints on the t-axis.
    F = sin(pi.*(x-.25)).^100';
    U(:,1) = F;
    [t,y] = ode15s(@(t,y) ode_sys(t,y),[0 T],U(:,1));
    U = y';  
    
end

function dydt = ode_sys(t,y)
        h = 1/(length(y)-1);
        e = ones(length(y),1);
        A = spdiags([e -2*e e], -1:1, length(y), length(y));
        A(1,1) = -2;
        A(1,2) = 2;
        A(end,end) = -2;
        A(end,end-1) = 2;
        r = 0.01/h^2;
        A = r.*A;
        dydt = A*y + y.*(1-y);
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
    %saveas(gca,'Solution_2.png')
end