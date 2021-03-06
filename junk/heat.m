function heat
 plotcase = 'time';
%plotcase = 'space';
 plotcase = 'Fishers Problem';
switch plotcase
    case 'time'
    M = 20000;
   for stp=1:8
       N = 2^stp;
       T = .1;
    
    % Solve the equation.
     %[x,t,U] = forward_euler(M,N,T);
     %[x,t,U] = backward_euler(M,N,T);
      %[x,t,U] = crank_nich(M,N,T);
     %[x,t,U] = CN_Periodic(20,40,.5);
      [x,t,U] = diffprob(20,400,T);

%      for i=1:length(U)
%          plot(x,U(:,i))
%          pause(0.05)
%      end
    k = T/2^stp;
    h = 1/M;
    r = k/h^2;
    e = time_conv(x,t,U);
    err(stp) = e(end-1)*sqrt(k);
   end
   stp = 1:8;
   dt = T./2.^stp;
   loglog(dt,err)
   % Creating 3D-plot.
    %plot_solution(x,t,U)
    
    case 'space'
        N = 1000;
   for stp=1:8
       M = 2^stp;
       T = .1;
    
    % Solve the equation.
     %[x,t,U] = forward_euler(20,400,0.5);
     %[x,t,U] = backward_euler(M,N,T);
      [x,t,U] = crank_nich(M,N,T);
     %[x,t,U] = CN_Periodic(20,40,.5);
     %[x,t,U] = diffprob(20,400,T);

%      for i=1:length(U)
%          plot(x,U(:,i))
%          pause(0.05) 
%      end
    k = T/N;
    h = 1/M;
    r = k/h^2;
    e = space_conv(x,t,U);
    err(stp) = norm(e)*sqrt(h);
   end
   stp = 1:8;
   H = 1./2.^stp;
   loglog(H,err)
   % Creating 3D-plot.
    %plot_solution(x,t,U)
    case 'Fishers Problem'
     [x,t,U] = diffprob(20,400,10);
      plot_solution(x,t,U)

end
    
end

function error = time_conv(x,t,U)
    exact = @(t,x) sin(pi.*x).*exp(-pi^2.*t);     
     for i=1:length(t)
        error(i) = norm(U(:,i) - exact(t(i),x)');
     end
end
function error = space_conv(x,t,U)
    exact = @(t,x) sin(pi.*x).*exp(-pi^2.*t);    
    error = U(:,end-1) - exact(t(end-1),x)';
%      for i=1:length(x)
%         error(i) = norm(U(i,end-1) - exact(t(end-1),x(i))');
%      end
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
    F = sin(pi.*(x-.25)).^100';
    U(:,1) = F;
    [t,y] = ode23s(@(t,y) ode_sys(t,y),[0 T],U(:,1));
    
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
        dydt = 0.01*A*y +y.*(1-y);
        %y.*(ones(length(y),1) - y)
end
function Y = f(x)
    % Initial value.
    Y = sin(pi*x);    
end