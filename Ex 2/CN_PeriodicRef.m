function [x,t,U] = CN_PeriodicRef(M,N,T)
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

function Y = f(x)
    % Initial value.
    Y = sin(pi*x);    
end