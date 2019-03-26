clear

    % scheme = 'FTBS'
      scheme = 'Lax-Wendrof'
     %scheme = 'Wendrof'
%% Problem 1
% r = .5;
% h = 1/160;
% k = r*h;
% x = 0:h:2.2;
% t = 0:k:2;
% %%Boundary
% g = 1;
% f = 0 ;
%% Problem 2 
k = 1/400;
h = 1/200;
r = k/h;
x = 0:h:2.2;
t = 0:k:2;
%Boundary
g = 1;
f = exp(-64*(x-.5).^2).*sin(32*pi*x); 
%% Exact 
[X,T]=meshgrid(x,t);
exact = X<T;
exact = exact';
%% System and BC setup
nx = length(x);
nt = length(t);
u = zeros(nx,nt);
u(1,1:end)=g;
u(:,1) = f';
%% FTBS

switch scheme
    case 'FTBS'
        e = ones(nx-1,1);
        A = spdiags([r*e (1-r)*e 0*e], -1:1, nx-1, nx-1);
        A(end,end-1) = 2;
        A(end,end-2) = -1;
        bound = zeros(nx-1,1);
        
        for i=1:nt-1
           bound(1) = r*u(1,i);
           u(2:end,i+1) = A*u(2:end,i) + bound;
        end
        
    case 'Lax-Wendrof'
        
        e = ones(nx-2,1);
        A = spdiags([(r/2+r^2/2)*e (1-r^2)*e (-r/2 + r^2/2)*e], -1:1, nx-2, nx-2);
        bound = zeros(nx-2,1);       
        for i=1:nt-1 
            u(end,i) = 2*u(end-1,i) - u(end-2,i);
            bound(1) = (r/2+r^2/2)*u(1,i);
            u(2:end-1,i+1) = A*u(2:end-1,i) + bound;
        end
        
    case 'Wendrof'
        AW = (1-r)/(1+r);
        e = ones(nx-2,1);
        A = spdiags([e AW*e 0*e], -1:1, nx-2, nx-2);
        M = spdiags([AW*e e 0*e], -1:1, nx-2, nx-2);
        bound = zeros(nx-2,1);
        
        for i=1:nt-1
           bound(1) = u(1,i) - AW*u(1,i+1);
           u(2:end-1,i+1) = M\(A*u(2:end-1,i) + bound);
           u(end,i+1) = 2*u(end-1,i+1) -u(end-1,i+1);
        end
        
end

for i=1:30:nt
plot(x,u(:,i))
hold on
pause(0.1)
end

        