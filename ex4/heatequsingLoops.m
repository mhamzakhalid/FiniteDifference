clear

  %  scheme = 'FTBS'
     scheme = 'Lax-Wendrof'
  %   scheme = 'Wendrof'
%% Problem 1
r = .5;
h = 1/160;
k = r*h;
x = 0:h:2.2;
t = 0:k:2;
%%Boundary
g = 1;
f = 0 ;
%%Exact 
[X,T]=meshgrid(x,t);
exact = X<T;
exact = exact';
%% Problem 2 
% k = 1/400;
% h = 1/200;
% r = k/h;
% x = 0:h:2.2;
% t = 0:k:2;
% %%Boundary
% g = 1;
% f = exp(-64*(x-.5).^2).*sin(32*pi*x); 

%% System and BC setup
nx = length(x);
nt = length(t);
u = zeros(nx,nt);
u(1,2:end)=g;
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
        
        for j=1:nt-1
            for i=2:nx-1
                u(i,j+1) = u(i,j) - r/2*(u(i+1,j)-u(i-1,j))+...
                            r^2/2*(u(i+1,j)-2*u(i,j)+u(i-1,j));
            end
        end
        
    case 'Wendrof'
        for j=1:nt-1
            for i=2:nx-1
                u(i,j+1) = u(i-1,j) - (1-r)/(1+r)*(u(i-1,j+1) - u(i,j));                
            end
        end
        
end
% norm(u-exact)
for i=1:30:nt
plot(u(:,i))
hold on
pause(0.1)
end


        