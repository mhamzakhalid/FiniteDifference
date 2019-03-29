clear

N = 64;
n = N+1;
L = 1 ;
dy = L/(N-1);
dx = L/(N-1); %step-length                                                   
k = sqrt(8/dx^2)+1; %Parameter in Helmotz Equation
%k = 1e2;
x = linspace(0,L,n);
y = linspace(0,L,n);
U = zeros(n);
for i=1:n
    for j=1:n
       %    F(i,j)= -exp(-50*((x(j) - 1/2)^2))*exp(-50*(y(i) - 1/2)^2);
        F(i,j) =   -20*pi^2*sin(2*pi*x(j)).*cos(4*pi*y(i));
    end 
end 
F = dx^2*F;
for i = 1:N-1
    rhs(i,:) = F(i+1,2:N);
end
rhs = reshape(rhs,[(N-1)*(N-1),1]);  % Right Hand side
%% Matrix setup
%L-Matrix
e = ones((N-1)*(N-1),1);
f = ones((N-1)*(N-1),1);
f((N-1):(N-1):(N-1)*(N-1)) = 0;
g = [1;f(1:end-1)];
Lap = spdiags([e,f, -4*e, g, e],[-(N-1),-1,0,1,(N-1)],(N-1)*(N-1),(N-1)*(N-1));
% Setting up A-Matrix 
A = Lap + dx^2*k^2*eye((N-1)*(N-1));
A = sparse(A);
if eig(A)>0
    %Postive Definite Case
    disp('A is Positive Definite Matrix')
else
    %A not Positive Definite 
    disp('Not PDM')
end
%CG test
u0 = zeros((N-1)^2,1);
tol = 1e-8;
%u= A \rhs;
%[u, iter,err_symGS] = PCG_SymGS(A,rhs,u0,tol);
[u,iter,error] = conj_grad(u0,A,rhs,tol);
u = reshape(u,N-1,N-1)';
 %semilogy(error)
for i = 1:N-1
    U(i+1,2:N) = u(i,:);
end
mesh(x,y,U)
% figure
% contour(U,100)

% for i=1:101
% plot(x,U(i,:))
% pause(0.1)
% end
%%  Functions
function  [x,i,error]= conj_grad(x0,A,b,tol)
    r0 = b - A*x0;
    p0 = r0;
    error = 1; 
    i = 0;
    while error>tol   
        Ap0 = A*p0;
        alpha = dot(r0,r0)/dot(Ap0,p0);
        x = x0 + alpha*p0;
        r = r0 - alpha*Ap0;
        beta = dot(r,r)/dot(r0,r0);
        p = r + beta*p0;  
        i = i + 1;  
        error(i) = norm(b - A*x);
        x0 = x;
        r0 = r;
        p0 = p;
        
    end
end
% Symmetric Gauss-Seidal

