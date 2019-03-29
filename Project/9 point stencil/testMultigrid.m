clear 

load('Reference_Sol.mat')
% load('Reference_Solinhom.mat')
%% TASK 1 
 task = 'CG';

%% Initial Values
for L = 2:9
    N = 2^L;
    h = 1/N;
    x = 0:h:1;
    tol = 1e-6;
    maxit = 1000;
    nu1 = 3;
    nu2 = 3;
    level = 1;
    e = 1.;
    it = 0;
    k = 1e2;
%% Programs
switch task
    case 'CG'
        f = @(x,y) 20*pi^2*sin(2*pi*x).*cos(4*pi*y); 
%       f = @(x,y) 0; 
          g = @(x,y) 0;
end
% Loading Geometry
[U rhs X Y] = geometryHelmholtz(N,f,g,task);

switch task
    case 'CG'
        [U,res,iter]=my_cg(U,rhs,N,tol,maxit,k);

%           semilogy(1:iter,res)
%         tc = sprintf('Convergence CG for N: %d',N);
%         title(tc)
%         xlabel('Iterations')
%         ylabel('Residual')
%         grid on

end
    c=find(ismember(x,x_ref));
    d=find(ismember(x_ref,x));
            % Calculating Error
    for pt = 1:length(c)
        e(pt)=U(c(pt),c(pt))-U_ref(d(pt),d(pt)); 
    end
    err(L-1) = norm(e)*sqrt(h);
 end

stp = 2:9;
dh = 1./(2.^stp);
loglog(dh,err)
xlabel('Step-size')
ylabel('Error')
grid on
title('Convergence in Space')
set(gca, 'FontName', 'Times New Roman')
p = polyfit(log10(dh),log10(err),1);
fprintf('Order of Congergence: %.2f \n',p(1))

