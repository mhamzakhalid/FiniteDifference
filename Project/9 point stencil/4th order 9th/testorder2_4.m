clear 
% close all
% load('Reference_Sol.mat')
load('Reference_SolV2.mat')
%% TASK  
 %  task = 'Order2'
  task = 'Order4'
%% Initial Values
 for L = 3:9
     tic
    N = 2^L;
    h = 1/N;
    x = 0:h:1;
    tol = 1e-12;
    maxit = 1000;
    nu1 = 3;
    nu2 = 3;
    level = 1;
    e = 1.;
    it = 0;
    k = sqrt(8/h^2)+1;
%% Programs
    f = @(x,y) 20*pi^2*sin(2*pi*x).*cos(4*pi*y); 
    g = @(x,y) 0;
    switch task 
        case {'Order4','Order2'}
            % Loading Geometry
            [U rhs X Y] = geometryHelmholtz(task,N,f,g,k);
            % Solving the Problem
            [U,res,iter] = my_cg(task,U,rhs,N,tol,maxit,k);
            
        case 'Preconditioners'
            
    end
    
    % Error Calculation
    c=find(ismember(x,x_ref));
    d=find(ismember(x_ref,x));
            % Calculating Error
    for pt = 1:length(c)
        e(pt)=U(c(pt),c(pt))-U_ref(d(pt),d(pt)); 
    end
    err(L-2) = norm(e)*sqrt(h);
    t_end(L-2) = toc;
 end

%% Convergence Plots
stp = 3:9;
dh = 1./(2.^stp);
p = polyfit(log10(dh),log10(err),1);
fprintf('Order of Congergence: %.2f \n',p(1))
%% Error
loglog(dh,err,'-*')
%% Time 
% loglog(dh,t_end,'-*')
xlabel('Step-size')
% ylabel('Error')
grid on
str = sprintf('Convergence in Space: %.2f',p(1));
title(str)
legend(task)
set(gca, 'FontName', 'Times New Roman')
%% Iteration Plots
% semilogy(1:iter,res,'-*')
% tc = sprintf('Convergence CG for N: %d',N);
% title(tc)
% xlabel('Iterations')
% ylabel('Residual')
% grid on

