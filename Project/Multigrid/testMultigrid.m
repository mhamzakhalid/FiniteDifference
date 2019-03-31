clear 
% close all

%% TASK 1 
 task = 'CG';
%% TASK 2
% task = 'MGV';
%% TASK 3
% task = 'PCG';
% load('Reference_Sol.mat')
%% Initial Values
% for L = 2:7;
L = 9;
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
switch task
    case 'CG'
%         f = @(x,y) 20*pi^2*sin(2*pi*x).*cos(4*pi*y); %CG
%         g = @(x,y) sin(2*pi*x).*cos(4*pi*y); %CG
          f = @(x,y) -20*pi^2*sin(2*pi*x).*cos(4*pi*y); 
          g = @(x,y) 0;
%           f = @(x,y) -1; 
%           g = @(x,y) 4*y.*(1-y);

    case {'MGV', 'PCG'}
          f = @(x,y) 20*pi^2*sin(2*pi*x).*cos(4*pi*y); 
          g = @(x,y) 0;

%                 f = @(x,y) 20*pi^2*sin(2*pi*x).*cos(4*pi*y); %CG
%         g = @(x,y) sin(2*pi*x).*cos(4*pi*y); %CG

%         f = @(x,y) -1; 
%         g = @(x,y) 4*y.*(1-y);
    otherwise 
        disp('Enter Correct case')
end
% Loading Geometry
[U rhs X Y] = geometryHelmholtz(N,f,g,task);

switch task
    case 'CG'
        [U,res,iter]=my_cg(U,rhs,N,tol,maxit,k);
%         mesh(X,Y,U)
%         ts = sprintf('Conjugate Gradient Solution for N: %d',N);
%         title(ts)
%         figure
%           semilogy(res)
%         tc = sprintf('Convergence CG for N: %d',N);
%         title(tc)
%         xlabel('Iterations')
%         ylabel('Residual')
%         grid on
    case 'MGV' 
        max_level = L;
        Au0 = matvec(U,N,h,k);

%         subplot(2,3,1) 
%         surf(X,Y,U)
%         title ('Initial guess')

        Au0 = 1/h^2.*matvec(U,N,h,k);
       while e>tol  
        U  = mgv(U,rhs,N,nu1,nu2,level,max_level,k);             
        Au = 1/h^2.*matvec(U,N,h,k); 
        res(it + 1) = norm(rhs(2:N,2:N) - Au(2:N,2:N));
        e = norm(rhs(2:N,2:N) - Au(2:N,2:N))/norm(rhs(2:N,2:N) - Au0(2:N,2:N));
        it = it +1;
%             if it<6
%                 subplot(2,3,it+1)
%                 surf(X,Y,U)
%                 t = sprintf('V-cycle Iteration:%d',it);
%                 title (t)
%             end
%          
       end
%        figure
%        surf(X,Y,U)
%        figure
       semilogy(res)
       tc = sprintf('Convergence Multigrid V-cycle for N: %d',N);
       title(tc)
       xlabel('Iterations')
       ylabel('Residual')
       grid on
%        
    case 'PCG'
        max_level = L;
        Au0 = matvec(U,N,h,k);

%         subplot(2,3,1) 
%         surf(X,Y,U)
        title ('Initial guess')
        
       [U res it]  = precondCG(U,rhs,N,tol,maxit,nu1,nu2,X,Y,k);             

%        figure
       semilogy(res)
       tc = sprintf('Convergence PCG for N: %d',N);
       title(tc)
       xlabel('Iterations')
       ylabel('Residual')
       grid on
end
   
%  c=find(ismember(x,x_ref));
%     d=find(ismember(x_ref,x));
%             % Calculating Error
%     for pt = 1:length(c)
%         e(pt)=U(c(pt),c(pt))-U_ref(d(pt),d(pt)); 
%     end
%     err(L-1) = norm(e)*sqrt(h)  ;
%  end
% 
% stp = 2:7;
% dh = 1./(2.^stp);
% loglog(dh,err)
% xlabel('Step-size')
% ylabel('Error')
% grid on
% title('Convergence in Space')
% set(gca, 'FontName', 'Times New Roman')
% 
% p = polyfit(log10(dh),log10(err),1);
% fprintf('Order of Congergence: %.2f \n',p(1))