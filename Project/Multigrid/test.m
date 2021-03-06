clear


%Initial values

grids =[10,20,50];
 for gr = 1:length(grids)
    N = grids(gr);
%for L=5:9   
    
 %   N = 2^L;
    h = 1/N;
    U = zeros(N+1);
    tol =  1e-12;
    maxit = 100;
    %grid
    x = 0:h:1;
    y = 0:h:1;
    [X,Y] = meshgrid(x,y);

    %% CG
    %Defining Right Hand side and boundary
    f = @(x,y) 20*pi^2*sin(2*pi*x).*cos(4*pi*y); %CG
    g = @(x,y) sin(2*pi*x).*cos(4*pi*y); %CG
    %Boundary on U
    U(:,1) = g(X(:,1),Y(:,1));
    U(:,end) = g(X(:,end),Y(:,end));
    U(1,:) = g(X(1,:)',Y(1,:)');
    U(end,:) = g(X(end,:)',Y(end,:)');
    rnd = rand(N-1);
    U(2:N,2:N) = rnd(1:end,1:end);
    %Right Hand side
    rhs = f(X,Y);
    %Implementing cg
    [U,res,iter]=my_cg(U,rhs,N,tol,maxit);
    %surf(X,Y,U)
    semilogy(res)
    hold on
    
    
 end
    grid on
    legend('N = 10','N = 20','N = 50')
    

%     %% MGV
%     f = @(x,y) -1; 
%     g = @(x,y) 4*y.*(1-y);
%     %Boundary Conditions
%     U(:,1) = g(X(:,1),Y(:,1));
%     rnd = rand(N-1);
%     U(2:N,2:N) = rnd(1:end,1:end);
%     %%Right Hand side
%     rhs = zeros(N+1);
%     rhs(:,:) = f(X,Y);
% 
%     %% MGV
%     nu1 = 3;
%     nu2 = 3;
%     max_level = L;
%     level = 1;
%     e =1.;
%     it = 0;
%     %Au0 = matvec(U,N);
% %     [U,res,iter]=my_cg(U,rhs,N,tol,maxit);
%     %surf(X,Y,U)
%     %[U res] = jacobi(U,rhs,2/3,N,5);
% 
% %       subplot(2,3,1) 
% %       surf(X,Y,U)
% %       title ('Initial guess')
% %  
% %     Au0 = 1/h^2.*matvec(U,N);
% %      while e>tol
% %     %  for cycle = 1:5    
% %        U  = mgv(U,rhs,N,nu1,nu2,level,max_level);     
% %        %subplot(2,3,cycle+1)
% %        Au = 1/h^2.*matvec(U,N); 
% %        %res(cycle) = norm(rhs(2:N,2:N) - Au(2:N,2:N));
% %        e = norm(rhs(2:N,2:N) - Au(2:N,2:N))/norm(rhs(2:N,2:N) - Au0(2:N,2:N));
% %        er(it+1) = norm(rhs(2:N,2:N) - Au(2:N,2:N));
% %        it = it +1;
% % %        surf(X,Y,U)
% % %        t = sprintf('V-cycle Iteration:%d',cycle);
% % %        title (t)
% %       %[U res it] = precondCG(U,rhs,N,tol,maxit,nu1,nu2);
% %      % end
% %      end
% %      fprintf('Iteration:%d\n nu1:%d\n nu2:%d\n ',it,nu1,nu2);
%       %figure
% %     semilogy(er)
%      title('Conjugate Gradient Solution')
%      grid on
%      hold all
%      %xlabel('Iteration')
%      %ylabel('Residual')
%      t = sprintf('N = 2^%d',L);
%      legend(t)
%      get(legend(gca),'String'); 
%     
%      
% end
%     %figure
%    
% 
% 
% 
