clear
%% Cases

plotcase = 'time';
%plotcase = 'space';

%% Reference Solution
%[xref,tref,Uref] = CN_PeriodicRef(1000,1000^2,.1);

%% Convergence Plots
switch plotcase
    case 'time'
    load('RefPeriodic.mat')

    M = 1000;
    for stp=1:8
           N = 2^stp;
           T = .1;
           k = T/N;
        % Solve the equation.
        [x,t,U] = CN_Periodic(M,N,T);
        % Calculating Error
        e = norm(Uendref(:,end) - U(:,end));
        err(stp) = e*sqrt(k);
    end
   % Convergence Plot
   stp = 1:8;
   dt = T./2.^stp;
   loglog(dt,err)
   xlabel('Time-Step')
   ylabel('Error')
   grid on
   title('Convergence in Time')
   set(gca, 'FontName', 'Times New Roman')
    case 'space'
        load('RefPeriodic.mat')

        N = 100;
        for stp=1:8
               M = 2^stp-1;
               T = .1;
               h = 1/M;
            % Solve the equation.
            [x,t,U] = CN_Periodic(M,N,T);
            % Finding same vertices
            c=find(ismember(x,xref));
            d=find(ismember(xref,x));
            % Calculating Error
            e=U(c,end)-Uendref(d,end);       
            err(stp) = norm(e)*sqrt(h);
        end
   stp = 1:8;
   dt = 1./(2.^stp-1);
   loglog(dt,err)
   xlabel('Step-size')
   ylabel('Error')
   grid on
   title('Convergence in Space')
   set(gca, 'FontName', 'Times New Roman')

end

%% ********** Periodic Boundary Conditions Problem Solution

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

%% Function
function Y = f(x)
    % Initial value.
    Y = sin(pi*x);    
end