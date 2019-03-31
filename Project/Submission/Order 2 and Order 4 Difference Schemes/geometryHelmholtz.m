function [U rhs X Y] = geometryHelmholtz(task,N,f,g,k)
    
    switch task
        case 'Order4'
            h = 1/N;
            U = zeros(N+1);
            %grid
            x = 0:h:1;
            y = 0:h:1;
            [X,Y] = meshgrid(x,y);
            U(:,1) = g(X(:,1),Y(:,1));            
            U(:,end) = g(X(:,end),Y(:,end));
            U(1,:) = g(X(1,:)',Y(1,:)');
            U(end,:) = g(X(end,:)',Y(end,:)');          
            rnd = rand(N-1);
            U(2:N,2:N) = rnd(1:end,1:end);
            %Right Hand side
            rhs = zeros(N+1);
            rhs = f(X,Y);
            index = 2:N;
            rhs(index,index) = (h^2/12 - 3 - h^2*k^2/12)*rhs(index,index) +...
                                h^2/12*(rhs(index-1,index) + rhs(index+1,index)+...
                                rhs(index,index-1) + rhs(index,index+1));
        case 'Order2'
            h = 1/N;
            U = zeros(N+1);
            %grid
            x = 0:h:1;
            y = 0:h:1;
            [X,Y] = meshgrid(x,y);
            U(:,1) = g(X(:,1),Y(:,1));           
            U(:,end) = g(X(:,end),Y(:,end));
            U(1,:) = g(X(1,:)',Y(1,:)');
            U(end,:) = g(X(end,:)',Y(end,:)');
            
            rnd = rand(N-1);
            U(2:N,2:N) = rnd(1:end,1:end);
            %Right Hand side
            rhs = zeros(N+1);
            rhs(:,:) = f(X,Y);
        case 'Solution'
%             x = linspace(0,pi,N+1);      
%             h = abs(x(2)-x(1))
            h = 1/N;
            x = 0:h:1;
            U = zeros(N+1);
            %grid
            y = x;
            [X,Y] = meshgrid(x,y);
            U(:,1) = g(X(:,1),Y(:,1));            
            U(:,end) = g(X(:,end),Y(:,end));
            U(1,:) = g(X(1,:)',Y(1,:)');
            U(end,:) = g(X(end,:)',Y(end,:)');          
            rnd = rand(N-1);
            U(2:N,2:N) = rnd(1:end,1:end);
            %Right Hand side
            rhs = zeros(N+1);
            rhs = f(X,Y);
            index = 2:N;
            rhs(index,index) = (h^2/12 - 3 - h^2*k^2/12)*rhs(index,index) +...
                                h^2/12*(rhs(index-1,index) + rhs(index+1,index)+...
                                rhs(index,index-1) + rhs(index,index+1));
    end
end