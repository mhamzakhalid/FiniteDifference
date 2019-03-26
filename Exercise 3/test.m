clear 


%% Testcases 
 %testcase = 'Convergence'
 testcase = 'Solution'

%% Code starts here

switch testcase
    case 'Convergence'
        load('refsol.mat')
        %% Setting up the Domain
        for stp = 3:7
             N = 2^stp-1;
            d.x=linspace(-3,3,N);
            d.y=linspace(-3,3,N);
            [X,Y]=meshgrid(d.x,d.y);
            d.omega=((abs(X)<3)&(Y<3)&(Y>1))|((abs(X)<3)&(Y>-3)&(Y<-1))|...
                    ((abs(Y)<=1)&(X<3)&(X>1))|((abs(Y)<=1)&(X>-3)&(X<-1));
            d.bound=((abs(Y)<=3)&((X==3)|(X==-3)))|((abs(X)<=3)&((Y==3)|(Y==-3))); 
            r=0;
            d.rr=zeros(N,N);
                for i=1:N
                    for j=1:N
                        if d.omega(i,j)
                            r=r+1;
                            d.ii(r)=i;
                            d.jj(r)=j;
                            d.rr(i,j)=r;    
                        end 

                    end
                end
            d.NU=r;

            %% Setting Up Matrix A
            d.L = sparse(d.NU);
            d.b = zeros(d.NU,1);
            for r=1:d.NU
                d.L(r,r) = -4;
                i = d.ii(r);
                k = d.jj(r);
                if d.bound(i-1,k)
                    d.b(r) = 1;
                end
                if d.bound(i+1,k)
                    d.b(r) = 1;
                end
                if d.bound(i,k-1)
                    d.b(r) = 1;
                end
                if d.bound(i,k+1)
                    d.b(r) = 1;
                end
                if d.omega(1+i,k)
                    d.L(r,d.rr(i+1,k))=1;
                end
                if d.omega(-1+i,k)
                    d.L(r,d.rr(i-1,k))=1;
                end
                if d.omega(i,1+k)
                    d.L(r,d.rr(i,k+1))=1;
                end
                if d.omega(i,-1+k)
                    d.L(r,d.rr(i,k-1))=1;
                end
                 perc =r/d.NU*100;
        %          fprintf('Percentage: %f\n',perc);
             end
        h(stp-2)= abs(d.x(2)-d.x(1));

        d.b(1)=2;
        d.b(N-2)=2;
        d.b(d.NU-N+3)=2;
        d.b(d.NU)=2;

        d.b=d.b;
        d.L=-d.L;

        d.U = d.L\d.b;
        d.U1 = zeros(size(d.rr,1));
        for i=1:N
                for j=1:N
                    if d.omega(i,j)
                        d.U1(i,j)=d.U(d.rr(i,j));
                    end
                    if d.bound(i,j)
                        d.U1(i,j) = 1;
                    end
                    d.U1;
                end
        end
        % surf(d.x,d.y,d.U1)

        %% Error Check
        dx = d.x;
        U1 = d.U1;
         c=find(ismember(dx,dxref));
         d=find(ismember(dxref,dx));
         e= U1(2,c)-Uref(2,d);
        clear c d
        error(stp-2)=h(stp-2)^(1/2)*norm(e);
        end
        loglog(h,error,'-*')
        grid on
        xlabel('stepsize')
        ylabel('Error')
        
    case 'Solution'
         N = 2^5-1;
            d.x=linspace(-3,3,N);
            d.y=linspace(-3,3,N);
            [X,Y]=meshgrid(d.x,d.y);
            d.omega=((abs(X)<3)&(Y<3)&(Y>1))|((abs(X)<3)&(Y>-3)&(Y<-1))|...
                    ((abs(Y)<=1)&(X<3)&(X>1))|((abs(Y)<=1)&(X>-3)&(X<-1));
            d.bound=((abs(Y)<=3)&((X==3)|(X==-3)))|((abs(X)<=3)&((Y==3)|(Y==-3))); 
            r=0;
            d.rr=zeros(N,N);
                for i=1:N
                    for j=1:N
                        if d.omega(i,j)
                            r=r+1;
                            d.ii(r)=i;
                            d.jj(r)=j;
                            d.rr(i,j)=r;    
                        end 

                    end
                end
            d.NU=r;

            %% Setting Up Matrix A
            d.L = sparse(d.NU);
            d.b = zeros(d.NU,1);
            for r=1:d.NU
                d.L(r,r) = -4;
                i = d.ii(r);
                k = d.jj(r);
                if d.bound(i-1,k)
                    d.b(r) = 1;
                end
                if d.bound(i+1,k)
                    d.b(r) = 1;
                end
                if d.bound(i,k-1)
                    d.b(r) = 1;
                end
                if d.bound(i,k+1)
                    d.b(r) = 1;
                end
                if d.omega(1+i,k)
                    d.L(r,d.rr(i+1,k))=1;
                end
                if d.omega(-1+i,k)
                    d.L(r,d.rr(i-1,k))=1;
                end
                if d.omega(i,1+k)
                    d.L(r,d.rr(i,k+1))=1;
                end
                if d.omega(i,-1+k)
                    d.L(r,d.rr(i,k-1))=1;
                end
                 perc =r/d.NU*100;
        %          fprintf('Percentage: %f\n',perc);
             end
        h= abs(d.x(2)-d.x(1));

        d.b(1)=2;
        d.b(N-2)=2;
        d.b(d.NU-N+3)=2;
        d.b(d.NU)=2;

        d.b=d.b;
        d.L=-d.L;

        d.U = d.L\d.b;
        d.U1 = zeros(size(d.rr,1));
        for i=1:N
                for j=1:N
                    if d.omega(i,j)
                        d.U1(i,j)=d.U(d.rr(i,j));
                    end
                    if d.bound(i,j)
                        d.U1(i,j) = 1;
                    end
                    d.U1;
                end
        end
        surf(d.x,d.y,d.U1)
        colormap summer;
end


