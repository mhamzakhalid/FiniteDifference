clear
N = 5;
n=2^N;
h = 1/n;
k = 1;
% Returns uniform triangulation of unit square
[X,Y] = meshgrid(0:1/n:1,0:1/n:1);
X = reshape(X,[],1);
Y = reshape(Y,[],1);
% we have a vector of X- and Y- coordinates of our vertices
P = [X,Y];
%Doing the triangulation
tri = delaunay(X,Y);
TR = triangulation(tri,P);

if N==2 
%Displyaing Vertices
triplot(TR)
hold on
vxlabels = arrayfun(@(n) {sprintf('V%d', n)}, (1:length(X))');
Hpl = text(X, Y, vxlabels, 'FontWeight', 'bold', 'HorizontalAlignment',...
   'center', 'BackgroundColor', 'none');
ic = incenter(TR);
numtri = size(TR,1);
trilabels = arrayfun(@(x) {sprintf('K%d', x)}, (1:numtri)');
Htl = text(ic(:,1), ic(:,2), trilabels, 'FontWeight', 'bold', ...
   'HorizontalAlignment', 'center', 'Color', 'blue');
end


%Initial Values
S=zeros(length(P));
M=zeros(length(P));
F=zeros(length(P),1);


%Defining F
f = @(x,y) 20*pi^2*sin(2*pi*x).*cos(4*pi*y); 

GradPhi = [-ones(2,1),eye(2)]';

for i=1:length(tri)
    v1=P(tri(i,1),:);
    v2=P(tri(i,2),:);
    v3=P(tri(i,3),:);
    J = [v2(1) - v1(1) , v3(1) - v1(1);
        v2(2) - v1(2) , v3(2) - v1(2)];
    G=GradPhi*mldivide(J,eye(2));
    %Setting up the stiffness matrix
    Ak=det(J)/2*(G*G');
    %Setting up the mass Matrix
    Mk = k^2*det(J)*computeMiniMassMatrix(v1',v2',v3');
    %Setting up the right hand side
    Fk = getrhs(v1',v2',v3',f)*det(J)/2;
    
   %Transporting local matrix to the global matrix
    S(tri(i,:),tri(i,:)) = S(tri(i,:),tri(i,:)) + Ak;
    M(tri(i,:),tri(i,:)) = M(tri(i,:),tri(i,:)) + Mk;
    F(tri(i,:),1) = F(tri(i,:),1) + Fk;
    
          clear Ak Fk Mk;
end
% System Setup
A = S - M;

%% Implementing the Boundary Conditions

%Finding boundary nodes for Dirichlet Left and Right
boundary = freeBoundary(TR);
p = [boundary(:,1);boundary(:,2)];
for i=1:length(p)
    A(p(i,1),:)=0;
    A(p(i,1),p(i,1))=1;
    F(p(i,1))=0;
end
%%  Solution
U = A\F;

%% plots

% figure
pn = [P,U];
TRn = triangulation(tri,pn);
trisurf(TRn)
view([20 30])
axis off  