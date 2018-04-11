%% Synthetic Robust Regression Experiment
clear all
m=100;
dim=50;
E=randn(m,dim);
b=randn(m,1);
% solve min_x (1/m)*norm(E*x-b,1)

v=eig(E'*E);
inds=find(v>1e-3);
lambdamin=v(inds(1));
c = lambdamin/(2*m);
x_start=randn(dim,1);
%xopt=E\b;
%fstar = (0.5/m)*norm(E*xopt-b,2)^2;
D=100;
n_epochs=1;
T=n_epochs*m;
T=1e3;

%% CVX

cvx_begin
   variable x(dim);
   maximize -(1/m)*norm(E*x-b,1);
cvx_end
xopt=x;
%% Subgradient with decaying stepsize
tic
A=2e-1;
d=1;
decay=0;

x=x_start;
alpha=A;

for k=1:T
   %i=randsample(m,1); 
   subg=0;
   for i=1:m
      subg=subg+(1/m)*sign(E(i,:)*x - b(i))*E(i,:)';
   end
   if(decay)
   alpha = A*(k+1)^(-d);
   end
   before_proj=x-alpha*subg;
   if(norm(before_proj,2)>=D)
       x=before_proj/norm(before_proj,2);
   else
       x=before_proj;
   end
   f(k)=(1/m)*norm(E*x-b,1);
   distsquare(k) = norm(x-xopt,2)^2;
    
end
toc
%%
plot(distsquare)
semilogy(max(10^(-5),distsquare-max(distsquare(600:end))))
%%
fopt = (1/m)*norm(E*xopt-b,1);
loglog(f-fopt)

%%
loglog(distsquare)

%%

w = ASSG_r(K,t,eps_0,eps,theta,c,A,b)

plot(1:dim, xopt,1:dim,x)


