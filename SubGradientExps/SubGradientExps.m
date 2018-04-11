%% Synthetic Experiments for Subgradient Method 

%% Regular Regression
clear all
m=100;
d=200;
E=randn(m,d);
b=randn(m,1);
% solve min_x (0.5/m)*norm(E*x-b,2)^2

v=eig(E'*E);
inds=find(v>1e-3);
lambdamin=v(inds(1));
c = lambdamin/(2*m);
x_start=randn(d,1);
xopt=E\b;
fstar = (0.5/m)*norm(E*xopt-b,2)^2;
D=norm(xopt,2)*10;
n_epochs=100;
T=n_epochs*m;
%=================
%% SG
%=================


A=10/c;
d=2;

x=x_start;
for k=1:T
   i=randsample(m,1);
   gradx=(E(i,:)*x - b(i))*E(i,:)';
   alpha = A*(k+1)^(-d);
   before_proj=x-alpha*gradx;
   if(norm(before_proj,2)>=D)
       x=before_proj/norm(before_proj,2);
   else
       x=before_proj;
   end
   f(k)=(1/(2*m))*norm(E*x-b,2)^2;
   dists(k)=norm(x-xopt,2);
end


%% function values
loglog(1:T,f-fstar);

%% SG plot distance to optimal

plot(log(1:T),log(dists.^2))
loglog(1:T,dists.^2)


%% Robust Regression

%% Regbust Regression