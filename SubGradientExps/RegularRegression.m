% regular regression
clear all
close all
%goodbyeCruelWorld
m=550;
d=100;
B=randn(m,d);
b=randn(m,1);

v=svd(B);
c_theta=v(end)^2;
%min_w (0.5)*norm(Aw-b,2)^2


max_so_far=0;
for i=1:m
    normNow=norm(B(i,:),2);
   if(normNow>max_so_far)
       max_so_far=normNow;
   end
end
A_2 = max_so_far;

D=10;
w_opt = B\b;

w_0=randn(d,1);
eps_0 = (1/(2*m))*norm(B*w_0 - b,2)^2;

%% Stochastic Approximation
tic
T=1e4;
w=w_0;
A=1e-1;
d=0.7;
for k=1:T
   j=randsample(m,1);
   g = B(j,:)'*(B(j,:)*w - b(j));
   alpha = (2*k+1)/(2*c_theta*(k+1)^2);
   alpha = A*k^(-d);
   w = w-alpha*g; 
   f(k)=(1/(2*m))*norm(B*w-b,2)^2;
   dsquare(k) = norm(w-w_opt,2)^2;
end
toc
%%
plot(f)
loglog(1:T,f)

%%
plot(dsquare)
loglog(1:T,dsquare)

%% ASSG_r
theta=0.5;
eps=1e0;
eps_0 = (1/(2*m))*norm(B*w_0 - b,2)^2;
K=log2(eps_0/eps);
c=(1/c_theta)^theta;
lambda_1= eps^(2*(1-theta))/(2*c^2*eps_0);
G= D*A_2+norm(b,2);
t_0=136*G^2;
delta_tilde = 1/(2*exp(1)*K);
t=136*G^2*(1+log(4*log(t_0/delta_tilde))+log(t_0));

[w,err] = ASSG_r(K,t,eps_0,eps,theta,c_theta,lambda_1,B,b,w_0);

%%
x=1.1:.01:10;

plot(x,x./log(x))

%%
beta=1:0.1:10;
F=sqrt(2*beta).*log(2*beta);
plot(beta,F);


