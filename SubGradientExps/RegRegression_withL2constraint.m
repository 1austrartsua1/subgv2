%% Regular regression with an $\ell_2$ constraint
%solve min_x (1/(2*m))*norm(E*x-b,2)^2: norm(x,2)<=tau
clear all
close all
at_CSL=1;
if(at_CSL)
   addpath(genpath('C:\Users\prjohns2.UOFI\Dropbox\MATLAB_F14on\Work_By_Semester\Fall16\SubGradientExps'));
else
   addpath(genpath('C:\Users\Pat\Dropbox\MATLAB_F14on\Work_By_Semester\Fall16\SubGradientExps')); 
end
m=100;
n=50;
tau=1;
E=randn(m,n);
b=randn(m,1);
v=svd(E);
c_theta = (1/(2*m))*v(end);
normE = norm(E,'fro');
B_sq=(1/m)*(normE)^(3/2)*((tau^2)*normE^(5/2) + 2*tau*norm(b,1)+(normE^2)*norm(b,2)^2);

x_init = randn(n,1);
x_init = ProjL2_ball(x_init,tau);
%% CVX to get x_opt

cvx_begin
   variable x_opt(n);
   maximize -(1/(2*m))*norm(E*x_opt - b,2);
   subject to
      norm(x_opt,2)<=tau;

cvx_end
f_opt = cvx_optval^2;

%% decaying subgradient descent

tic
T=1e5;
A=1/c_theta;
d=1;
best_step=1;

[x,f_decay1,d_sq1] = decayStochasticSG_RegRegr_L2Constraint(T,x_init,A,d,E,b,tau,x_opt,best_step,c_theta);

best_step=0;
[x,f_decay2,d_sq2] = decayStochasticSG_RegRegr_L2Constraint(T,x_init,A,d,E,b,tau,x_opt,best_step,c_theta);

A=1;
d=0.5;
[x,f_decay3,d_sq3,fhat,d_sq3Hat] = decayStochasticSG_RegRegr_L2Constraint(T,x_init,A,d,E,b,tau,x_opt,best_step,c_theta);

A=1;
d=0.75;
[x,f_decay4,d_sq4] = decayStochasticSG_RegRegr_L2Constraint(T,x_init,A,d,E,b,tau,x_opt,best_step,c_theta);

%plot(1:T,f_decay1-f_opt)
%plot(d_sq)
toc

%%
loglog(1:T,d_sq1,1:T,d_sq2,1:T,d_sq3Hat)
legend('(A,d)=(1/c_\theta,1)','Prop. 1','(A,d)=(1,0.5)');
grid on
axis([1 T 10^(-3) 10^1]);
xlabel('Subgradient evals');
ylabel('$d(x_k,\mathcal{X})^2$',...
'Interpreter', 'Latex', 'FontSize', 15, 'Color', 'k');
%% constant stepsize up to some tolerance?

