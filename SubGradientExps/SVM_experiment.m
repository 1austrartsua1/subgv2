%% SVM experiment
clear all
close all


%% Generate a data set
n = 50;m=100;
C=randn(n,m);
y=2*(rand(1,m)>0.5)-1;
Y=diag(y);

c=C(:,1);
y_0=1;
%%
tau = 2;
x_1=randn(n,1);
x_1 = L1BallProjection(x_1,tau);

%% CVX
cvx_begin
   variables x(n) z(m) d(m);
   maximize -sum(d);
   subject to
      z==ones(m,1)-Y*C'*x;
      d>=z;
      d>=0;
      norm(x,1)<=tau;
cvx_end
    
x_opt = x;
f_opt  = -cvx_optval;

%%
%f_optCheck = sum((ones(m,1)-(y').*(C'*x_opt) >0).*(ones(m,1)-(y').*(C'*x_opt)))
%% Subgradient Method with decaying stepsizes
T=40000;
fprintf('p=0.5 for %d iterations\n',T);
p=0.5;
alpha_1=0.01;
[f_dec05,dsq_dec05] = subg_desc_svm(m,x_opt,C,y,tau,p,alpha_1,x_1,T);
f_dec05(end)-f_opt
p=1;
alpha_1=0.1;
fprintf('p=1 for %d iterations\n',T);
[f_dec1,dsq_dec1] = subg_desc_svm(m,x_opt,C,y,tau,p,alpha_1,x_1,T);
f_dec1(end)-f_opt

%%
loglog(1:T,dsq_dec1,1:T,dsq_dec05)
legend('\alpha_k=0.1 k^{-1}','\alpha_k=0.01 k^{-0.5}')

%% DS-SG
%compute G
runDSSG=0;
if(runDSSG)
tic
G=0;
for i=1:m
    G=G+norm(C(:,i),2);
end
c = 21.5; 
beta=2;
Omega_1 = 4*tau^2;
epsilon=1e-1;
%M = ceil(log(Omega_1/epsilon)/log(beta));
%M=20;
mdim=m;

[dsq_DSSG,f_DSSG] = DS_SG_SVM(beta,M,x_1,Omega_1,G,c,x_opt,C,y,tau,mdim);
toc
semilogy(dsq_DSSG)
semilogy(f_DSSG-f_opt)
end

%% RSG
runRSG=0;
if(runRSG)
g_1 = svm_subg(C,y,x_1);
eps_0 = G*sqrt(Omega_1);
epsilon = 1e-1;
alpha = 2;
c_rsg=c;
t = ceil((1/c_rsg^2) * alpha^2 * G^2);
K_rsg = ceil(log(eps_0/epsilon)/log(alpha));
eta = eps_0/(alpha*G^2);
fprintf('Number of RSG iterations %d\n',K_rsg*t);

[dsq_RSG,f_RSG] = RSG4SVM(K_rsg,alpha,eta,t,x_1,x_opt,C,y,tau,m);
plot(dsq_RSG)
end
%% Shor

G=0;
for i=1:m
    G=G+norm(C(:,i),2);
end
Omega_1 = 4*tau^2;
runShor=0;
if(runShor)
T=40000;
c_shor =40;
[dsq_shor,f_shor] = Shor_SVM(c_shor,G,Omega_1,tau,x_1,T,x_opt,C,y,m);
end
%% Shor2
c_0 = G/2;
M = 5;
eps = 1e-8;
           
[dshor2,fshor2] = metaShorSVM(c_0,G,Omega_1,C,y,m,tau,x_1,x_opt,M,eps);
semilogy(fshor2 - f_opt)

%% DS2-SG
G=0;
for i=1:m
    G=G+norm(C(:,i),2);
end
c_1 = G;
beta=4;
Omega_1 = 4*tau^2;
epsilon=1e-8;
mdim=m;
M = ceil(log(Omega_1/epsilon)/log(beta));
M=20;
L=5;
f_DS2SG=[];
dsq_DS2SG=[];
stepLength=[];
c=c_1;
x=x_1;
for l=1:L
    
   [dsqs,fs,x] = DS_SG_SVM(beta,M,x,Omega_1,G,c,x_opt,C,y,tau,mdim); 
   dsq_DS2SG=[dsq_DS2SG dsqs];
   f_DS2SG = [f_DS2SG fs];
   c=c/2;
end

plot(f_DS2SG)
%% R2SG
thetahat=0;
S=9;
c=G;
g_1 = svm_subg(C,y,x_1);
eps_0 = G*sqrt(Omega_1);
epsilon = 1e-1;
alpha = 2;
c_rsg=c;
t = ceil((1/c_rsg^2) * alpha^2 * G^2);
K_rsg = ceil(log(eps_0/epsilon)/log(alpha));
K_rsg=20;
eta = eps_0/(alpha*G^2);
x_R2SG=x_1;
f_R2SG=[];
dsq_R2SG=[];
total=0;
for s=1:S
    total=total+ K_rsg*t;
    fprintf('Number of R2SG so far iterations %d\n',total);
   [dsqs,fs,x_R2SG] = RSG4SVM(K_rsg,alpha,eta,t,x_R2SG,x_opt,C,y,tau,m); 
   f_R2SG = [f_R2SG fs];
   dsq_R2SG=[dsq_R2SG dsqs];
   t = ceil(t*2^(1-thetahat));
end



%% plot distances
%Tmin = min([length(dsq_DSSG) length(dsq_dec1) length(dsq_dec05) length(dsq_RSG)]);
%semilogy((1:Tmin),dsq_DSSG(1:Tmin),1:Tmin,dsq_dec1(1:Tmin),...
%    1:Tmin, dsq_dec05(1:Tmin),...
%    1:Tmin, dsq_RSG(1:Tmin),1:Tmin,dsq_shor(1:Tmin),1:Tmin,dsq_DS2SG(1:Tmin))
%legend('DSSG','p=1','p=0.5','RSG','Shor','DS2-SG')

%% plot function values
Tmin = min([length(f_dec05) length(f_dec1) ...
     length(f_DS2SG) length(f_R2SG) length(fshor2)]);
f_dec05_Up = mycummin(f_dec05(1:Tmin));
f_dec1_Up = mycummin(f_dec1(1:Tmin));
f_DS2SG_Up = mycummin(f_DS2SG(1:Tmin));
f_R2SG_Up = mycummin(f_R2SG(1:Tmin));
f_shor2_Up = mycummin(fshor2(1:Tmin));
f_opt = min([f_opt f_DS2SG fshor2]);
target = f_opt;
h = semilogy(1:Tmin,f_dec05_Up-target,...
         1:Tmin,f_dec1_Up-target,...
         1:Tmin,f_DS2SG_Up-target,...
         1:Tmin,f_shor2_Up-target,...
         1:Tmin,f_R2SG_Up-target,'linewidth',3);
grid on 
axis([0 Tmin 1e-9 1e1])
legend('\alpha_k=0.01k^{-0.5}','\alpha_k=0.1 k^{-1}'...
    ,'DS2-SG','Shor2','R2SG');
xlabel('Number of iterations')
ylabel('$h(x_k) - h^*$',...
'Interpreter', 'Latex', 'FontSize', 15, 'Color', 'k');

%axis([0 2e4 1e-9 1e1])
