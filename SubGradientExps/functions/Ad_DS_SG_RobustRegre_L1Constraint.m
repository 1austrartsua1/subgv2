function [xtilde,f,d,sparsity] = ...
Ad_DS_SG_RobustRegre_L1Constraint(beta,G,M,c_1,Omega,x_init,T,E,b,tau,x_opt,K)

alpha = (2/sqrt(3)) * c_1/(G^2) * sqrt(Omega/beta);

xtilde = x_init;
f=[];
d=[];
sparsity=[];
for l=1:T
    [xtilde,fupdate,d_update,ferg,sparsity_update] =...
        DS_SG_RobustRegre_L1Constraint(beta,M,alpha,K,xtilde,E,b,tau,x_opt,0);
    f=[f fupdate];
    d=[d d_update];
    sparsity = [sparsity sparsity_update];
    K=4*K;
    alpha=alpha/2;
end