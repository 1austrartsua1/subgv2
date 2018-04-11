function [dsq_r2sg,f_r2sg] = R2SG(beta,M,alpha_1,K_init,x_init,E,b,tau,x_opt,num_rounds,increase_factor)

K=K_init;
dsq_r2sg=[];f_r2sg=[];
for i=1:num_rounds
   [x,f,d_sq,f_ergs,sparsity0,sparsity_erg,d_sq_rsg,xhat] = ...
   DS_SG_RobustRegre_L1Constraint(beta,M,alpha_1,K,x_init,E,b,tau,x_opt,1);
   x_init=xhat;
   K=ceil(K*increase_factor);
   dsq_r2sg = [dsq_r2sg d_sq_rsg];
   f_r2sg = [f_r2sg f_ergs];
end

end