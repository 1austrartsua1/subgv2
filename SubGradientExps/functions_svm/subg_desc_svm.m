function [f,dsq] = subg_desc_svm(m,x_opt,C,y,tau,p,alpha_1,x_1,T)

x=x_1;
for i=1:T
   alpha_i = alpha_1*i^(-p); 
   g_i = svm_subg(C,y,x);
   x = L1BallProjection(x-alpha_i*g_i,tau);
   %x = x-alpha_i*g_i;
   f(i) = sum((ones(m,1)-(y').*(C'*x) >0).*(ones(m,1)-(y').*(C'*x)));
   dsq(i) = norm(x-x_opt,2)^2;
    
end

