function [x,f,d_sq,fhat,d_sq_av] = decayStochasticSG_RegRegr_L2Constraint(T,x_init,A,d,E,b,tau,x_opt,best_step,c_theta)
x=x_init;xhat=x_init;
n=length(x);
m=length(b);
for i=1:T
    if(best_step)
       alpha = (2*i+1)/(2*c_theta*(i+1)^2);
    else
       alpha = A*i^(-d); 
    end
   j=randsample(m,1);
   stoc_subg = (E(j,:)*x - b(j))*E(j,:)';
   x = ProjL2_ball(x-alpha*stoc_subg,tau); 
   xhat = (i/(i+1))*xhat + (1/(i+1))*x;
   f(i) = (1/(2*m))*norm(E*x-b,2)^2; 
   fhat(i) = (1/(2*m))*norm(E*xhat-b,2)^2; 
   d_sq(i)= norm(x-x_opt,2)^2;
   d_sq_av(i) = norm(xhat-x_opt,2)^2;
end

end