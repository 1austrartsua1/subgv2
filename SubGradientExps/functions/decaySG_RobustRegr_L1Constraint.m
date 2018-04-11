function [x,f,d_sq,d_sq_hat] = decaySG_RobustRegr_L1Constraint(T,x_init,A,d,E,b,tau,x_opt)
x=x_init;
xhat=x;
for i=1:T
   alpha = A*i^(-d); 
   subg = E'*sign(E*x-b);
   x = L1BallProjection(x-alpha*subg,tau); 
   xhat = (i/(i+1))*xhat + (1/(i+1))*x;
   f(i) = norm(E*x-b,1); 
   d_sq(i)= norm(x-x_opt,2)^2;
   d_sq_hat(i)=norm(xhat-x_opt,2)^2;
end

end