function [dsq,f_RSG,x] = RSG4SVM(K_rsg,alpha,eta,t,x_1,x_opt,C,y,tau,m)
x=x_1;
dsq=[];
f_RSG=[];
for k=1:K_rsg
   [x,dsq_s,f_s] = SG(eta,t,x,x_opt,C,y,tau,m); 
   eta = eta/alpha;
   dsq=[dsq dsq_s];
   f_RSG=[f_RSG f_s];
end
end

function [xhat,dsq_s,f_s] = SG(eta,t,x_1,x_opt,C,y,tau,m)
x=x_1;
xhat=x;
   for i=1:t
      g = svm_subg(C,y,x);
      x = L1BallProjection(x-eta*g,tau);
      dsq_s(i) = norm(x-x_opt,2)^2;
      f_s(i) = sum((ones(m,1)-(y').*(C'*x) >0).*(ones(m,1)-(y').*(C'*x))); 
      xhat = xhat+x;
   end
   xhat=xhat/t;
   

end