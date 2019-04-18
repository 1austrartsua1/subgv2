function [d,f] = metaShor(c_0,B,p_1,E,b,tau,x_init,x_opt,M,eps)
   c = c_0;
   d = [];
   f = [];
   x = x_init;
   phi = acos(c/B);
   T = ceil((log((sin(phi))^(-1)))^(-1)*log(p_1/(cos(phi)*eps)));
   for i=1:M
       T
       fprintf('running iteration %d\n',i);
       
       [dnew,fnew,xnew] = Shor(c,B,p_1,E,b,tau,x,T,x_opt);
          
       d = [d, dnew];
       f = [f, fnew];
       c = c/2;
       phi = acos(c/B);
       T = ceil((log((sin(phi))^(-1)))^(-1)*log(p_1/(cos(phi)*eps)));
       x = xnew;  
       
   
   end

