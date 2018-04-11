function [x,f,d_opt_sq,f_erg,sparsity,sparsity_erg,d_erg_sq,xhat] = ...
DS_SG_RobustRegre_L1Constraint(beta,M,alpha_1,K,x_init,E,b,tau,x_opt,Ergodic)
                              
x=x_init;
xhat=x;
alpha=alpha_1;
i_sgs=1;
f_erg=[];
sparsity_erg=[];
d_erg_sq=[];
for m=1:M
    if(Ergodic)
       x=xhat; 
    end
    
    for k=1:K
       subg = E'*sign(E*x-b);
       x = L1BallProjection(x-alpha*subg,tau);   
       if(Ergodic)
          xhat = (k/(k+1))*xhat + (1/(k+1))*x; 
          f_erg(i_sgs)=norm(E*xhat-b,1);
          sparsity_erg(i_sgs) = sum(abs(xhat)>1e-4);
          d_erg_sq(i_sgs) = norm(xhat-x_opt,2)^2;
       end
       
       d_opt_sq(i_sgs) = norm(x-x_opt,2)^2;
       sparsity(i_sgs) = sum(abs(x)>1e-4);
       f(i_sgs) = norm(E*x-b,1);
       i_sgs=i_sgs+1;
       %if(f(i_sgs)-f_LB<=sqrt((p_1*(c_theta^2))/(beta^m)))
       %    flag(i_sgs)=1;
       %    break;
       %else
       %    flag(i_sgs)=0;
       %end
       %start=0;
      % if((d_opt_sq(i_sgs-1)<=p_1/(beta^m)))
       %   flag(i_sgs)=1;
       %   if(Early_stop && i_sgs>start)
       %      break; 
       %   end
       %else
       %   flag(i_sgs)=0; 
       %end
    end
    alpha=alpha/sqrt(beta);
    
end
