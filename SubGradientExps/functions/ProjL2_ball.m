function [out]=ProjL2_ball(x,tau)

if(norm(x,2)<=tau)
    out=x;
else
    out = tau*x/norm(x,2);
end
end
