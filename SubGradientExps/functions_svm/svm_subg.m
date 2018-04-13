function [g] = svm_subg(C,y,x)
m=length(y);
Y=diag(y);
g = -C*Y*(ones(m,1) - Y*C'*x>0);
end

