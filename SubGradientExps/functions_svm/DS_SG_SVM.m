function [dsq,f,xhat, stepLength] = DS_SG_SVM(beta,M,x_1,Omega_1,G,c,x_opt,C,y,tau,mdim)

K = ceil(sqrt(3/2)*(G/c)^2 * sqrt(beta) * log(3*beta));
alpha=(2/sqrt(3)) * c/(G^2) * sqrt(Omega_1/beta);

%K = ceil((G/c)^2 * sqrt(beta) * log(2*beta));
%alpha=2*c/(G^2) * sqrt(Omega_1/(2*beta));
fprintf('Total DS-SG Iterations is %d\n', M*K);

xhat = x_1;
dsq=[];
f=[];
stepLength=[];

for m=1:M
   [xhat dsq_new, f_new,stepLengths] = fixedSG(K,alpha,xhat,x_opt,C,y,tau,mdim);
   alpha = alpha/sqrt(beta);
   dsq=[dsq dsq_new];
   f = [f f_new];
   stepLength = [stepLength stepLengths];
end

end

function [x dsqs fs stepLengths] = fixedSG(K,alpha,xhat,x_opt,C,y,tau,mdim)
x=xhat;
for i=1:K
    g = svm_subg(C,y,x);
    xold = x;
    x = L1BallProjection(x-alpha*g,tau);
    stepLengths(i) = norm(x-xold,2);
    dsqs(i) = norm(x-x_opt,2)^2;
    fs(i) = sum((ones(mdim,1)-(y').*(C'*x) >0).*(ones(mdim,1)-(y').*(C'*x)));
end

end