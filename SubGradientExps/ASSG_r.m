function [w,err] = ASSG_r(K,t,eps_0,eps,theta,c,A,b,w_0)
[m,n]=size(A);
lambda=eps^(2*(1-theta))/(2*eps_0*c^2);
w=[];
err=[];
w_opt = A\b;
w=w_0;
for k=1:K
   [w,err_new] = SSGS(w,lambda,t,m,n,w_opt);    
   err=[err err_new]; 
end
   
end


function [w_out,err] = SSGS(w_in,lambda,t,A,b,m,n)
    w=w_in
    for i=1:t
        w_prime = (1-2/i)*w + (2/i)*w_in - (2/(lambda*i))*subg(w,A,b,m,n);
        w_out = ((i-1)/i)*w_out + w_prime/i;
        err(i)=norm(w_out-w_opt)^2;
    end
end

function g=subg(w,A,b,m,n)
   j=randsample(m,1);
   g = A(j,:)'*(A(j,:)*w - b(j));
end

%what the funk