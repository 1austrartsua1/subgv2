function [xbar,f_chamPock] = ChamPock_PD(tauStep,sigma,theta,x_0,y_0,T,E,b,tau)

x=x_0;y=y_0;xbar=x;
for i=1:T
    y = prox_sigafstar(y+sigma*E*xbar,sigma,b);
    xn = L1BallProjection(x-tauStep*E'*y,tau);
    xbar = xn+theta*(xn-x);
    x=xn;
    f_chamPock(i) = norm(E*xbar-b,1); 
end

end



function [y] = prox_sigafstar(yin,sigma,b)
    y = yin-sigma*b-sigma*prox_l1(yin/sigma - b,1/sigma);
end

