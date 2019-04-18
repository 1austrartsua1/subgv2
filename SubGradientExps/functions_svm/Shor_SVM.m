function [d_sq,f,xout] = Shor_SVM(c,G,Omega_1,tau,x_1,T,x_opt,C,y,m)
phi = acos(c/G);
%step2 compute h_1
% and step3 compute r(phi)
if(phi<pi/4)
    h_1 = G*sqrt(Omega_1)/(2*c);
    r_phi = 0.5*c/G;
else
    h_1 = sqrt(Omega_1)*c/G;
    r_phi = sin(phi);
end

h_k=h_1;
x_k=x_1;
for k=1:(T)
   h_kp1 = h_k*r_phi; 
   
   g_k = svm_subg(C,y,x_k);
   x_kp1 = L1BallProjection(x_k - h_kp1*g_k/norm(g_k,2),tau);       
   x_k = x_kp1;  
   f(k) = sum((ones(m,1)-(y').*(C'*x_k) >0).*(ones(m,1)-(y').*(C'*x_k)));
   d_sq(k)= norm(x_k-x_opt,2)^2;
  hs(k)=h_kp1;
  h_k = h_kp1;
end
xout = x_k;