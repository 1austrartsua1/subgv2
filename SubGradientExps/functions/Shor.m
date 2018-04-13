function [d_sq,f,xout] = Shor(c_theta,B,p_1,E,b,tau,x_init,T,x_opt)
                           
   phi = acos(c_theta/B);


%step2 compute h_1
% and step3 compute r(phi)
if(phi<pi/4)
    h_1 = B*sqrt(p_1)/(2*c_theta);
    r_phi = 0.5*c_theta/B;
else
    h_1 = sqrt(p_1)*c_theta/B;
    r_phi = sin(phi);
end

h_k=h_1;
x_k=x_init;
for k=1:(T)
   h_kp1 = h_k*r_phi; 
   
   g_k = E'*sign(E*x_k-b);
   x_kp1 = L1BallProjection(x_k - h_kp1*g_k/norm(g_k,2),tau);       
   x_k = x_kp1;  
   f(k) = norm(E*x_k-b,1); 
  d_sq(k)= norm(x_k-x_opt,2)^2;
  hs(k)=h_kp1;
  h_k = h_kp1;
end
xout = x_k;