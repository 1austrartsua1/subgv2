function [v_out] = L1BallProjection(v_in,z)
   if(norm(v_in,1)<=z)
       v_out=v_in;
       return;
   end
   v=abs(v_in);
   temp=size(v_in);n=temp(1);
   U=true(n,1);
   s=0;
   rho=0;
   u_count=n;
   %inds=1:n;
   while(u_count>0)      
      k=1; 
      while(U(k)~=1)
          k=k+1;
      end
      %U_inds = inds(U);
      %k=U_inds(k_ind);
      %G=logical(zeros(n,1));
      %L=logical(zeros(n,1));
      %g_count=1;
      %l_count=1;
     % for i=1:u_count
      G=(U & (v>=v(k)));
      L=(U & (v<v(k)));
      %    if(v(U_inds(i))>=v(k))
             %G(g_count)=U(i);
      %       G(U_inds(i))=1;
      %       g_count=g_count+1;
      %    else
             %L(l_count)=U(i);
      %       L(U_inds(i))=1;
      %       l_count=l_count+1;
      %    end
      %end
      delrho = sum(G);
      if(delrho>=1)
         %dels = sum(v(G(1:(g_count-1))));
         %G_ind = inds(G);
         dels = sum(v(G));
      else
         dels = 0;
      end
      
      if( (s+dels)-(rho+delrho)*v(k) < z)
          s=s+dels;
          rho=rho+delrho;
          U = L;
      else
          %U = setdiff(G(1:(g_count-1)),k);
          U=G;
          U(k)=0;
      end
      u_count=sum(U); 
          
   end
   theta = (s-z)/rho;
   v=(v-theta>0).*(v-theta);
   v_out = v.*sign(v_in);
end