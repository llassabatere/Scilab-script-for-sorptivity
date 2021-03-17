// modified March 2020 for paper on sorptivity

// Saturation degree

  function y = Se_(ttheta,theta_r,theta_s)          // calcul de Se correspondant à une valeur de theta
    y = (ttheta-theta_r)/(theta_s-theta_r);
  endfunction

//////////////////////////////////////////////////////////////////////////////

//// Dimensionless functions for !! BC model !!         // 

    // Water Retention curve & Hydraulic conductivity
  
      function y = Se_hstar_BC(h,lambda) 
        for i = 1:size(h,1)
            if abs(h(i)) < 1 then y(i) = 1;
                else y(i) = (abs(h(i)))^(-lambda);
            end
        end
      endfunction
    
      function y = hstar_Se_BC(Se,lambda) 
        y=Se.^(-1/lambda);
       endfunction
      
      function  y = Kr_BC(Se,eta)
        y = Se.^eta;
      endfunction
  
     // Hydraulic diffusivity
   
      function y = D_etoile_BC(xS,lambda,eta)
        y = -1/lambda*xS.^(eta-(1/lambda+1));
      endfunction
  
    // Exact analytical expressions for cp
  
      function y=cp_function_BC_psi(psi_bc)
        y = 2+(1-psi_bc)/(5*psi_bc+1)+(1-psi_bc)/(7*psi_bc+1);
      endfunction
  
  
//////////////////////////////////////////////////////////////////////////////
  
//// Dimensionless functions for !! vGB model !!         // 

    // Water Retention curve & Hydraulic conductivity
  
      function y = Se_hstar_vGB(h,n)
        m = 1-2/n;
        y = (1+h.^n).^(-m);
      endfunction

      function y = hstar_Se_vGB(Se,n)
        y=(Se.^(-n/(n-2))-1).^(1/n)
       endfunction
        
      function  y = Kr_vGB(SSe,eta)
        y = SSe.^eta;
      endfunction
  
     // Hydraulic diffusivity
   
      function y = D_etoile_vGB(xS,n,eta)
          m = 1-2/n;
          y = (m-1)/(2*m)*xS.^(eta-(m+1)/(2*m)).*(1-xS.^(1/m)).^(-(m+1)/2);
      endfunction
      
    // exact analytical expressions for cp
 
     function y=cp_function_vGB_m(m)
       y = gamma(3/2-m/2)*mtlb_a(gamma(5/2*m+1/2)/gamma(2*m+1),gamma(7/2*m+1/2)/gamma(3*m+1));
     endfunction
        
  
//////////////////////////////////////////////////////////////////////////////
  
//// Dimensionless functions for !! vGM model !!         // 

    // Water Retention curve & Hydraulic conductivity
  
      function y = Se_hstar_vGM(h,n)             // model vG-M
        m = 1-1/n;
        y = (1+h.^n).^(-m);
      endfunction
 
      function y = hstar_Se_vGM(Se,n)             // model vG-B
        y=(Se.^(-n/(n-1))-1).^(1/n)
       endfunction
    
      function  y = Kr_vGM(SSe,n,l)
        m = 1-1/n;  
        y = SSe.^l.*(1-(1-SSe.^(1/m)).^m).^2;
      endfunction

     // Hydraulic diffusivity
   
      function y = D_etoile_vGM(xS,n,l)
          m = 1-1/n;
          y = (m-1)/m*xS.^(l-1/m).*((1-xS.^(1/m)).^(-m)+(1-xS.^(1/m)).^m-2);
      endfunction

    // exact analytical expressions for cp
    
     function y=cp_function_vGM_m(m,l)
       y1 = gamma(2-m)*(gamma(m*(l+1))/(gamma(m*l)*(m*(l+1)-1))+gamma(m*(l+2))/(gamma(m*(l+1))*(m*(l+2)-1)));
       y2 = gamma(m*(l+1))*gamma(1+m)/(gamma(m*(l+2))*(m*(l+1)-1))+gamma(m*(l+2))*gamma(1+m)/(gamma(m*(l+3))*(m*(l+2)-1));
       y3 = 1/((l+1)*m-1)+1/((l+2)*m-1);
       y = y1+(1-m)*y2-2*(1-m)*y3;
     endfunction
         
  
//////////////////////////////////////////////////////////////////////////////
  
//// Dimensionless functions for !! KG model !!         // 

    // Water Retention curve & Hydraulic conductivity
      
      function y = Se_hstar_KG(h,sigma)             // model vG-B
        for i = 1:size(h,1)
            y(i) = 1/2*erfc(log(abs(h(i)))/(sqrt(2)*sigma));
        end
      endfunction
    
      function y = hstar_Se_KG(Se,sigma)             // model vG-B
        for i = 1:size(Se,1)
            y(i) = exp(sqrt(2)*sigma*erfinv(1-2*Se(i)));
        end
       endfunction
  
      function  y = Kr_KG(Se,sigma,l)
        for i = 1:size(Se,1)
            y(i) = Se(i)^l*(1/2*erfc((erfinv(1-(2*Se(i)))+sigma/sqrt(2))))^2;
        end
      endfunction

     // Hydraulic diffusivity
  
      function y = D_etoile_KG(xS,sigma,l)
        for i = 1:size(xS,1)
            y(i) = -1/2*sqrt(%pi/2)*sigma*xS(i)^l*(erfc(erfinv(1-2*xS(i))+sigma/sqrt(2)))^2*exp((erfinv(1-2*xS(i)))^2+sqrt(2)*sigma*erfinv(1-2*xS(i)));
        end
      endfunction
  
    // Numerical computation for cp
    
      function y = S_etoile_2_int_KG(xS,Se0,Sef,sigma,l)
        y=(Sef+xS-2*Se0)*D_etoile_KG(xS,sigma,l);
      endfunction

      function y = S_etoile_2_Dt_KG_comb1(Se0,Sef,Se0b,Sefb,sigma,l)                   
        // integration between Se0 and Sef of S_etoile_2_int..
        y = intg(Se0,Sef,list(S_etoile_2_int_KG,Se0b,Sefb,sigma,l))
      endfunction
      
      function y = S_etoile_2_int_Kh2_KG(xh,Se0,Sef,sigma,l)
        Se = Se_hstar_KG(xh,sigma);
        y=(Sef+Se-2*Se0).*Kr_KG(Se,sigma,l);
      endfunction
     
      function y = S_etoile_2_Kh_KG_comb2(h0,hf,Se0,Sef,sigma,l)                   
        // integration between h0 and hf of S_etoile_2_int_Kh2..
        y = intg(h0,hf,list(S_etoile_2_int_Kh2_KG,Se0,Sef,sigma,l));
      endfunction
       
      function y = S_etoile_2_KG_comb(Se0,Sef,sigma,l)
          hec = 1;
          Sec = Se_hstar_KG(hec,sigma);
          hf = hstar_Se_KG(Sef,sigma);
          y1 = S_etoile_2_Dt_KG_comb1(Se0,Sec,Se0,Sef,sigma,l);
          y2 = S_etoile_2_Kh_KG_comb2(hec,hf,Se0,Sef,sigma,l);
          y = y1+y2;
      endfunction
  
  
  

