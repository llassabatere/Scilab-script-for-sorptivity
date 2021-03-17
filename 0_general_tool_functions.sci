// to erase picture
    
    function y = close_windows(pi)
      for i = 1:pi
        close;
      end;
      y = pi
    endfunction
    
// to erase files
    
    function y = erase_file(data_name)
      if MSDOS then unix('del '+data_name);       
      else unix('rm -f '+data_name); end
      y = 0;     
    endfunction
    
// to filter data (delete Nan and Inf values

    function y = nn_finite(X)
        na = find(isnan(X) == %f)';
        nb = find(isinf(X) == %f)';
        y = intersect(na,nb)';
    endfunction
    

    function y = nn_finite_pos(X)
        na = find(isnan(X) == %f)';
        nb = find(isinf(X) == %f)';
        nc = find(X>0)';
        y = intersect(intersect(na,nb),nc)';
    endfunction
    
//general tool functions 

  function y = mean_end(x,n)              // mean of the n last points
      n_end = nc_end(x,n);
      y = mean(x(n_end,:));  
    endfunction
      
  function y = nc_end(x,n)                // rank of the n last points
      nx = size(x,1);
      n_end = [nx-n+1:1:nx]
      y = n_end;  
  endfunction
  
  function y = reduc(x,n1,n2,nc,V)                // sequence
      n_ = [n1:nc:n2];
      y = x(n_,:)-V*x(n1,:);
  endfunction
   
  
  function y = deriv_sqr(t,x)                    // estimation of the derivative in regards to suqre root of time
      y = zeros(size(x,1)-1,2)
      for i = 1:(size(x,1)-1)                                        
      y(i,2) = (x(i+1)-x(i))/((t(i+1))^0.5-(t(i))^0.5);
      y(i,1) = (t(i+1)*t(i))^0.5;           
      end;
    endfunction
  
  function y = deriv(t,x)                    // estimation of the derivative in regards to suqre root of time
      y = zeros(size(x,1)-1,2)
      for i = 1:(size(x,1)-1)                                        
      y(i,2) = (x(i+1)-x(i))/(t(i+1)-t(i));
      y(i,1) = (t(i+1)+t(i))/2;           
      end;
    endfunction
    
  function y = ab_end(t,x,n)
      zx = x(nc_end(x,n),:);
      zt = t(nc_end(x,n),:);
      [a,b,sig] = reglin(zt',zx');
      y = [b;a];
  endfunction
        
  function y = section(r)
      y = %pi*r^2;
  endfunction
    
// transformation of hypermatrix into matrix 
   
     function y = matr(zz)  
       for j = 1:size(zz,3)
         y(:,j) = zz(:,:,j);
       end;
     endfunction   
     
// Interpol and discretization
   
     function y = interpol(a,x,y)
         d = splin(x,y);
         y = interp(a,x,y,d);
     endfunction
     
   function y = discret_log(Imin,Imax,Rc,nc)
     I_num__1 = logspace(log(Imin)/log(10),log(Imax/Rc)/log(10),nc)';
     I_num__2 = linspace(Imax/Rc,Imax,nc)';
     y = [0;I_num__1;I_num__2];
   endfunction
     
// reduction of data with maximum (for CI)

  function y = compteur(X,Y)
      for i = 1:size(X,2)
          k = 1;
          while X(k,i) < Y(i,:),
              k = k+1;
              if k == size(X,1) then break
              end
          end
          y(i,:) = k-1;
       end
   endfunction
   

  function y = recompose_Z(X,Y,Z)
    aa = X(1:Z(1,:),:);
    aaa = Y(1:Z(1,:),1);
    k = 1;
    
    while k < size(Y,2)
        zz = max(aa)+X(2:Z(k+1,:),:);
        zzz = max(aaa)+Y(2:Z(k+1,:),k+1);
        aa = [aa;zz];
        aaa = [aaa;zzz];
        k = k+1;
    end

    y = [aa aaa];
   endfunction
   
   function y = filtre_0(Y)
       k = 2;
        while Y(k,:) <> 0
            k = k+1;
            if k == size(Y,1)+1 then break end
        end
        y = k-1;
   endfunction

  function y = filtre_0_V(Y)
      for i = 1:size(Y,2)
          y(i,:) = filtre_0(Y(:,i))
      end
   endfunction
   
// to reject values

  function y = reject_val(X,val_reject)
    
    a = cell();
    
    for i = 1:size(val_reject,1)
        a{i} = find(X == val_reject(i))';
    end
    
    b = [];
    
    for i = 1:size(a,1)
        b = union(b,a{i})';
    end
    
    nn_ = linspace(1,size(X,1),size(X,1))';
    
    for ii = 1:size(b,1)
        nn_(b(ii)) = 0;
    end
    
    y = find(nn_ > 0)';
    
   endfunction
   
  function y = reject_val_approx(X,val_reject)      // for approached values
    
    for i = 1:size(val_reject,1)
        [a b(i)] = min(abs(X - val_reject(i)));
    end
    
    nn_ = linspace(1,size(X,1),size(X,1))';
    
    for ii = 1:size(b,1)
        nn_(b(ii)) = 0;
    end
    
    y = find(nn_ > 0)';
    
   endfunction
   
//general tool functions 

  function y = Er(a,b)                      // error between two vectors
    delta = (a-b)'*(a-b);
    y = delta/(a'*a);
  endfunction
  
  function y = Er_V(a,b)                      // error between two vectors
    y = (a-b);
  endfunction
  
  function y = Er_Vr(a,b)                      // error between two vectors
    y = (a-b)./b;
  endfunction
  
  function y = R_Nash(a,b)                      // error between two vectors
    delta = (a-b)'*(a-b);
    vara = (a-mean(a))'*(a-mean(a));
    y = 1 - delta/vara;
  endfunction
  
  function y = R_det(a,b)                      // error between two vectors
    cov_2 = ((a-mean(a))'*(b-mean(b)))^2;
    vara = (a-mean(a))'*(a-mean(a));
    varb = (b-mean(b))'*(b-mean(b));
    y = cov_2/(vara*varb);
  endfunction  
   
  function y = RMSE(a,b)                      // error between two vectors
    delta = (a-b)'*(a-b);
    y = sqrt(delta/size(a,1));
  endfunction
  
  function y = NRMSE(a,b)                      // error between two vectors
    y = RMSE(a,b)/(max(b)-min(b));
  endfunction
  
  function y = CV_RMSE(a,b)                      // error between two vectors
    y = RMSE(a,b)/mean(b);
  endfunction

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
 
