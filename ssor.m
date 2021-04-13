% iterative method in component form

function x = ssor(A,b,xprev,omega,tol,max_iter)
  n = length(A);
  iter = 0;
  epsi = norm(A*xprev-b);
  while (epsi > tol) || (iter < max_iter)
    for i = 1:n
      sigma = 0;
      for j = 1:i-1
        sigma = sigma + A(i,j)*xprev(j);
      endfor
      for k = i+1:n
        sigma = sigma + A(i,k)*xprev(k);
      endfor
      sigma = (b(i)-sigma)/A(i,i);
      x(i) = omega*sigma + (1-omega)*xprev(i);
    endfor
    
    for i=n:-1:1
      sigma = 0;
      for j = 1:i-1
        sigma = sigma + A(i,j)*x(j);
      endfor
      for k = i+1:n
        sigma = sigma + A(i,k)*x(k);
        
      endfor
      sigma = (b(i)-sigma)/A(i,i);
      x(i) = omega*sigma + (1-omega)*x(i);
    endfor
    
    xprev = x';
    epsi = norm(A*xprev-b);
    iter = iter+1;
    %fprintf('%11.0f %4.10f %4.10f %4.10f %4.10f\n',[iter;x'] )
    iter
    x'
  endwhile
endfunction
