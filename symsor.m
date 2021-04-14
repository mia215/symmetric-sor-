% iterative method in matrix form

function x = symsor(A,b,xprev,omega,tol,max_iter) 
  
  D = diag(diag(A));
  L = tril(A)-D;
  U = triu(A)-D;
  
  epsi = norm(A*xprev-b);
  k = 0; 
  B1 = inv(D+(omega*U))*(-omega*L+((1-omega)*D));
  B2 = inv(D+(omega*L))*(-omega*U+((1-omega)*D));
  while (epsi > tol) || (k < max_iter)
     x = B1*B2*xprev + (omega*(2-omega)*inv(D+(omega*U))*D*inv(D+(omega*L))*b);
     xprev = x;
     epsi = norm(A*xprev-b);
     k = k+1;   
     %fprintf('%11.0f %4.10f %4.10f %4.10f %4.10f\n',[k;x] )
     k
     x
  endwhile
endfunction
