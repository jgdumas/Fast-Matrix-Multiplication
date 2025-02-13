function C = FMMa(A, B, nmin)
  if nargin < 3, nmin = min([6,3,3]); end
  [m,k] = size(A); 
  [k,n] = size(B);
  if (m >= max(k,n))
    C=FMMa_6_3_3(A,B,nmin);
  elseif (n >= max(m,k))
    C=FMMa_3_3_6(A,B,nmin);
  else
    C=FMMa_3_6_3(A,B,nmin);
  end
end
