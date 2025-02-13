function C = FMM(A, B, nmin, level)
  if nargin < 3, nmin = 3; end
  if nargin < 4, level = 3; end
  [m,k] = size(A); 
  [k,n] = size(B);
  if (m >= max(k,n))
    C=FMM_6_3_3(A,B,nmin, level);
  elseif (n >= max(m,k))
    C=FMM_3_3_6(A,B,nmin, level);
  else
    C=FMM_3_6_3(A,B,nmin, level);
  end
end
