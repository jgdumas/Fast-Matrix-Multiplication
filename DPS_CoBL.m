function M = DPS_CoBL(A, nmin)
if nargin < 2, nmin = 8; end
n = length(A);
if n ~= 2^( log2(n) )
   error('The matrix dimension must be a power of 2.')
end
if n <= nmin
   M=A;
else
   m = n/2; i = 1:m; j = m+1:n;

   SQRT3o2=sqrt(3)/2;
   SQRT3o3=sqrt(3)/3;
   SQRT3o6=sqrt(3)/6;

   M1 = DPS_CoBL(A(i,i), nmin);
   M2 = DPS_CoBL(A(j,i), nmin);
   M3 = DPS_CoBL(A(i,j), nmin);
   M4 = DPS_CoBL(A(j,j), nmin);

   T1 = (M3-M2)/2;
   T2 = M4*SQRT3o3;
   T3 = T1+M1*SQRT3o2;

   S1 = M4*SQRT3o2+T3;
   S2 = M2-T2;
   S3 = M3+T2;
   S4 = M4*SQRT3o6-T3;

   M = [ S1 S2; S3 S4];
end
