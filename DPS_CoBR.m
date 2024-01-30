function M = DPS_CoBR(A, nmin)
%          Right change of basis of the sparsification
%          of DPS's algorithm via alternate basis.

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
   SQRT3o6=sqrt(3)/6;
   SQRT3f2o3=2*sqrt(3)/3;

   M1 = DPS_CoBR(A(i,i), nmin);
   M2 = DPS_CoBR(A(j,i), nmin);
   M3 = DPS_CoBR(A(i,j), nmin);
   M4 = DPS_CoBR(A(j,j), nmin);

   T1 = (M1-M4)/2;
   T2 = M3*SQRT3o6-M2*SQRT3o2;

   S1 = T1+T2;
   S2 = (M1+M4)/2+(M2-M3)*SQRT3o2;
   S3 = M3*SQRT3f2o3;
   S4 = T2-T1;

   M = [ S1 S2; S3 S4];
end
