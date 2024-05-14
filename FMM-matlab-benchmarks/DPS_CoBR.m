function M = DPS_CoBR(A, nmin)
%          Right change of basis of the sparsification
%          of DPS's algorithm via alternative basis.

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
   SQRT3f2o3=sqrt(3)*2/3;

   M1 = DPS_CoBR(A(i,i), nmin);
   M2 = DPS_CoBR(A(i,j), nmin);
   M3 = DPS_CoBR(A(j,i), nmin);
   M4 = DPS_CoBR(A(j,j), nmin);

T1=M2*SQRT3o3;
T2=M1+M4;
T3=M2-M3;

S1=M2*SQRT3f2o3;
S2=M1-T1;
S3=T1-M4;
S4=T3*SQRT3o2-T2/2;

   M = [ S1 S2; S3 S4];
end
