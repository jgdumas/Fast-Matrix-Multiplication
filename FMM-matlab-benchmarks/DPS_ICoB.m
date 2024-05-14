function M = DPS_ICoB(A, nmin)
%          Inverse change of basis of the sparsification
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

   M1 = DPS_ICoB(A(i,i), nmin);
   M2 = DPS_ICoB(A(i,j), nmin);
   M3 = DPS_ICoB(A(j,i), nmin);
   M4 = DPS_ICoB(A(j,j), nmin);

T1=M4/2;
T2=M2-M3;

S4=M4*SQRT3o2;
S1=S4+T2*SQRT3o3-M1*SQRT3f2o3;
S2=-M2-T1;
S3=T1-M3;

   M = [ S1 S2; S3 S4];
end
