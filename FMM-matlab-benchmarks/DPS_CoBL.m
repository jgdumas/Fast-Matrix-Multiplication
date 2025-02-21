function M = DPS_CoBL(A, nmin)
SQRT3o2=sqrt(3)/2;
SQRT3o3=sqrt(3)/3;
SQRT3f2o3=sqrt(3)*2/3;

if nargin < 2, nmin = 8; end     % Threshold to conventional
[m,n] = size(A);
if (m <= nmin)||(n <= nmin)
   M=A;
else
[m,n] = size(A);
m0 = 0; m1 = 1*m/2; m2 = m;
 r0 = m0+1:m1; r1 = m1+1:m2;
n0 = 0; n1 = 1*n/2; n2 = n;
 c0 = n0+1:n1; c1 = n1+1:n2;
tA0 = A(r0,c0);
tA1 = A(r0,c1);
tA2 = A(r1,c0);
tA3 = A(r1,c1);


iM0 = DPS_CoBL( tA0, nmin);
iM1 = DPS_CoBL( tA1, nmin);
iM2 = DPS_CoBL( tA2, nmin);
iM3 = DPS_CoBL( tA3, nmin);
r4 = iM3*SQRT3o3;
oM0 = (iM2-iM1)/2-(iM0+iM3)*SQRT3o2;
oM1 = iM1+r4;
oM2 = iM2-r4;
oM3 = iM3*SQRT3f2o3;

M = [ oM0 oM1 ; oM2 oM3 ] ;
end
