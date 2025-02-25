function M = DPS_ICoB(A, nmin, peeling, level)
SQRT3o2=sqrt(3)/2;
SQRT3o3=sqrt(3)/3;
SQRT3f2o3=sqrt(3)*2/3;

  if nargin < 3, nmin = 4; end    % Threshold to conventional
  if nargin < 4, peeling = 1; end % Static (1) or Dynamic (2) peeling
  if nargin < 5, level = 8; end   % Verbose level
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


iM0 = DPS_ICoB( tA0, nmin, peeling, level);
iM1 = DPS_ICoB( tA1, nmin, peeling, level);
iM2 = DPS_ICoB( tA2, nmin, peeling, level);
iM3 = DPS_ICoB( tA3, nmin, peeling, level);

oM3 = iM2*SQRT3o2;
t5 = iM2/2;
oM2 = t5-iM1;
oM1 = -t5-iM0;
oM0 = oM3-iM3*SQRT3f2o3+(iM0-iM1)*SQRT3o3;

M = [ oM0 oM1 ; oM2 oM3 ] ;
end
