function M = FMM633_CoBR(A, nmin, peeling, level)
  if nargin < 3, nmin = 4; end    % Threshold to conventional
  if nargin < 4, peeling = 1; end % Static (1) or Dynamic (2) peeling
  if nargin < 5, level = 8; end   % Verbose level
[m,n] = size(A);
if (m <= nmin)||(n <= nmin)
   M=A;
else
[m,n] = size(A);
m0 = 0; m1 = 1*m/3; m2 = 2*m/3; m3 = m;
 r0 = m0+1:m1; r1 = m1+1:m2; r2 = m2+1:m3;
n0 = 0; n1 = 1*n/3; n2 = 2*n/3; n3 = n;
 c0 = n0+1:n1; c1 = n1+1:n2; c2 = n2+1:n3;
tA0 = A(r0,c0);
tA1 = A(r0,c1);
tA2 = A(r0,c2);
tA3 = A(r1,c0);
tA4 = A(r1,c1);
tA5 = A(r1,c2);
tA6 = A(r2,c0);
tA7 = A(r2,c1);
tA8 = A(r2,c2);


iM0 = FMM633_CoBR( tA0, nmin, peeling, level);
iM1 = FMM633_CoBR( tA1, nmin, peeling, level);
iM2 = FMM633_CoBR( tA2, nmin, peeling, level);
iM3 = FMM633_CoBR( tA3, nmin, peeling, level);
iM4 = FMM633_CoBR( tA4, nmin, peeling, level);
iM5 = FMM633_CoBR( tA5, nmin, peeling, level);
iM6 = FMM633_CoBR( tA6, nmin, peeling, level);
iM7 = FMM633_CoBR( tA7, nmin, peeling, level);
iM8 = FMM633_CoBR( tA8, nmin, peeling, level);
t9 = iM7+iM6;
t10 = t9-iM2;
t11 = -iM4-iM5;
t12 = iM4-iM5;
t13 = iM0-t12;
t15 = -t10-iM0;
t16 = iM6-iM0-iM7;
oM0 = iM0-t10+t12;
oM1 = iM2+t9+t13;
oM2 = iM2-t11-t16;
oM3 = t11-iM7-iM8;
oM4 = -t11-t15;
oM5 = t10-t13;
oM6 = iM1+t16;
oM7 = t15-t11;
oM8 = t9-iM3-iM4;

M = [ oM0 oM1 oM2 ; oM3 oM4 oM5 ; oM6 oM7 oM8 ] ;
end
