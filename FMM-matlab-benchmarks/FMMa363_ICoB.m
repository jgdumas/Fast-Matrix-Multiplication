function M = FMMa363_ICoB(A, nmin, peeling, level)
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


iM0 = FMMa363_ICoB( tA0, nmin, peeling, level);
iM1 = FMMa363_ICoB( tA1, nmin, peeling, level);
iM2 = FMMa363_ICoB( tA2, nmin, peeling, level);
iM3 = FMMa363_ICoB( tA3, nmin, peeling, level);
iM4 = FMMa363_ICoB( tA4, nmin, peeling, level);
iM5 = FMMa363_ICoB( tA5, nmin, peeling, level);
iM6 = FMMa363_ICoB( tA6, nmin, peeling, level);
iM7 = FMMa363_ICoB( tA7, nmin, peeling, level);
iM8 = FMMa363_ICoB( tA8, nmin, peeling, level);
t14 = iM6-iM5;
t15 = iM5-iM3;
b2 = iM8-iM7-iM6-iM1;
b4 = iM8+iM3+iM0;
t13 = -t15;
b6 = iM2+iM1+b4;
oM4 = iM6+iM5-iM1+b4;
b7 = iM4-iM0+b2;
oM5 = iM0+t15+b2;
oM7 = b6-t14;
oM6 = t14-iM4+b6;
oM2 = t13-b7;
oM0 = iM2+t13+b7;
oM1 = iM2;
oM3 = iM7;
oM8 = iM4;

M = [ oM0 oM1 oM2 ; oM3 oM4 oM5 ; oM6 oM7 oM8 ] ;
end
