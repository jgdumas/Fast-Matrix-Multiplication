function M = FMMa633_ICoB(A, nmin, peeling, level)
  if nargin < 3, nmin = 4; end    % Threshold to conventional
  if nargin < 4, peeling = 1; end % Static (1) or Dynamic (2) peeling
  if nargin < 5, level = 8; end   % Verbose level
[m,n] = size(A);
if (m <= nmin)||(n <= nmin)
   M=A;
else
[m,n] = size(A);
m0 = 0; m1 = 1*m/6; m2 = 2*m/6; m3 = 3*m/6; m4 = 4*m/6; m5 = 5*m/6; m6 = m;
 r0 = m0+1:m1; r1 = m1+1:m2; r2 = m2+1:m3; r3 = m3+1:m4; r4 = m4+1:m5; r5 = m5+1:m6;
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
tA9 = A(r3,c0);
tA10 = A(r3,c1);
tA11 = A(r3,c2);
tA12 = A(r4,c0);
tA13 = A(r4,c1);
tA14 = A(r4,c2);
tA15 = A(r5,c0);
tA16 = A(r5,c1);
tA17 = A(r5,c2);


iM0 = FMMa633_ICoB( tA0, nmin, peeling, level);
iM1 = FMMa633_ICoB( tA1, nmin, peeling, level);
iM2 = FMMa633_ICoB( tA2, nmin, peeling, level);
iM3 = FMMa633_ICoB( tA3, nmin, peeling, level);
iM4 = FMMa633_ICoB( tA4, nmin, peeling, level);
iM5 = FMMa633_ICoB( tA5, nmin, peeling, level);
iM6 = FMMa633_ICoB( tA6, nmin, peeling, level);
iM7 = FMMa633_ICoB( tA7, nmin, peeling, level);
iM8 = FMMa633_ICoB( tA8, nmin, peeling, level);
iM9 = FMMa633_ICoB( tA9, nmin, peeling, level);
iM10 = FMMa633_ICoB( tA10, nmin, peeling, level);
iM11 = FMMa633_ICoB( tA11, nmin, peeling, level);
iM12 = FMMa633_ICoB( tA12, nmin, peeling, level);
iM13 = FMMa633_ICoB( tA13, nmin, peeling, level);
iM14 = FMMa633_ICoB( tA14, nmin, peeling, level);
iM15 = FMMa633_ICoB( tA15, nmin, peeling, level);
iM16 = FMMa633_ICoB( tA16, nmin, peeling, level);
iM17 = FMMa633_ICoB( tA17, nmin, peeling, level);
t73 = iM17/2;
t50 = (iM13+iM6)/2;
t49 = (iM4-iM5)/2;
t48 = (-iM15-iM3)/2;
t47 = (iM8-iM3)/2;
t46 = (iM15+iM6)/2;
t45 = (iM9-iM14)/2;
t44 = (iM16-iM17)/2;
t43 = (-iM5-iM4)/2;
t42 = (iM0-iM2)/2;
t41 = (iM7-iM14)/2;
t40 = (iM9-iM1)/2;
t37 = (iM11+iM0)/2;
t36 = (iM10-iM13)/2;
t35 = (iM12+iM10-iM8)/2;
t30 = (-iM13-iM11-iM3)/2;
t28 = (iM7-iM16+iM2-iM1)/2;
b21 = (iM6-iM4)/2-t47;
b23 = iM2/2+t44;
b29 = t44-iM10/2+t41;
b30 = t43-t48+t42+t41;
b31 = t50-iM16/2-t40;
b32 = (-iM12-iM8)/2+t42+t40;
t34 = t73-t49;
b34 = iM12/2+t46;
t29 = t45-t37;
t26 = t45+t37;
b35 = (iM15-iM17+iM12-iM5)/2-t36;
oM15 = (iM17+iM11)/2+t35;
oM4 = (iM14+iM9+iM0)/2+t49+t35;
b37 = (iM7+iM1)/2+t46-t30;
oM14 = t49+t46+t30;
oM1 = t34+t29;
t25 = -t28;
oM9 = t34+t28;
b38 = t47+t36-t26;
oM0 = t43-t73+t26;
oM7 = t25-t29;
oM2 = t50+t48-t35+t25;
oM8 = b21-b35;
oM6 = b21+b35;
oM13 = b29-b32;
oM12 = b29+b32;
oM5 = b37+b23;
oM3 = b37-b23;
oM11 = b38+b34;
oM10 = b38-b34;
oM17 = b31-b30;
oM16 = b31+b30;

M = [ oM0 oM1 oM2 ; oM3 oM4 oM5 ; oM6 oM7 oM8 ; oM9 oM10 oM11 ; oM12 oM13 oM14 ; oM15 oM16 oM17 ] ;
end
