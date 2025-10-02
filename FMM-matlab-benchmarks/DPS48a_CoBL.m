function M = DPS48a_CoBL(A, nmin, peeling, level)
  if nargin < 2, nmin = 4; end    % Threshold to conventional
  if nargin < 3, peeling = 1; end % Static (1) or Dynamic (2) peeling
  if nargin < 4, level = 8; end   % Verbose level
[m,n] = size(A);
if (m<=nmin)||(n<=nmin)||(m<4)||(n<4)
   M=A;
else
[m,n] = size(A);
m0 = 0; m1 = 1*m/4; m2 = 2*m/4; m3 = 3*m/4; m4 = m;
 r0 = m0+1:m1; r1 = m1+1:m2; r2 = m2+1:m3; r3 = m3+1:m4;
n0 = 0; n1 = 1*n/4; n2 = 2*n/4; n3 = 3*n/4; n4 = n;
 c0 = n0+1:n1; c1 = n1+1:n2; c2 = n2+1:n3; c3 = n3+1:n4;
tA0 = A(r0,c0);
tA1 = A(r0,c1);
tA2 = A(r0,c2);
tA3 = A(r0,c3);
tA4 = A(r1,c0);
tA5 = A(r1,c1);
tA6 = A(r1,c2);
tA7 = A(r1,c3);
tA8 = A(r2,c0);
tA9 = A(r2,c1);
tA10 = A(r2,c2);
tA11 = A(r2,c3);
tA12 = A(r3,c0);
tA13 = A(r3,c1);
tA14 = A(r3,c2);
tA15 = A(r3,c3);


iM0 = DPS48a_CoBL( tA0, nmin, peeling, level);
iM1 = DPS48a_CoBL( tA1, nmin, peeling, level);
iM2 = DPS48a_CoBL( tA2, nmin, peeling, level);
iM3 = DPS48a_CoBL( tA3, nmin, peeling, level);
iM4 = DPS48a_CoBL( tA4, nmin, peeling, level);
iM5 = DPS48a_CoBL( tA5, nmin, peeling, level);
iM6 = DPS48a_CoBL( tA6, nmin, peeling, level);
iM7 = DPS48a_CoBL( tA7, nmin, peeling, level);
iM8 = DPS48a_CoBL( tA8, nmin, peeling, level);
iM9 = DPS48a_CoBL( tA9, nmin, peeling, level);
iM10 = DPS48a_CoBL( tA10, nmin, peeling, level);
iM11 = DPS48a_CoBL( tA11, nmin, peeling, level);
iM12 = DPS48a_CoBL( tA12, nmin, peeling, level);
iM13 = DPS48a_CoBL( tA13, nmin, peeling, level);
iM14 = DPS48a_CoBL( tA14, nmin, peeling, level);
iM15 = DPS48a_CoBL( tA15, nmin, peeling, level);

t16 = iM4+iM11;
t17 = iM1+iM14;
t18 = iM7-iM13;
t19 = iM2+iM8;
t20 = iM6+iM9;
t21 = iM3+iM12;
t22 = t16+t20;
t23 = iM5-iM15;
t24 = iM0+iM10;
t25 = iM1-iM14;
t26 = t18+t19;
t27 = t21-t17;
t28 = t17+t21;
t29 = t16-t20;
t30 = t25+t26;
t31 = t18-t19;
t32 = iM4-iM11;
t35 = iM12-iM3-t22;
t37 = t23+iM10-iM0;
t38 = t25-t26;
oM11 = t27-t29;
t40 = t32-t31;
oM10 = t22+t28;
t43 = t24-t23;
t45 = t18+iM8-iM2;
t46 = iM5+t24+iM15;
oM8 = t29+t27;
oM7 = t22-t28;
t49 = iM6-iM9;
t50 = iM13+iM7+t19;
t51 = t23+t24;
oM0 = t16+t38;
oM1 = t43+t49-t27-t40;
oM2 = t30-t35+t51;
oM3 = t46+t50-oM10;
oM4 = t37-t45+oM8;
oM5 = t40-t17;
oM6 = t17-t31-t32;
oM9 = t16+t30;
oM12 = t21+t30-t43+t49-t16;
oM13 = t50+oM11-t46;
oM14 = t35-t38+t51;
oM15 = t37+t45-oM7;

M = [ oM0 oM1 oM2 oM3 ; oM4 oM5 oM6 oM7 ; oM8 oM9 oM10 oM11 ; oM12 oM13 oM14 oM15 ] ;
end
