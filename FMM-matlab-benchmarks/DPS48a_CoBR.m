function M = DPS48a_CoBR(A, nmin, peeling, level)
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


iM0 = DPS48a_CoBR( tA0, nmin, peeling, level);
iM1 = DPS48a_CoBR( tA1, nmin, peeling, level);
iM2 = DPS48a_CoBR( tA2, nmin, peeling, level);
iM3 = DPS48a_CoBR( tA3, nmin, peeling, level);
iM4 = DPS48a_CoBR( tA4, nmin, peeling, level);
iM5 = DPS48a_CoBR( tA5, nmin, peeling, level);
iM6 = DPS48a_CoBR( tA6, nmin, peeling, level);
iM7 = DPS48a_CoBR( tA7, nmin, peeling, level);
iM8 = DPS48a_CoBR( tA8, nmin, peeling, level);
iM9 = DPS48a_CoBR( tA9, nmin, peeling, level);
iM10 = DPS48a_CoBR( tA10, nmin, peeling, level);
iM11 = DPS48a_CoBR( tA11, nmin, peeling, level);
iM12 = DPS48a_CoBR( tA12, nmin, peeling, level);
iM13 = DPS48a_CoBR( tA13, nmin, peeling, level);
iM14 = DPS48a_CoBR( tA14, nmin, peeling, level);
iM15 = DPS48a_CoBR( tA15, nmin, peeling, level);

t16 = iM12-iM4;
t17 = iM0+iM5;
t18 = iM13+t17;
t19 = iM1+t16;
t20 = iM8+iM9;
t21 = iM14-iM6;
t22 = iM7-iM15;
t23 = iM8-iM9;
t24 = iM4+iM12;
t25 = iM2-iM10;
t26 = iM7+iM15;
t27 = t18-t19;
t28 = t17-iM13;
t29 = iM7+iM11;
t30 = iM5-iM13;
t31 = t19-t23;
t32 = iM2-iM14;
t34 = iM1+t18-t24;
oM0 = t21+t16;
oM1 = t18+t19+t23+(iM3+iM15)*2;
oM2 = t20+t27+t29*2;
oM3 = t28-t31-t21*2;
oM4 = t20-t25*2-t34;
oM5 = t18+t31-(iM6+iM10)*2;
oM6 = t16-t22;
oM7 = t27+t32*2-t20;
oM8 = t20+t22*2-t28-t19;
oM9 = t22-t30;
oM10 = t30-t21;
oM11 = t20+t26*2+t34;
oM12 = t24-t26;
oM13 = iM4+iM8+t29;
oM14 = iM0-iM8+t25;
oM15 = iM0-iM12+t32;

M = [ oM0 oM1 oM2 oM3 ; oM4 oM5 oM6 oM7 ; oM8 oM9 oM10 oM11 ; oM12 oM13 oM14 oM15 ] ;
end
