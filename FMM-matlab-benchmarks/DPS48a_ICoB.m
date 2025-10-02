function M = DPS48a_ICoB(A, nmin, peeling, level)
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


iM0 = DPS48a_ICoB( tA0, nmin, peeling, level);
iM1 = DPS48a_ICoB( tA1, nmin, peeling, level);
iM2 = DPS48a_ICoB( tA2, nmin, peeling, level);
iM3 = DPS48a_ICoB( tA3, nmin, peeling, level);
iM4 = DPS48a_ICoB( tA4, nmin, peeling, level);
iM5 = DPS48a_ICoB( tA5, nmin, peeling, level);
iM6 = DPS48a_ICoB( tA6, nmin, peeling, level);
iM7 = DPS48a_ICoB( tA7, nmin, peeling, level);
iM8 = DPS48a_ICoB( tA8, nmin, peeling, level);
iM9 = DPS48a_ICoB( tA9, nmin, peeling, level);
iM10 = DPS48a_ICoB( tA10, nmin, peeling, level);
iM11 = DPS48a_ICoB( tA11, nmin, peeling, level);
iM12 = DPS48a_ICoB( tA12, nmin, peeling, level);
iM13 = DPS48a_ICoB( tA13, nmin, peeling, level);
iM14 = DPS48a_ICoB( tA14, nmin, peeling, level);
iM15 = DPS48a_ICoB( tA15, nmin, peeling, level);

t60 = iM12/8;
t59 = iM10/8;
t58 = iM6/8;
t57 = iM5/8;
t40 = (iM2+iM0)/8;
t39 = (iM3+iM1)/8;
t38 = (iM12+iM1)/8;
t37 = (iM8-iM7)/8;
t36 = (iM13+iM10)/8;
t34 = (iM15+iM11)/8;
t33 = (iM11-iM15)/8;
t31 = (iM3-iM9)/8;
t29 = (iM14+iM4)/8;
b16 = iM9/8+t37;
b17 = iM14/8+t36;
b20 = t57-t37-t33;
b21 = t58+t34+t38;
b22 = t31-iM4/8;
t27 = t58+t40;
b25 = t29-t34-t57;
t19 = t39-b16;
oM6 = t59+t38+t40+b16;
oM15 = t31+b17-(iM8+iM0+iM5)/8;
t18 = t40+b17;
b27 = t60+t33+b17;
oM7 = (iM6-iM7+iM5+iM1-iM0)/8+b22;
oM14 = (iM9-iM13-iM12)/8+t29+t27;
oM2 = b25-(iM13+iM3+iM6)/8;
oM1 = t36-t60+b25;
oM10 = t39-t59+b20;
oM9 = t58-t38+b20;
oM5 = t27-t19;
oM4 = t27+t19;
oM13 = t18-b22;
oM12 = b22+t18;
oM11 = (iM7-iM4+iM2)/8+b21;
oM8 = t57+t37+b21;
oM3 = (iM8+iM2)/8+b27;
oM0 = b27-(iM4+iM5)/8;

M = [ oM0 oM1 oM2 oM3 ; oM4 oM5 oM6 oM7 ; oM8 oM9 oM10 oM11 ; oM12 oM13 oM14 oM15 ] ;
end
