function M = FMM336_ICoB(A, nmin, peeling, level)
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
n0 = 0; n1 = 1*n/6; n2 = 2*n/6; n3 = 3*n/6; n4 = 4*n/6; n5 = 5*n/6; n6 = n;
 c0 = n0+1:n1; c1 = n1+1:n2; c2 = n2+1:n3; c3 = n3+1:n4; c4 = n4+1:n5; c5 = n5+1:n6;
tA0 = A(r0,c0);
tA1 = A(r0,c1);
tA2 = A(r0,c2);
tA3 = A(r0,c3);
tA4 = A(r0,c4);
tA5 = A(r0,c5);
tA6 = A(r1,c0);
tA7 = A(r1,c1);
tA8 = A(r1,c2);
tA9 = A(r1,c3);
tA10 = A(r1,c4);
tA11 = A(r1,c5);
tA12 = A(r2,c0);
tA13 = A(r2,c1);
tA14 = A(r2,c2);
tA15 = A(r2,c3);
tA16 = A(r2,c4);
tA17 = A(r2,c5);


iM0 = FMM336_ICoB( tA0, nmin, peeling, level);
iM1 = FMM336_ICoB( tA1, nmin, peeling, level);
iM2 = FMM336_ICoB( tA2, nmin, peeling, level);
iM3 = FMM336_ICoB( tA3, nmin, peeling, level);
iM4 = FMM336_ICoB( tA4, nmin, peeling, level);
iM5 = FMM336_ICoB( tA5, nmin, peeling, level);
iM6 = FMM336_ICoB( tA6, nmin, peeling, level);
iM7 = FMM336_ICoB( tA7, nmin, peeling, level);
iM8 = FMM336_ICoB( tA8, nmin, peeling, level);
iM9 = FMM336_ICoB( tA9, nmin, peeling, level);
iM10 = FMM336_ICoB( tA10, nmin, peeling, level);
iM11 = FMM336_ICoB( tA11, nmin, peeling, level);
iM12 = FMM336_ICoB( tA12, nmin, peeling, level);
iM13 = FMM336_ICoB( tA13, nmin, peeling, level);
iM14 = FMM336_ICoB( tA14, nmin, peeling, level);
iM15 = FMM336_ICoB( tA15, nmin, peeling, level);
iM16 = FMM336_ICoB( tA16, nmin, peeling, level);
iM17 = FMM336_ICoB( tA17, nmin, peeling, level);
t48 = iM12+iM11;
b1 = iM11-iM16;
t40 = iM17+iM9;
t46 = iM14+iM8;
t47 = iM7+iM6;
t39 = iM17-iM6;
b8 = iM3-iM15-iM9-iM4;
b9 = iM6-iM10-iM4-iM3;
t44 = iM13+iM2;
b11 = iM4+iM2;
t45 = iM16-iM1;
b13 = iM12-iM5+iM1;
b14 = iM7+iM0;
b15 = iM13-iM10+iM0;
t32 = iM14-t48;
b17 = t47-iM5;
b18 = iM10+iM0+t45;
b21 = iM5+t48+t44;
b22 = iM9-iM8-b11;
b23 = b1-iM15;
t26 = (iM8-iM14+b14)/8;
t23 = (iM12-iM3-iM1-t40)/8;
t22 = (-iM10-t44-t39)/8;
b28 = b23-iM14-iM3;
b29 = iM13-b14+t40+b23;
oM12 = (iM16+iM11-t46+b15)/8;
b34 = iM8-iM17+t45+b17;
t20 = (iM15-b11+b17)/8;
oM8 = iM17+t47+b15+b13;
oM5 = t46-iM2+b13;
oM3 = b8-b15-iM2;
oM7 = iM7+t39+t46-b8;
oM16 = (-b8-b1-b13)/8;
oM14 = b22-b28;
oM2 = b22+b28;
oM10 = t26-t22;
oM4 = t26+t22;
oM13 = b21+b18;
oM1 = b21-b18;
oM6 = t20-t23;
oM0 = t23+t20;
oM15 = t32-b34;
oM9 = t32+b34;
oM17 = b29-b9;
oM11 = b9+b29;

M = [ oM0 oM1 oM2 oM3 oM4 oM5 ; oM6 oM7 oM8 oM9 oM10 oM11 ; oM12 oM13 oM14 oM15 oM16 oM17 ] ;
end
