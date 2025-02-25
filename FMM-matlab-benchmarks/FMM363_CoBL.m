function M = FMM363_CoBL(A, nmin, peeling, level)
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


iM0 = FMM363_CoBL( tA0, nmin, peeling, level);
iM1 = FMM363_CoBL( tA1, nmin, peeling, level);
iM2 = FMM363_CoBL( tA2, nmin, peeling, level);
iM3 = FMM363_CoBL( tA3, nmin, peeling, level);
iM4 = FMM363_CoBL( tA4, nmin, peeling, level);
iM5 = FMM363_CoBL( tA5, nmin, peeling, level);
iM6 = FMM363_CoBL( tA6, nmin, peeling, level);
iM7 = FMM363_CoBL( tA7, nmin, peeling, level);
iM8 = FMM363_CoBL( tA8, nmin, peeling, level);
iM9 = FMM363_CoBL( tA9, nmin, peeling, level);
iM10 = FMM363_CoBL( tA10, nmin, peeling, level);
iM11 = FMM363_CoBL( tA11, nmin, peeling, level);
iM12 = FMM363_CoBL( tA12, nmin, peeling, level);
iM13 = FMM363_CoBL( tA13, nmin, peeling, level);
iM14 = FMM363_CoBL( tA14, nmin, peeling, level);
iM15 = FMM363_CoBL( tA15, nmin, peeling, level);
iM16 = FMM363_CoBL( tA16, nmin, peeling, level);
iM17 = FMM363_CoBL( tA17, nmin, peeling, level);
t18 = iM2-iM14;
t20 = iM1+iM13;
t21 = iM9-iM15;
t23 = iM2+iM14;
t24 = iM9+iM15;
t25 = -iM11-iM17;
r26 = iM12/8;
t27 = r26-iM7-iM5;
r28 = (iM4+iM10)/8;
r29 = (iM0-iM6)/8;
t29 = iM5+r29;
t31 = t18-t21;
t32 = iM11-iM17;
r33 = iM16/8;
t34 = iM8-r33;
t36 = iM3-iM7-r33;
t37 = iM1-iM13;
t38 = iM3+t20+r28;
t39 = t24-t23;
t40 = t24-r28;
t41 = iM5+t34;
t42 = -t25-r29;
r43 = (iM0+iM6)/8;
t43 = t18-r43;
t44 = t32-r26;
t45 = t20+t29;
t47 = t37-t25-r26;
r50 = (iM4-iM10)/8;
oM0 = t37+t41-r29-t21;
oM1 = t38+t44-iM8;
oM2 = t36-r29-t18-t32;
oM3 = r33-t31-t47;
oM4 = iM7-iM8+t40+t42;
oM5 = t23-t36+t42;
oM6 = t23-t29-t38;
oM7 = t20-t21-r43-t41;
oM8 = iM3-t18-t45+r50;
oM9 = t31-r50-t27;
oM10 = t34+t45-t24;
oM11 = t40-t18-t27;
oM12 = t21+r28-t23-t27;
oM13 = t38+t43-iM5;
oM14 = t20+t31+t44+r33;
oM15 = t39-t47-r33;
oM16 = t27+t39-r50;
oM17 = t25-t36-t43;

M = [ oM0 oM1 oM2 oM3 oM4 oM5 ; oM6 oM7 oM8 oM9 oM10 oM11 ; oM12 oM13 oM14 oM15 oM16 oM17 ] ;
end
