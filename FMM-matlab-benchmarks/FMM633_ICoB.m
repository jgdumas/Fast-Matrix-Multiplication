function M = FMM633_ICoB(A, nmin, peeling, level)
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


iM0 = FMM633_ICoB( tA0, nmin, peeling, level);
iM1 = FMM633_ICoB( tA1, nmin, peeling, level);
iM2 = FMM633_ICoB( tA2, nmin, peeling, level);
iM3 = FMM633_ICoB( tA3, nmin, peeling, level);
iM4 = FMM633_ICoB( tA4, nmin, peeling, level);
iM5 = FMM633_ICoB( tA5, nmin, peeling, level);
iM6 = FMM633_ICoB( tA6, nmin, peeling, level);
iM7 = FMM633_ICoB( tA7, nmin, peeling, level);
iM8 = FMM633_ICoB( tA8, nmin, peeling, level);
iM9 = FMM633_ICoB( tA9, nmin, peeling, level);
iM10 = FMM633_ICoB( tA10, nmin, peeling, level);
iM11 = FMM633_ICoB( tA11, nmin, peeling, level);
iM12 = FMM633_ICoB( tA12, nmin, peeling, level);
iM13 = FMM633_ICoB( tA13, nmin, peeling, level);
iM14 = FMM633_ICoB( tA14, nmin, peeling, level);
iM15 = FMM633_ICoB( tA15, nmin, peeling, level);
iM16 = FMM633_ICoB( tA16, nmin, peeling, level);
iM17 = FMM633_ICoB( tA17, nmin, peeling, level);
b1 = -iM15-iM14;
b2 = iM14-iM15;
t48 = iM13-iM10;
b3 = -iM13-iM10;
b5 = -iM11-iM8;
b6 = iM8-iM11;
t58 = iM7/8;
t45 = iM17+iM7;
t47 = iM16+iM6;
b8 = iM5-iM12;
t41 = iM12+iM5;
t57 = iM4/8;
t53 = iM3+iM2;
b10 = iM2-iM3;
t52 = iM9-iM0;
b13 = iM4-iM17-iM0;
b14 = t53/8-t58;
t25 = t52+b10;
t27 = b1-t48;
b15 = b6-iM4;
t33 = iM1+b3;
t37 = (iM17-b3)/8;
b17 = b8-iM1;
t38 = (iM1+b8)/8;
t26 = t41+b2;
t28 = t53+b5;
t24 = (iM6-iM16-b5)/8;
b19 = b17-iM16;
b20 = b1+t48+b17;
oM6 = (t41-iM4-iM1)/8+b14;
t31 = iM9+b13;
b22 = b13-iM9;
oM16 = (-t47-t31)/8;
t22 = iM6-b19;
oM14 = t47-iM7+t28;
t21 = t58+t28/8;
oM2 = b15+t27;
t23 = t27/8;
oM12 = t33+t31+t26;
b25 = t25-t45;
b26 = (t52+b2)/8+t24;
oM9 = t57-t38+t24;
oM10 = (iM9-iM17+iM0)/8+b14-t23;
oM3 = t47/8+t38+t23;
oM1 = t45+t25+t22;
oM11 = t57-(b5+t22)/8;
oM15 = (t26-t33)/8-t21;
oM4 = t21-b22/8;
oM0 = iM6+b19+b25;
oM5 = (b25-b15)/8;
oM13 = b22+b20;
oM17 = t58+(b10-b6-b20)/8;
oM8 = b26-t37;
oM7 = t37+b26;

M = [ oM0 oM1 oM2 ; oM3 oM4 oM5 ; oM6 oM7 oM8 ; oM9 oM10 oM11 ; oM12 oM13 oM14 ; oM15 oM16 oM17 ] ;
end
