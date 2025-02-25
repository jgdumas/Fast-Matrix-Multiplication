function M = FMMa336_ICoB(A, nmin, peeling, level)
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


iM0 = FMMa336_ICoB( tA0, nmin, peeling, level);
iM1 = FMMa336_ICoB( tA1, nmin, peeling, level);
iM2 = FMMa336_ICoB( tA2, nmin, peeling, level);
iM3 = FMMa336_ICoB( tA3, nmin, peeling, level);
iM4 = FMMa336_ICoB( tA4, nmin, peeling, level);
iM5 = FMMa336_ICoB( tA5, nmin, peeling, level);
iM6 = FMMa336_ICoB( tA6, nmin, peeling, level);
iM7 = FMMa336_ICoB( tA7, nmin, peeling, level);
iM8 = FMMa336_ICoB( tA8, nmin, peeling, level);
iM9 = FMMa336_ICoB( tA9, nmin, peeling, level);
iM10 = FMMa336_ICoB( tA10, nmin, peeling, level);
iM11 = FMMa336_ICoB( tA11, nmin, peeling, level);
iM12 = FMMa336_ICoB( tA12, nmin, peeling, level);
iM13 = FMMa336_ICoB( tA13, nmin, peeling, level);
iM14 = FMMa336_ICoB( tA14, nmin, peeling, level);
iM15 = FMMa336_ICoB( tA15, nmin, peeling, level);
iM16 = FMMa336_ICoB( tA16, nmin, peeling, level);
iM17 = FMMa336_ICoB( tA17, nmin, peeling, level);
t71 = iM16/2;
t70 = iM4/2;
t50 = (-iM11-iM3)/2;
t49 = (iM10-iM7)/2;
t48 = (iM14-iM2)/2;
t47 = (-iM12-iM0)/2;
t46 = (iM11-iM3)/2;
t45 = (iM10+iM7)/2;
t44 = (iM15-iM12)/2;
t41 = (-iM17-iM9)/2;
t40 = (iM2-iM1)/2;
t38 = (iM15+iM8-iM1)/2;
t34 = (iM17+iM14+iM5)/2;
t33 = (-iM14-iM8-iM0)/2;
t32 = (iM13+iM9-iM8)/2;
t29 = (iM13+iM6+iM5)/2;
t28 = (iM15-iM16+iM6)/2;
b19 = (iM0-iM12)/2+t49;
t42 = -t48;
b20 = t46-t70;
t36 = t70+t46;
t30 = t71-t45;
b21 = t45-t47;
b22 = (iM5-iM17-iM1)/2-t44;
b24 = (iM6-iM4+iM0)/2+t40;
t35 = -t38;
oM5 = t47+t36+t34;
b26 = t49+t40+t32;
b27 = t50+t44+t32;
oM7 = t38+t36;
oM12 = (iM9-iM13+iM6+iM2)/2-t36+t30;
oM8 = t34+t30;
t25 = t41-t29;
t24 = t41+t29;
t23 = t71+b21;
oM3 = t35-t71+b21;
oM4 = t71+b20+b19;
t18 = t42-t25;
oM16 = t42-t35+t25;
oM15 = t48-b20-t24;
oM13 = t48+t24+t23;
oM10 = t23-t50-t70;
oM6 = t33-b22;
oM0 = t33+b22;
oM17 = t28-b26;
oM11 = t28+b26;
oM14 = b27-b24;
oM2 = b24+b27;
oM9 = t70-t50-t18;
oM1 = b19-t71+t18;

M = [ oM0 oM1 oM2 oM3 oM4 oM5 ; oM6 oM7 oM8 oM9 oM10 oM11 ; oM12 oM13 oM14 oM15 oM16 oM17 ] ;
end
