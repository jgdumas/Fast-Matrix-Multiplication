function M = FMMa336_CoBR(A, nmin, peeling, level)
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


iM0 = FMMa336_CoBR( tA0, nmin, peeling, level);
iM1 = FMMa336_CoBR( tA1, nmin, peeling, level);
iM2 = FMMa336_CoBR( tA2, nmin, peeling, level);
iM3 = FMMa336_CoBR( tA3, nmin, peeling, level);
iM4 = FMMa336_CoBR( tA4, nmin, peeling, level);
iM5 = FMMa336_CoBR( tA5, nmin, peeling, level);
iM6 = FMMa336_CoBR( tA6, nmin, peeling, level);
iM7 = FMMa336_CoBR( tA7, nmin, peeling, level);
iM8 = FMMa336_CoBR( tA8, nmin, peeling, level);
iM9 = FMMa336_CoBR( tA9, nmin, peeling, level);
iM10 = FMMa336_CoBR( tA10, nmin, peeling, level);
iM11 = FMMa336_CoBR( tA11, nmin, peeling, level);
iM12 = FMMa336_CoBR( tA12, nmin, peeling, level);
iM13 = FMMa336_CoBR( tA13, nmin, peeling, level);
iM14 = FMMa336_CoBR( tA14, nmin, peeling, level);
iM15 = FMMa336_CoBR( tA15, nmin, peeling, level);
iM16 = FMMa336_CoBR( tA16, nmin, peeling, level);
iM17 = FMMa336_CoBR( tA17, nmin, peeling, level);
t18 = iM3-iM7;
t19 = iM6+iM17;
t20 = iM5-iM13;
t21 = iM10+iM15;
t22 = iM4+iM13;
t23 = iM8-iM14;
t25 = iM2-t19;
t26 = iM9+iM16;
t27 = iM11+iM2+iM12;
t28 = t18+t20;
t29 = iM16+iM17;
t30 = iM0+iM5;
t31 = iM12+iM15;
t32 = iM8+iM14;
t33 = iM0-iM11;
r59 = (-t18+t22)/4;
r60 = -(t25+t26)/4;
r61 = (t21+t27)/4;
r62 = (t23-iM9-iM10)/4;
r63 = (t23+t28)/4;
r64 = (iM0+iM7-t20)/4;
r65 = (iM15-iM3-iM6)/4;
r66 = (iM4-iM3-t30)/4;
r67 = (iM1-t21-t25)/4;
r68 = (iM13+t30-iM7)/4;
r69 = (t29-t31)/4;
r70 = (t26+iM2+t19)/4;
r71 = (t21-t27)/4;
r72 = (t29+t31)/4;
r73 = (t28+t32)/4;
r74 = (t18+t22)/4;
r75 = (iM16+iM1+t33)/4;
oM0 = r63+r72;
oM1 = r63-r72;
oM2 = r59-r61;
oM3 = r70-r64;
oM4 = t23/4-r65-r75;
oM5 = r60-r68;
oM6 = r61+r59;
oM7 = r73-r69;
oM8 = r64+r70;
oM9 = r69+r73;
oM10 = r62+(iM1-iM4+iM5+iM12+iM17)/4;
oM11 = t32/4-r65+r75;
oM12 = -r60-r68;
oM13 = r71+r74;
oM14 = r62+(iM6+iM7-t22-t33)/4;
oM15 = r66-r67;
oM16 = r66+r67;
oM17 = r71-r74;

M = [ oM0 oM1 oM2 oM3 oM4 oM5 ; oM6 oM7 oM8 oM9 oM10 oM11 ; oM12 oM13 oM14 oM15 oM16 oM17 ] ;
end
