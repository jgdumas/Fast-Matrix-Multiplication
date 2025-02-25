function M = FMMa633_CoBL(A, nmin, peeling, level)
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


iM0 = FMMa633_CoBL( tA0, nmin, peeling, level);
iM1 = FMMa633_CoBL( tA1, nmin, peeling, level);
iM2 = FMMa633_CoBL( tA2, nmin, peeling, level);
iM3 = FMMa633_CoBL( tA3, nmin, peeling, level);
iM4 = FMMa633_CoBL( tA4, nmin, peeling, level);
iM5 = FMMa633_CoBL( tA5, nmin, peeling, level);
iM6 = FMMa633_CoBL( tA6, nmin, peeling, level);
iM7 = FMMa633_CoBL( tA7, nmin, peeling, level);
iM8 = FMMa633_CoBL( tA8, nmin, peeling, level);
iM9 = FMMa633_CoBL( tA9, nmin, peeling, level);
iM10 = FMMa633_CoBL( tA10, nmin, peeling, level);
iM11 = FMMa633_CoBL( tA11, nmin, peeling, level);
iM12 = FMMa633_CoBL( tA12, nmin, peeling, level);
iM13 = FMMa633_CoBL( tA13, nmin, peeling, level);
iM14 = FMMa633_CoBL( tA14, nmin, peeling, level);
iM15 = FMMa633_CoBL( tA15, nmin, peeling, level);
iM16 = FMMa633_CoBL( tA16, nmin, peeling, level);
iM17 = FMMa633_CoBL( tA17, nmin, peeling, level);
t18 = iM13-iM17;
t19 = iM7+iM15;
t20 = iM4-iM9;
t21 = iM8+iM12;
t22 = iM2-t19;
t23 = iM3+t18;
t24 = iM1-iM11;
t25 = iM5+iM11;
t26 = iM0-iM9;
t27 = t23-t24;
t29 = iM8+iM17;
t30 = iM12+iM15;
t31 = iM0+iM5;
t33 = iM7+iM16;
t34 = t20-t25;
t35 = iM6+iM2+iM16;
r56 = iM13/4;
r57 = t22/4;
r58 = t23/4;
r59 = -t26/4;
r60 = (iM10-t21)/4;
r61 = (iM6-t27)/4;
r62 = (iM14-t29)/4;
r63 = (t20-t35)/4;
r64 = (iM10+t21)/4;
r65 = (iM8-t24)/4;
r66 = (iM17+iM1+t31)/4;
r67 = (t26-t30)/4;
r68 = (iM2+t19)/4;
r69 = (iM4+iM14-iM15+iM6+iM10)/4;
r70 = (iM11+iM13-iM5-iM12)/4;
r71 = (t22+t34)/4;
r72 = (iM3+iM14-t33)/4;
r73 = (iM3-t18)/4;
oM0 = r64-r57-r58;
oM1 = (iM6+t27+t30)/4-r59;
oM2 = r61+r67;
oM3 = r66+r69;
oM4 = r69-r66;
oM5 = r60-r73-r57;
oM6 = r72-(iM0+iM9)/4-r65;
oM7 = r59-r65-r72;
oM8 = r61-r67;
oM9 = r56+(iM12-t34-t35)/4;
oM10 = r63+r70;
oM11 = r68+r73-r64;
oM12 = r60-r68-r58;
oM13 = (-iM14-t20-t25-t29)/4-r57;
oM14 = r63-r70;
oM15 = r62-r71;
oM16 = -r62-r71;
oM17 = r56-r60+(iM4+t31+t33-iM1)/4;

M = [ oM0 oM1 oM2 ; oM3 oM4 oM5 ; oM6 oM7 oM8 ; oM9 oM10 oM11 ; oM12 oM13 oM14 ; oM15 oM16 oM17 ] ;
end
