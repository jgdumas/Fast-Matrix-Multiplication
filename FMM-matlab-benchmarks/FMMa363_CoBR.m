function M = FMMa363_CoBR(A, nmin, peeling, level)
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


iM0 = FMMa363_CoBR( tA0, nmin, peeling, level);
iM1 = FMMa363_CoBR( tA1, nmin, peeling, level);
iM2 = FMMa363_CoBR( tA2, nmin, peeling, level);
iM3 = FMMa363_CoBR( tA3, nmin, peeling, level);
iM4 = FMMa363_CoBR( tA4, nmin, peeling, level);
iM5 = FMMa363_CoBR( tA5, nmin, peeling, level);
iM6 = FMMa363_CoBR( tA6, nmin, peeling, level);
iM7 = FMMa363_CoBR( tA7, nmin, peeling, level);
iM8 = FMMa363_CoBR( tA8, nmin, peeling, level);
iM9 = FMMa363_CoBR( tA9, nmin, peeling, level);
iM10 = FMMa363_CoBR( tA10, nmin, peeling, level);
iM11 = FMMa363_CoBR( tA11, nmin, peeling, level);
iM12 = FMMa363_CoBR( tA12, nmin, peeling, level);
iM13 = FMMa363_CoBR( tA13, nmin, peeling, level);
iM14 = FMMa363_CoBR( tA14, nmin, peeling, level);
iM15 = FMMa363_CoBR( tA15, nmin, peeling, level);
iM16 = FMMa363_CoBR( tA16, nmin, peeling, level);
iM17 = FMMa363_CoBR( tA17, nmin, peeling, level);
t18 = iM7+iM8;
t19 = iM12+iM13;
t20 = iM9-iM11;
t21 = iM4-iM5;
t22 = iM0-iM1;
t23 = iM12-iM13;
t24 = iM0+iM1;
t25 = iM15-iM17;
t26 = iM9+iM11;
t27 = iM7-iM8;
t28 = iM4+iM5;
t29 = iM15+iM17;
t31 = t18+t29;
r51 = iM10/4;
r52 = -iM14/4;
r53 = -t18/4;
r54 = -t20/4;
r55 = t27/4;
r56 = -(iM2+iM6)/4;
r57 = (iM3+iM14-iM16)/4;
r58 = (iM3+iM2-iM10)/4;
r59 = (t28+iM16+t19)/4;
r60 = (t19+t25)/4;
r61 = (iM16-t19)/4;
r62 = (iM2-t31)/4;
r63 = (t23+iM16-t21)/4;
r64 = (iM3-iM6)/4;
r65 = (iM10-t18)/4;
r66 = (t23-t25)/4;
r67 = (t20-t22)/4;
r68 = (t21-t26)/4;
r69 = (t24-t27)/4;
r70 = (t20+t24)/4;
r71 = (t22+t26)/4;
oM0 = r55-r58-r66;
oM1 = (t23+t31)/4-r58;
oM2 = t22/4+r59-r65;
oM3 = r57+r70-r53;
oM4 = r52-r62-r68;
oM5 = r51-r63-r69;
oM6 = r56-r61-r68;
oM7 = r59+r54-r56;
oM8 = r54+r56-r63;
oM9 = r53-r58-r60;
oM10 = r57-r71+r53;
oM11 = r61+r65+(t21-t24)/4;
oM12 = -r51-r59-r69;
oM13 = (t19-t29)/4-r64-r70;
oM14 = r60-r67+(iM3+iM6)/4;
oM15 = r64+r66-r71;
oM16 = r55-r57+r67;
oM17 = -r52-r54-r62-t28/4;

M = [ oM0 oM1 oM2 ; oM3 oM4 oM5 ; oM6 oM7 oM8 ; oM9 oM10 oM11 ; oM12 oM13 oM14 ; oM15 oM16 oM17 ] ;
end
