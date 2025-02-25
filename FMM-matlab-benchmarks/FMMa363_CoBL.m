function M = FMMa363_CoBL(A, nmin, peeling, level)
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


iM0 = FMMa363_CoBL( tA0, nmin, peeling, level);
iM1 = FMMa363_CoBL( tA1, nmin, peeling, level);
iM2 = FMMa363_CoBL( tA2, nmin, peeling, level);
iM3 = FMMa363_CoBL( tA3, nmin, peeling, level);
iM4 = FMMa363_CoBL( tA4, nmin, peeling, level);
iM5 = FMMa363_CoBL( tA5, nmin, peeling, level);
iM6 = FMMa363_CoBL( tA6, nmin, peeling, level);
iM7 = FMMa363_CoBL( tA7, nmin, peeling, level);
iM8 = FMMa363_CoBL( tA8, nmin, peeling, level);
iM9 = FMMa363_CoBL( tA9, nmin, peeling, level);
iM10 = FMMa363_CoBL( tA10, nmin, peeling, level);
iM11 = FMMa363_CoBL( tA11, nmin, peeling, level);
iM12 = FMMa363_CoBL( tA12, nmin, peeling, level);
iM13 = FMMa363_CoBL( tA13, nmin, peeling, level);
iM14 = FMMa363_CoBL( tA14, nmin, peeling, level);
iM15 = FMMa363_CoBL( tA15, nmin, peeling, level);
iM16 = FMMa363_CoBL( tA16, nmin, peeling, level);
iM17 = FMMa363_CoBL( tA17, nmin, peeling, level);
t18 = iM1-iM9;
t19 = iM2+iM11;
t20 = iM6-iM13;
t21 = iM10-iM14;
t22 = iM13-iM17;
t23 = iM0-iM15;
t24 = iM8-iM16;
t25 = iM5-t20;
t27 = iM14-iM15;
t28 = iM1-iM4;
t29 = t18+t19-iM12;
t30 = iM2-iM9;
t32 = iM5+iM7-iM12;
t33 = iM8+iM10;
t34 = iM4-iM15;
r59 = iM17/2;
r60 = (t22-t27)/2;
r61 = (t18+t23)/2;
r62 = (iM11+t28)/2;
r63 = (iM12+t33-iM3)/2;
r64 = (t22+t27)/2;
r65 = (-iM0-t19+iM17)/2;
r66 = (iM16+t29)/2;
r67 = (iM16-t29)/2;
r68 = (t24+t25)/2;
r69 = -(iM2+iM0+t28)/2;
r70 = (iM6+iM14+iM3-iM16-iM7)/2;
r71 = (t20+iM5+t24)/2;
r72 = (-t30+t32)/2;
r73 = (t18-t23)/2;
r74 = (t25+t21-iM3)/2;
r75 = (t21+t34)/2;
oM0 = r73-r71;
oM1 = r71+r73;
oM2 = r65-r70;
oM3 = (t21-t30-t32-t34)/2;
oM4 = r70+r65;
oM5 = r75+r72;
oM6 = r69-r74;
oM7 = r60+r66;
oM8 = r64-r67;
oM9 = r59+(iM6+iM7-iM9+iM11-t23+t33-iM4)/2;
oM10 = t22/2-r62+r63;
oM11 = r61+r68;
oM12 = r64+r67;
oM13 = r60-r66;
oM14 = r72-r75;
oM15 = -r69-r74;
oM16 = iM13/2+r59+r62+r63;
oM17 = r68-r61;

M = [ oM0 oM1 oM2 oM3 oM4 oM5 ; oM6 oM7 oM8 oM9 oM10 oM11 ; oM12 oM13 oM14 oM15 oM16 oM17 ] ;
end
