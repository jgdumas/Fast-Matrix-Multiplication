function M = FMM363_CoBR(A, nmin, peeling, level)
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


iM0 = FMM363_CoBR( tA0, nmin, peeling, level);
iM1 = FMM363_CoBR( tA1, nmin, peeling, level);
iM2 = FMM363_CoBR( tA2, nmin, peeling, level);
iM3 = FMM363_CoBR( tA3, nmin, peeling, level);
iM4 = FMM363_CoBR( tA4, nmin, peeling, level);
iM5 = FMM363_CoBR( tA5, nmin, peeling, level);
iM6 = FMM363_CoBR( tA6, nmin, peeling, level);
iM7 = FMM363_CoBR( tA7, nmin, peeling, level);
iM8 = FMM363_CoBR( tA8, nmin, peeling, level);
iM9 = FMM363_CoBR( tA9, nmin, peeling, level);
iM10 = FMM363_CoBR( tA10, nmin, peeling, level);
iM11 = FMM363_CoBR( tA11, nmin, peeling, level);
iM12 = FMM363_CoBR( tA12, nmin, peeling, level);
iM13 = FMM363_CoBR( tA13, nmin, peeling, level);
iM14 = FMM363_CoBR( tA14, nmin, peeling, level);
iM15 = FMM363_CoBR( tA15, nmin, peeling, level);
iM16 = FMM363_CoBR( tA16, nmin, peeling, level);
iM17 = FMM363_CoBR( tA17, nmin, peeling, level);
t18 = iM0-iM1;
t19 = iM12+iM13;
t20 = iM7-iM8;
t24 = iM15+iM17;
t26 = iM0+iM1;
t27 = iM12-iM13;
r29 = iM6/8;
r30 = (iM9+iM11)/8;
t30 = t19+r30;
r31 = (iM15-iM17)/8;
t31 = t26+r29-r31;
r32 = iM10/8;
r33 = iM16/8;
r34 = iM3/8;
t34 = iM2-r32+r34;
r35 = t20/8;
t35 = t18-r33-r35;
t36 = iM14+r32;
r38 = (iM4+iM5)/8;
t39 = t30-r29;
r40 = (iM7+iM8)/8;
t42 = iM2+r33;
t43 = iM2-r38+r30;
t44 = t18+r40+r33;
t45 = t36+r29;
t46 = t19+r38;
t47 = iM14+r34;
r52 = t24/8;
r53 = (iM4-iM5)/8;
r55 = (iM9-iM11)/8;
r56 = (t20+t24)/8;
oM0 = r32-r38-t27-t35;
oM1 = t43+r31+r40-iM14;
oM2 = t34+r56-t19;
oM3 = t39-t42-r53;
oM4 = iM14+t43+r56;
oM5 = t45-r52+r53-t26;
oM6 = t47-r30-t44;
oM7 = t31-r34-t30;
oM8 = r55-t35-t47;
oM9 = r53-t19-t35-r32;
oM10 = r53-t18-t45-r31;
oM11 = r34+r52-t18-t39;
oM12 = t44+t46-r32;
oM13 = t27+t31+r34+r55;
oM14 = r29+t42+t46-r55;
oM15 = t27-t34+r52+r40;
oM16 = t31+t36+r38;
oM17 = r35-t27-t34+r31;

M = [ oM0 oM1 oM2 ; oM3 oM4 oM5 ; oM6 oM7 oM8 ; oM9 oM10 oM11 ; oM12 oM13 oM14 ; oM15 oM16 oM17 ] ;
end
