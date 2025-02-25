function M = FMM336_CoBR(A, nmin, peeling, level)
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


iM0 = FMM336_CoBR( tA0, nmin, peeling, level);
iM1 = FMM336_CoBR( tA1, nmin, peeling, level);
iM2 = FMM336_CoBR( tA2, nmin, peeling, level);
iM3 = FMM336_CoBR( tA3, nmin, peeling, level);
iM4 = FMM336_CoBR( tA4, nmin, peeling, level);
iM5 = FMM336_CoBR( tA5, nmin, peeling, level);
iM6 = FMM336_CoBR( tA6, nmin, peeling, level);
iM7 = FMM336_CoBR( tA7, nmin, peeling, level);
iM8 = FMM336_CoBR( tA8, nmin, peeling, level);
iM9 = FMM336_CoBR( tA9, nmin, peeling, level);
iM10 = FMM336_CoBR( tA10, nmin, peeling, level);
iM11 = FMM336_CoBR( tA11, nmin, peeling, level);
iM12 = FMM336_CoBR( tA12, nmin, peeling, level);
iM13 = FMM336_CoBR( tA13, nmin, peeling, level);
iM14 = FMM336_CoBR( tA14, nmin, peeling, level);
iM15 = FMM336_CoBR( tA15, nmin, peeling, level);
iM16 = FMM336_CoBR( tA16, nmin, peeling, level);
iM17 = FMM336_CoBR( tA17, nmin, peeling, level);
t18 = iM8+iM14;
t20 = -iM4-iM10;
t21 = iM0-iM6;
t24 = iM4-iM10;
r26 = iM1/8;
t27 = iM0+iM6;
r28 = iM11/8;
t28 = iM16+r26-r28;
t30 = iM5-iM17;
r33 = t18/8;
t33 = t21+r33;
r35 = (iM3-iM15)/8;
r36 = iM9/8;
t36 = iM12+r26-r36;
r37 = (iM7-iM13)/8;
t37 = r37-iM12-r35;
t38 = r37+t20+r28;
r39 = (iM5+iM17)/8;
r40 = iM2/8;
t40 = iM12+r40;
t43 = t27+r33;
t44 = t21-r39-r40;
t45 = iM16+r36;
r46 = (iM8-iM14)/8;
t46 = t21-r46;
r47 = (iM7+iM13)/8;
t47 = r28-t20+r47;
t48 = t24+r28;
r55 = (iM3+iM15)/8;
r56 = (t18+t30)/8;
oM0 = r26+t20-t44-r55;
oM1 = t45+r40+r37-r39-t27;
oM2 = t38-t46-r36;
oM3 = r35-t28-t46;
oM4 = t30/8-t36+r46-t24;
oM5 = t33+t47-r36;
oM6 = t37-t48-r40;
oM7 = t20-t36-r56;
oM8 = r39+r46-iM16-t37;
oM9 = r36+t38-t43;
oM10 = t40+t47-r35;
oM11 = t28+t43+r35;
oM12 = iM16-t37+r56;
oM13 = t45-r47-t44;
oM14 = r33+t24-t36+r39;
oM15 = t28-t33-r55;
oM16 = r55-t38-t40;
oM17 = t33+t48+r36-r37;

M = [ oM0 oM1 oM2 oM3 oM4 oM5 ; oM6 oM7 oM8 oM9 oM10 oM11 ; oM12 oM13 oM14 oM15 oM16 oM17 ] ;
end
