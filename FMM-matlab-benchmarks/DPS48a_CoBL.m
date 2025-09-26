function M = DPS48a_CoBL(A, nmin, peeling, level)
NTH2ROOT2once=nthroot(2,2);
NTH2ROOT2o2=nthroot(2,2)/2;
NTH2ROOT2t2o=2/nthroot(2,2);
NTH2ROOT2t1o=1/nthroot(2,2);
NTH2ROOT2o4=nthroot(2,2)/4;
NTH2ROOT2o8=nthroot(2,2)/8;
NTH2ROOT2o16=nthroot(2,2)/16;
NTH2ROOT2o3=nthroot(2,2)/3;
NTH2ROOT2f2o3=nthroot(2,2)*2/3;

  if nargin < 2, nmin = 4; end    % Threshold to conventional
  if nargin < 3, peeling = 1; end % Static (1) or Dynamic (2) peeling
  if nargin < 4, level = 8; end   % Verbose level
[m,n] = size(A);
if (m<=nmin)||(n<=nmin)||(m<4)||(n<4)
   M=A;
else
[m,n] = size(A);
m0 = 0; m1 = 1*m/4; m2 = 2*m/4; m3 = 3*m/4; m4 = m;
 r0 = m0+1:m1; r1 = m1+1:m2; r2 = m2+1:m3; r3 = m3+1:m4;
n0 = 0; n1 = 1*n/4; n2 = 2*n/4; n3 = 3*n/4; n4 = n;
 c0 = n0+1:n1; c1 = n1+1:n2; c2 = n2+1:n3; c3 = n3+1:n4;
tA0 = A(r0,c0);
tA1 = A(r0,c1);
tA2 = A(r0,c2);
tA3 = A(r0,c3);
tA4 = A(r1,c0);
tA5 = A(r1,c1);
tA6 = A(r1,c2);
tA7 = A(r1,c3);
tA8 = A(r2,c0);
tA9 = A(r2,c1);
tA10 = A(r2,c2);
tA11 = A(r2,c3);
tA12 = A(r3,c0);
tA13 = A(r3,c1);
tA14 = A(r3,c2);
tA15 = A(r3,c3);


iM0 = DPS48a_CoBL( tA0, nmin, peeling, level);
iM1 = DPS48a_CoBL( tA1, nmin, peeling, level);
iM2 = DPS48a_CoBL( tA2, nmin, peeling, level);
iM3 = DPS48a_CoBL( tA3, nmin, peeling, level);
iM4 = DPS48a_CoBL( tA4, nmin, peeling, level);
iM5 = DPS48a_CoBL( tA5, nmin, peeling, level);
iM6 = DPS48a_CoBL( tA6, nmin, peeling, level);
iM7 = DPS48a_CoBL( tA7, nmin, peeling, level);
iM8 = DPS48a_CoBL( tA8, nmin, peeling, level);
iM9 = DPS48a_CoBL( tA9, nmin, peeling, level);
iM10 = DPS48a_CoBL( tA10, nmin, peeling, level);
iM11 = DPS48a_CoBL( tA11, nmin, peeling, level);
iM12 = DPS48a_CoBL( tA12, nmin, peeling, level);
iM13 = DPS48a_CoBL( tA13, nmin, peeling, level);
iM14 = DPS48a_CoBL( tA14, nmin, peeling, level);
iM15 = DPS48a_CoBL( tA15, nmin, peeling, level);
t16 = iM5+iM12;
t17 = iM7+iM14;
t18 = iM1+iM8;
t19 = iM2+iM4;
t20 = iM3+iM10;
t21 = iM11-iM13;
t22 = iM3-iM10;
t23 = iM0+iM9;
t24 = iM6+iM15;
t25 = t17-t18;
t26 = t16+t23;
t27 = t20+t24;
t28 = t19-t21;
t29 = iM5-iM12;
t30 = t22+t28;
t31 = t17+t18;
t32 = iM11+iM13;
t34 = t26+t27;
oM9 = t25+t22-t29;
t39 = iM1-iM8;
t41 = t19+t21;
t43 = t22-t28;
t45 = t26+iM6-iM15;
t46 = t16-t20;
t47 = t18+iM14-iM7;
t48 = iM9-iM0-t17;
t49 = iM4-iM2-t21-t25;
t50 = t27-t26;
t51 = t19+t31+t32;
oM0 = t16+t20+t31;
oM1 = t49-t50;
oM2 = t27-t29-t39+t41-t48;
oM3 = t46-t25;
oM4 = t20-t29-t41;
oM5 = t30-t16;
oM6 = t25+t46;
oM7 = t50+t49;
oM8 = t19+t23+t24-t32-oM9;
oM10 = t51-t34;
oM11 = t45+t47-t43;
oM12 = t16-t43;
oM13 = t34+t51;
oM14 = t30-t45+t47;
oM15 = t16+t24+t30+t39-t48;

M = [ oM0 oM1 oM2 oM3 ; oM4 oM5 oM6 oM7 ; oM8 oM9 oM10 oM11 ; oM12 oM13 oM14 oM15 ] ;
end
