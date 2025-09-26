function M = DPS48a_CoBR(A, nmin, peeling, level)
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


iM0 = DPS48a_CoBR( tA0, nmin, peeling, level);
iM1 = DPS48a_CoBR( tA1, nmin, peeling, level);
iM2 = DPS48a_CoBR( tA2, nmin, peeling, level);
iM3 = DPS48a_CoBR( tA3, nmin, peeling, level);
iM4 = DPS48a_CoBR( tA4, nmin, peeling, level);
iM5 = DPS48a_CoBR( tA5, nmin, peeling, level);
iM6 = DPS48a_CoBR( tA6, nmin, peeling, level);
iM7 = DPS48a_CoBR( tA7, nmin, peeling, level);
iM8 = DPS48a_CoBR( tA8, nmin, peeling, level);
iM9 = DPS48a_CoBR( tA9, nmin, peeling, level);
iM10 = DPS48a_CoBR( tA10, nmin, peeling, level);
iM11 = DPS48a_CoBR( tA11, nmin, peeling, level);
iM12 = DPS48a_CoBR( tA12, nmin, peeling, level);
iM13 = DPS48a_CoBR( tA13, nmin, peeling, level);
iM14 = DPS48a_CoBR( tA14, nmin, peeling, level);
iM15 = DPS48a_CoBR( tA15, nmin, peeling, level);
t16 = iM0-iM13;
t17 = iM1+iM12;
t19 = iM9-t16;
t20 = iM10+iM11;
t21 = iM14+iM15;
t23 = iM14-iM15;
t24 = iM10-iM11;
t25 = iM2+iM3;
t26 = iM6-iM5*NTH2ROOT2once;
t27 = iM7-t19*NTH2ROOT2o2;
t28 = iM2-iM3;
r29 = (iM4-t17)*NTH2ROOT2o2;
t29 = iM6+r29;
t31 = iM7-t26;
t32 = t27+t29;
t33 = t21-iM8*NTH2ROOT2once;
r34 = (iM5-iM8)*NTH2ROOT2o2;
t34 = t28-r34;
t36 = iM7+t26;
t37 = t20-iM12*NTH2ROOT2once;
t38 = t29-t27;
r39 = (iM5+iM8)*NTH2ROOT2o2;
t39 = t25-r39;
r40 = iM1*NTH2ROOT2once;
t40 = t25-r40;
t41 = t23+r34;
r42 = iM9*NTH2ROOT2once;
r43 = iM13*NTH2ROOT2once;
r45 = (iM9+t16)*NTH2ROOT2o2;
oM0 = t23+t31-r43;
oM1 = t20+t39+r29-r45;
oM2 = t34-t38;
oM3 = t34-r29-r45-t24;
oM4 = t20+t40-r42;
oM5 = t24+t28+(iM0+iM8)*NTH2ROOT2once;
oM6 = t37-t33;
oM7 = r40-t28+t31;
oM8 = t36+t40;
oM9 = t33+t37;
oM10 = t32+t39;
oM11 = (iM4+t17+t19)*NTH2ROOT2o2-t24-t41;
oM12 = t24+r42-r43-t23;
oM13 = t38-t41;
oM14 = t21-r39-t32;
oM15 = t21+t36-r43;

M = [ oM0 oM1 oM2 oM3 ; oM4 oM5 oM6 oM7 ; oM8 oM9 oM10 oM11 ; oM12 oM13 oM14 oM15 ] ;
end
