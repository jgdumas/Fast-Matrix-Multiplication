function M = DPS48a_ICoB(A, nmin, peeling, level)
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


iM0 = DPS48a_ICoB( tA0, nmin, peeling, level);
iM1 = DPS48a_ICoB( tA1, nmin, peeling, level);
iM2 = DPS48a_ICoB( tA2, nmin, peeling, level);
iM3 = DPS48a_ICoB( tA3, nmin, peeling, level);
iM4 = DPS48a_ICoB( tA4, nmin, peeling, level);
iM5 = DPS48a_ICoB( tA5, nmin, peeling, level);
iM6 = DPS48a_ICoB( tA6, nmin, peeling, level);
iM7 = DPS48a_ICoB( tA7, nmin, peeling, level);
iM8 = DPS48a_ICoB( tA8, nmin, peeling, level);
iM9 = DPS48a_ICoB( tA9, nmin, peeling, level);
iM10 = DPS48a_ICoB( tA10, nmin, peeling, level);
iM11 = DPS48a_ICoB( tA11, nmin, peeling, level);
iM12 = DPS48a_ICoB( tA12, nmin, peeling, level);
iM13 = DPS48a_ICoB( tA13, nmin, peeling, level);
iM14 = DPS48a_ICoB( tA14, nmin, peeling, level);
iM15 = DPS48a_ICoB( tA15, nmin, peeling, level);
t61 = iM14/8;
t60 = iM12/8;
t59 = iM7*NTH2ROOT2o8;
t58 = iM5/8;
b2 = iM14*NTH2ROOT2o16+iM5*NTH2ROOT2o8;
t39 = (iM4-iM9)*NTH2ROOT2o8;
t38 = (iM10+iM8)/8;
t37 = (iM15-iM0)/8;
t36 = (iM11+iM9)/8;
t35 = (iM11+iM8)*NTH2ROOT2o8;
t34 = (iM7-iM4)/8;
t33 = (iM10+iM9)/8;
b15 = iM10*NTH2ROOT2o8+(iM12+iM2)*NTH2ROOT2o16;
b16 = (-iM8-iM4)*NTH2ROOT2o8+(iM1-iM15)*NTH2ROOT2o16;
b19 = iM3*NTH2ROOT2o16+t36*NTH2ROOT2once;
b20 = (iM3-iM15-iM1)*NTH2ROOT2o16+t35;
t31 = t35/NTH2ROOT2once;
b21 = t36-t58+t34;
b23 = (iM6+iM0)*NTH2ROOT2o16+b2;
b24 = (iM13+iM0)*NTH2ROOT2o16-b2;
t29 = t58+t31;
b25 = t58+t34-t31;
b26 = (-iM10-iM7)/8-t39/NTH2ROOT2once;
b29 = b15-iM13*NTH2ROOT2o16-t59;
b30 = iM6*NTH2ROOT2o16+t59+b15;
oM5 = b23-b16;
oM4 = b23+b16;
oM10 = (iM3-iM13+iM2)/8+t29-b26;
oM6 = (-iM6-iM1)/8+t29+b26;
oM1 = b24-b20;
oM0 = b24+b20;
oM9 = b29-b19;
oM8 = b19+b29;
b31 = b21-t38;
oM11 = t60+t38+b21;
oM13 = b30-t39;
oM12 = t39+b30;
t16 = t33-b25;
oM7 = t61+t37+t33+b25;
oM15 = t60-b31;
oM2 = (iM13+iM3+iM1)/8+b31;
oM14 = (iM2-iM6)/8+t16;
oM3 = t37-t61+t16;

M = [ oM0 oM1 oM2 oM3 ; oM4 oM5 oM6 oM7 ; oM8 oM9 oM10 oM11 ; oM12 oM13 oM14 oM15 ] ;
end
