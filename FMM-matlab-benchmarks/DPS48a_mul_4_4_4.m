function C = DPS48a_mul_4_4_4(A, B, nmin, peeling, level)
if nargin < 3, nmin = 4; end    % Threshold to conventional
if nargin < 4, peeling = 1; end % Static (1) or Dynamic (2) peeling
if nargin < 5, level = 8; end   % Verbose level
[m,k] = size(A); [k2,n] = size(B);
if (k2 ~= k), error('Incompatible matrix dimensions.'); end
% Recursively cuts into nmin*4^l x nmin*4^l x nmin*4^l blocks, with decreasing maximal l
if (m<=nmin)||(k<=nmin)||(n<=nmin)||(m<4)||(k<4)||(n<4)
  % fprintf("# MM Direct: %d x %d x %d\n",m,k,n)
  C = A*B;
else
  C = zeros(m,n);
  mu=m-rem(m,4);ku=k-rem(k,4);nu=n-rem(n,4);
  l=ceil(min([log(mu)/log(4),log(ku)/log(4),log(nu)/log(4)]));
  if (peeling == 1)
    l = min([floor(log(m/nmin)/log(4)),floor(log(k/nmin)/log(4)),floor(log(n/nmin)/log(4))]);
    mu=nmin*4^l; ku=nmin*4^l; nu=nmin*4^l;
  end
  if (mu < m) || (ku < k) || (nu < n)
    % fprintf("# Core SubMatrix[%d]: %d x %d x %d\n",l,mu,ku,nu)
    C(1:mu,1:nu)=DPS48a_mul_4_4_4(A(1:mu,1:ku),B(1:ku,1:nu),nmin, peeling, level);
    if (m > mu)
      % fprintf("# MM peel m: %d x %d x %d\n",m-mu,k,n)
      C(mu+1:m,1:n)=C(mu+1:m,1:n)+DPS48a_mul(A(mu+1:m,1:k),B,nmin, peeling, level);
    end
    if (k > ku) && (mu > 0) && (nu > 0)
      % fprintf("# MM peel k: %d x %d x %d\n",mu,k-ku,nu)
      C(1:mu,1:nu)=C(1:mu,1:nu)+DPS48a_mul(A(1:mu,ku+1:k),B(ku+1:k,1:nu),nmin, peeling, level);
    end
    if (n > nu) && (mu > 0)
      % fprintf("# MM peel n: %d x %d x %d\n",mu,k,n-nu)
      C(1:mu,nu+1:n)=C(1:mu,nu+1:n)+DPS48a_mul(A(1:mu,1:k),B(1:k,nu+1:n),nmin, peeling, level);
    end
  else
    if l>=level, fprintf("# Core<4;4;4>[%d]: %d x %d x %d\n",l,m,k,n); end

[m,n] = size(A);
m0 = 0; m1 = 1*m/4; m2 = 2*m/4; m3 = 3*m/4; m4 = m;
 r0 = m0+1:m1; r1 = m1+1:m2; r2 = m2+1:m3; r3 = m3+1:m4;
n0 = 0; n1 = 1*n/4; n2 = 2*n/4; n3 = 3*n/4; n4 = n;
 c0 = n0+1:n1; c1 = n1+1:n2; c2 = n2+1:n3; c3 = n3+1:n4;
t16 = A(r1,c3)-A(r2,c3);
t17 = A(r2,c2)-A(r1,c3);
t18 = A(r1,c2)+A(r2,c1);
t19 = A(r2,c0)+A(r2,c2);
t20 = A(r1,c1)+A(r2,c1);
t21 = A(r1,c1)-A(r0,c0);
t22 = A(r0,c0)-A(r1,c2);
t23 = A(r2,c0)+A(r2,c3);
t24 = A(r3,c0)+t16;
t27 = A(r2,c1)-A(r1,c1);
t28 = A(r2,c0)-A(r2,c3);
t30 = A(r1,c2)-A(r2,c1);
t32 = A(r1,c3)+A(r2,c3);
t33 = A(r2,c2)-A(r2,c0);
t34 = A(r1,c3)+A(r2,c2);
t35 = A(r0,c0)+A(r1,c1);
t36 = A(r0,c0)+A(r1,c2);
r37 = A(r2,c1)*2;
r38 = (t18-t21)/2;
r39 = (t19+t16)/2;
r40 = (t20+t22)/2;
r41 = (t17-t23)/2;
oA0 = A(r3,c1)-t20+t28;
oA3 = A(r0,c3)+t17+t22;
oA5 = A(r0,c2)-r37;
oA6 = t16-t35-A(r0,c1);
oA8 = t18+t32-A(r0,c2);
oA9 = A(r0,c0)*2+A(r3,c2);
oA11 = t21-t32-A(r3,c2);
oA12 = r37-t19-t24;
oA13 = t33-A(r0,c2)-t30;
oA15 = t24+t30;
oA16 = t17-t36-A(r3,c3);
oA17 = A(r3,c1)+t22-t23;
oA18 = A(r0,c3)+A(r2,c2)*2;
oA19 = t27-t34-A(r3,c3);
oA20 = A(r2,c0)*2-A(r1,c0);
oA21 = A(r1,c3)*2+A(r3,c3);
oA22 = t20-t34-A(r0,c3);
oA23 = A(r0,c1)+t19+t21;
oA24 = t18-t24;
oA26 = A(r3,c2)+t33+t35;
oA27 = A(r2,c3)*2-A(r3,c1);
oA28 = A(r0,c1)+A(r1,c1)*2;
oA29 = A(r1,c0)+t27-t28;
oA30 = t23+t36-A(r1,c0);
oA33 = r38-A(r2,c2);
oA34 = A(r1,c1)+r41;
oA37 = A(r2,c1)-r39;
oA38 = r39-A(r0,c0);
oA39 = r40-A(r1,c3);
oA40 = r40-A(r2,c0);
oA41 = r41-A(r1,c2);
oA43 = A(r2,c3)+r38;
oA1 = A(r3,c3);
oA2 = A(r0,c1);
oA4 = A(r3,c0);
oA7 = A(r3,c2);
oA10 = A(r1,c0);
oA14 = A(r3,c1);
oA25 = A(r0,c2);
oA31 = A(r0,c3);
oA32 = A(r2,c0);
oA35 = A(r1,c1);
oA36 = A(r2,c2);
oA42 = A(r2,c3);
oA44 = A(r1,c2);
oA45 = A(r1,c3);
oA46 = A(r0,c0);
oA47 = A(r2,c1);

[m,n] = size(B);
m0 = 0; m1 = 1*m/4; m2 = 2*m/4; m3 = 3*m/4; m4 = m;
 r0 = m0+1:m1; r1 = m1+1:m2; r2 = m2+1:m3; r3 = m3+1:m4;
n0 = 0; n1 = 1*n/4; n2 = 2*n/4; n3 = 3*n/4; n4 = n;
 c0 = n0+1:n1; c1 = n1+1:n2; c2 = n2+1:n3; c3 = n3+1:n4;
oB34 = B(r2,c3)-B(r0,c1)-B(r2,c0);
v48 = B(r0,c2)-B(r0,c1);
v49 = B(r2,c0)+oB34;
v50 = B(r1,c0)+B(r1,c1);
oB36 = v49+B(r0,c1);
g52 = v48/2;
v52 = B(r2,c0)-g52;
v53 = B(r0,c3)-B(r1,c0);
g54 = B(r0,c2)/2;
v54 = B(r3,c1)-g54;
oB47 = B(r0,c3)+v50;
g56 = oB34/2;
v56 = -g56-v52;
v57 = B(r1,c3)-B(r1,c1);
g58 = B(r0,c3)/2;
v58 = B(r0,c0)+g58;
oB39 = B(r1,c3)+v50;
g60 = B(r2,c0)/2;
v60 = B(r1,c2)+g60;
g61 = oB36/2;
v61 = B(r3,c0)+g61;
oB35 = v53-B(r1,c3);
g63 = B(r1,c0)/2;
v63 = B(r3,c2)+g63;
v64 = B(r2,c1)-g60;
v65 = B(r2,c2)-g58;
oB45 = B(r0,c2)-v49;
g67 = -B(r1,c3)/2;
v67 = B(r3,c3)+g67;
g69 = B(r0,c1)/2;
g70 = B(r1,c1)/2;
g73 = oB47/2;
g74 = -v57/2;
oB16 = g74-B(r0,c0);
oB3 = B(r3,c2)-g70+g67;
oB32 = B(r2,c0)-v48;
oB19 = B(r2,c2)-B(r0,c3)-g74;
oB28 = B(r3,c3)-v53/2;
oB23 = v56-v58;
oB26 = v65-g56;
oB13 = g69-v58;
oB24 = v61-g70;
oB25 = v61-g73;
oB33 = B(r0,c3)-v57;
oB22 = v54+oB39/2;
oB11 = v54+g56;
oB1 = v67+oB45/2;
oB38 = v56*2;
oB2 = v60+oB35/2;
oB8 = v63+g69;
oB29 = v54-g63;
oB10 = B(r1,c2)+g52;
oB12 = g54+v65;
oB31 = g61-v67;
oB15 = v64-g70;
oB6 = v63-v56;
oB20 = B(r2,c1)-v52;
oB5 = v64+g73;
oB21 = B(r3,c0)+v49-g52;
oB9 = g67-v60;
oB0 = B(r2,c2);
oB4 = B(r3,c1);
oB7 = B(r3,c3);
oB14 = B(r1,c2);
oB17 = B(r0,c0);
oB18 = B(r3,c0);
oB27 = B(r2,c1);
oB30 = B(r3,c2);
oB37 = B(r0,c2);
oB40 = B(r0,c3);
oB41 = B(r0,c1);
oB42 = B(r2,c0);
oB43 = B(r1,c0);
oB44 = B(r1,c1);
oB46 = B(r1,c3);

iC0 = DPS48a_mul( oA0, oB0, nmin, peeling, level);
iC1 = DPS48a_mul( oA1, oB1, nmin, peeling, level);
iC2 = DPS48a_mul( oA2, oB2, nmin, peeling, level);
iC3 = DPS48a_mul( oA3, oB3, nmin, peeling, level);
iC4 = DPS48a_mul( oA4, oB4, nmin, peeling, level);
iC5 = DPS48a_mul( oA5, oB5, nmin, peeling, level);
iC6 = DPS48a_mul( oA6, oB6, nmin, peeling, level);
iC7 = DPS48a_mul( oA7, oB7, nmin, peeling, level);
iC8 = DPS48a_mul( oA8, oB8, nmin, peeling, level);
iC9 = DPS48a_mul( oA9, oB9, nmin, peeling, level);
iC10 = DPS48a_mul( oA10, oB10, nmin, peeling, level);
iC11 = DPS48a_mul( oA11, oB11, nmin, peeling, level);
iC12 = DPS48a_mul( oA12, oB12, nmin, peeling, level);
iC13 = DPS48a_mul( oA13, oB13, nmin, peeling, level);
iC14 = DPS48a_mul( oA14, oB14, nmin, peeling, level);
iC15 = DPS48a_mul( oA15, oB15, nmin, peeling, level);
iC16 = DPS48a_mul( oA16, oB16, nmin, peeling, level);
iC17 = DPS48a_mul( oA17, oB17, nmin, peeling, level);
iC18 = DPS48a_mul( oA18, oB18, nmin, peeling, level);
iC19 = DPS48a_mul( oA19, oB19, nmin, peeling, level);
iC20 = DPS48a_mul( oA20, oB20, nmin, peeling, level);
iC21 = DPS48a_mul( oA21, oB21, nmin, peeling, level);
iC22 = DPS48a_mul( oA22, oB22, nmin, peeling, level);
iC23 = DPS48a_mul( oA23, oB23, nmin, peeling, level);
iC24 = DPS48a_mul( oA24, oB24, nmin, peeling, level);
iC25 = DPS48a_mul( oA25, oB25, nmin, peeling, level);
iC26 = DPS48a_mul( oA26, oB26, nmin, peeling, level);
iC27 = DPS48a_mul( oA27, oB27, nmin, peeling, level);
iC28 = DPS48a_mul( oA28, oB28, nmin, peeling, level);
iC29 = DPS48a_mul( oA29, oB29, nmin, peeling, level);
iC30 = DPS48a_mul( oA30, oB30, nmin, peeling, level);
iC31 = DPS48a_mul( oA31, oB31, nmin, peeling, level);
iC32 = DPS48a_mul( oA32, oB32, nmin, peeling, level);
iC33 = DPS48a_mul( oA33, oB33, nmin, peeling, level);
iC34 = DPS48a_mul( oA34, oB34, nmin, peeling, level);
iC35 = DPS48a_mul( oA35, oB35, nmin, peeling, level);
iC36 = DPS48a_mul( oA36, oB36, nmin, peeling, level);
iC37 = DPS48a_mul( oA37, oB37, nmin, peeling, level);
iC38 = DPS48a_mul( oA38, oB38, nmin, peeling, level);
iC39 = DPS48a_mul( oA39, oB39, nmin, peeling, level);
iC40 = DPS48a_mul( oA40, oB40, nmin, peeling, level);
iC41 = DPS48a_mul( oA41, oB41, nmin, peeling, level);
iC42 = DPS48a_mul( oA42, oB42, nmin, peeling, level);
iC43 = DPS48a_mul( oA43, oB43, nmin, peeling, level);
iC44 = DPS48a_mul( oA44, oB44, nmin, peeling, level);
iC45 = DPS48a_mul( oA45, oB45, nmin, peeling, level);
iC46 = DPS48a_mul( oA46, oB46, nmin, peeling, level);
iC47 = DPS48a_mul( oA47, oB47, nmin, peeling, level);

b5 = iC44+iC24+iC15;
b7 = iC40+iC17;
b8 = iC8-iC4;
b10 = iC28+iC25;
b11 = iC25+iC5+iC47;
v64 = iC3-iC16;
b12 = iC10+iC32;
v65 = iC17-iC30;
b16 = iC10-iC31-iC13;
b18 = iC26+iC23;
b19 = iC2-iC5;
b21 = iC13-iC8-iC41;
b22 = iC27-iC42;
b23 = iC21+iC11-iC27;
b31 = b19-iC3-iC39;
v58 = iC20+b12;
b36 = b11-b22;
b41 = b21+b7;
b42 = iC38-iC23+b23;
v48 = b31-b42;
oC10 = iC46+iC15*2+iC44-iC14+b22+v65+v64+v58;
oC5 = iC45+iC21-b8+b5-b36;
oC12 = iC43+iC30+b10+b41+b42;
oC2 = iC36+iC31-iC20*2+b18-b12-b11-b5;
oC6 = iC35+iC14+iC28+iC2+v58+b36;
oC3 = iC33+iC16-b21-v48;
oC11 = iC29+iC0-b10;
oC9 = iC22-iC19-b19;
oC0 = iC18+iC20-b18;
oC7 = iC12+b16;
oC13 = iC9-iC15-v64;
oC14 = iC7+iC24-v65;
oC1 = iC6-b23;
oC15 = iC1+b8-iC14;
oC8 = iC37+iC19-iC0-iC4-b7+b16+v48+oC3;
oC4 = iC34+iC19+iC0+iC26-iC11+b31+b41+oC3;

C = [ oC0 oC1 oC2 oC3 ; oC4 oC5 oC6 oC7 ; oC8 oC9 oC10 oC11 ; oC12 oC13 oC14 oC15 ] ;
  end
end
end
