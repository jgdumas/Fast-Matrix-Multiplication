function C = FMMa_3_6_3(A, B, nmin, peeling, level)
if nargin < 3, nmin = 3; end    % Threshold to conventional
if nargin < 4, peeling = 1; end % Static (1) or Dynamic (2) peeling
if nargin < 5, level = 3; end   % Verbose level
[m,k] = size(A); [k2,n] = size(B);
if (k2 ~= k), error('Incompatible matrix dimensions.'); end
% Recursively cuts into nmin*3^l x nmin*6^l x nmin*3^l blocks, with decreasing maximal l
if (m<=nmin)||(k<=nmin)||(n<=nmin)||(m<3)||(k<6)||(n<3)
  % fprintf("# MM Direct: %d x %d x %d\n",m,k,n)
  C = A*B;
else
  C = zeros(m,n);
  mu=m-rem(m,3);ku=k-rem(k,6);nu=n-rem(n,3);
  l=ceil(min([log(mu)/log(3),log(ku)/log(6),log(nu)/log(3)]));
  if (peeling == 1)
    l = min([floor(log(m/nmin)/log(3)),floor(log(k/nmin)/log(6)),floor(log(n/nmin)/log(3))]);
    mu=nmin*2^l; ku=nmin*2^l; nu=nmin*2^l;
  end
  if (mu < m) || (ku < k) || (nu < n)
    % fprintf("# Core SubMatrix[%d]: %d x %d x %d\n",l,mu,ku,nu)
    C(1:mu,1:nu)=FMMa_3_6_3(A(1:mu,1:ku),B(1:ku,1:nu),nmin, peeling, level);
    if (m > mu)
      % fprintf("# MM peel m: %d x %d x %d\n",m-mu,k,n)
      C(mu+1:m,1:n)=C(mu+1:m,1:n)+FMMa(A(mu+1:m,1:k),B,nmin, peeling, level);
    end
    if (k > ku) && (mu > 0) && (nu > 0)
      % fprintf("# MM peel k: %d x %d x %d\n",mu,k-ku,nu)
      C(1:mu,1:nu)=C(1:mu,1:nu)+FMMa(A(1:mu,ku+1:k),B(ku+1:k,1:nu),nmin, peeling, level);
    end
    if (n > nu) && (mu > 0)
      % fprintf("# MM peel n: %d x %d x %d\n",mu,k,n-nu)
      C(1:mu,nu+1:n)=C(1:mu,nu+1:n)+FMMa(A(1:mu,1:k),B(1:k,nu+1:n),nmin, peeling, level);
    end
  else
    if l>=level, fprintf("# Core<3;6;3>[%d]: %d x %d x %d\n",l,m,k,n); end

[m,n] = size(A);
m0 = 0; m1 = 1*m/3; m2 = 2*m/3; m3 = m;
r0 = m0+1:m1; r1 = m1+1:m2; r2 = m2+1:m3; 
n0 = 0; n1 = 1*n/6; n2 = 2*n/6; n3 = 3*n/6; n4 = 4*n/6; n5 = 5*n/6; n6 = n;
c0 = n0+1:n1; c1 = n1+1:n2; c2 = n2+1:n3; c3 = n3+1:n4; c4 = n4+1:n5; c5 = n5+1:n6; 
t18 = A(r0,c4)-A(r2,c5);
t19 = A(r0,c1)+A(r1,c5);
t20 = A(r0,c2)-A(r1,c3);
t21 = A(r1,c2)-t18;
t22 = A(r2,c0)-A(r2,c3);
t23 = A(r1,c4)+A(r2,c1);
t24 = A(r2,c2)-t20;
t25 = t19-t24;
t26 = A(r0,c0)-A(r1,c1);
t27 = A(r0,c3)-A(r2,c0);
t28 = A(r1,c0)-A(r1,c3);
t29 = A(r2,c4)-A(r2,c5);
t30 = A(r1,c4)+A(r1,c5);
t31 = A(r2,c3)+t21;
t33 = A(r1,c0)+A(r0,c2)+A(r2,c2);
t34 = A(r0,c4)+A(r0,c5);
r58 = (A(r2,c1)-t29)/2;
r59 = (t22-t25)/2;
r60 = -(t19-t23)/2;
r61 = (t28+t30)/2;
r62 = (t27-A(r1,c2)-t18)/2;
r63 = (A(r2,c0)+A(r2,c3)-A(r1,c1)-t34)/2;
r64 = (A(r1,c4)-A(r2,c2)-t20)/2;
r65 = (t28-t30)/2;
r66 = (A(r1,c5)+t33)/2;
r67 = (t21-t27)/2;
r68 = (t26-t31)/2;
r69 = (t22+t25)/2;
r70 = (A(r2,c4)+A(r2,c5)-A(r0,c3)-t26)/2;
r71 = (A(r2,c1)+t29)/2;
r72 = (t31-A(r0,c0)-A(r1,c1))/2;
r73 = (t19+t23)/2;
oA2 = r70-r66;
oA3 = r66+r70;
oA5 = r63+r64;
oA6 = (A(r0,c5)+A(r1,c1)+A(r1,c4)-t22-t24-A(r0,c4))/2;
oA7 = r64-r63;
oA8 = r67+r73;
oA9 = r62+r60;
oA10 = r73-r67;
oA11 = r60-r62;
oA16 = r69-r58;
oA19 = r59+r71;
oA21 = (A(r0,c1)-A(r0,c3)+t23-t33-t34-A(r0,c0))/2;
oA25 = r65-r72;
oA27 = -r65-r72;
oA28 = r58+r69;
oA31 = r71-r59;
oA36 = r61-r68;
oA37 = -r61-r68;
v40 = oA6-oA9;
v41 = v40+oA19;
v42 = oA8-oA28;
v43 = v41-oA16;
v44 = oA10-oA31;
oA0 = v42+oA25;
oA1 = v44-oA27;
oA4 = v43-oA11;
v48 = oA2-oA7;
v49 = -oA3-oA5;
oA14 = oA1+v40;
v51 = oA4-oA37;
v53 = oA21-oA8-oA10;
oA12 = v43-oA0;
v55 = oA36-oA6;
oA34 = oA10-v51;
oA20 = v41+oA1-oA11;
oA23 = oA6+oA19-oA0;
oA30 = v40+v44-oA7;
oA24 = oA4-v49-v42;
oA29 = oA5-v42+v43;
oA26 = -v44-v48-oA6;
oA32 = oA37-oA12-oA2;
oA22 = oA3+v53-oA2;
oA33 = v55-oA8;
oA38 = oA0-v48-v51;
oA13 = oA8+v49;
oA18 = v48-v53;
oA17 = oA21-oA2-oA5;
oA39 = v49+v55-oA1;
oA35 = oA36-oA14-oA3;
oA15 = oA10+v48;

[m,n] = size(B);
m0 = 0; m1 = 1*m/6; m2 = 2*m/6; m3 = 3*m/6; m4 = 4*m/6; m5 = 5*m/6; m6 = m;
r0 = m0+1:m1; r1 = m1+1:m2; r2 = m2+1:m3; r3 = m3+1:m4; r4 = m4+1:m5; r5 = m5+1:m6; 
n0 = 0; n1 = 1*n/3; n2 = 2*n/3; n3 = n;
c0 = n0+1:n1; c1 = n1+1:n2; c2 = n2+1:n3; 
t18 = B(r3,c0)-B(r3,c2);
t19 = B(r0,c1)+B(r5,c2);
t20 = B(r0,c0)+B(r2,c2);
t21 = B(r1,c1)-B(r5,c0);
t22 = B(r1,c2)+B(r4,c0);
t23 = B(r3,c0)+B(r3,c2);
t24 = B(r0,c1)+B(r2,c1);
t25 = B(r4,c1)-t18;
t28 = B(r4,c2)+B(r2,c0)+B(r3,c1);
t29 = B(r0,c0)+B(r1,c2);
t30 = B(r1,c1)-B(r4,c1);
t31 = B(r4,c2)+B(r1,c0)-B(r5,c1);
t33 = B(r0,c2)+B(r2,c0)+B(r5,c1);
t34 = B(r2,c1)+t23;
r60 = (t21+t29)/4;
r61 = (t24+t30)/4;
r62 = (t19-t28)/4;
r63 = -(t20-t31)/4;
r64 = (t21-t29)/4;
r65 = (B(r4,c0)+B(r2,c0)-B(r1,c0)-B(r0,c0))/4;
r66 = (t23-t30)/4;
r67 = (B(r5,c0)+t19-t25)/4;
r68 = (t19+t28)/4;
r69 = (B(r0,c2)+B(r4,c2)-t21)/4;
r70 = (B(r3,c1)+t20+t22-B(r5,c1))/4;
r71 = (B(r2,c2)+B(r1,c2)-B(r5,c2))/4;
r72 = (B(r1,c1)+t25)/4;
r73 = (t20+t31)/4;
r74 = (t22-t33)/4;
r75 = (t18+t24)/4;
r76 = (B(r0,c1)-t34)/4;
r77 = (t22+t33)/4;
oB0 = r72-r77;
oB1 = r66+r74;
oB2 = r74-r66;
oB3 = r72+r77;
oB4 = r64-r62;
oB5 = r68-r60;
oB6 = r62+r64;
oB7 = r60+r68;
oB8 = r75-r73;
oB9 = r76+r63;
oB10 = r63-r76;
oB11 = r73+r75;
oB17 = t34/4+r69-r71;
oB19 = r69+r71+(B(r2,c1)+t18)/4;
oB24 = r65-r67;
oB25 = r65+r67;
oB32 = r61+r70;
oB33 = r61-r70;
v40 = oB10-oB1;
v41 = oB2-oB9;
oB13 = v40-oB4;
v43 = oB5-oB25;
v44 = oB7-oB3;
v45 = oB11-oB19;
v46 = oB5-oB17;
v47 = oB0-oB6;
v48 = oB8-oB33;
v49 = oB6-v45;
v51 = oB4-oB24;
v52 = oB32-oB11-v41;
oB20 = v51-oB11;
oB18 = oB8-v49-oB7;
oB14 = v44-oB11;
oB12 = oB5+v41;
oB39 = v46+oB13;
oB29 = oB5-v52;
oB28 = v48+oB13;
oB30 = oB32-oB2-oB7;
oB27 = v44-v51-oB1;
oB36 = -v44-v49;
oB35 = oB3-v52;
oB23 = v41+v43-oB0;
oB16 = oB10+v46-oB4-oB9;
oB37 = oB2-oB10+oB17;
oB38 = v45-oB0;
oB34 = oB0+v40+v48;
oB22 = oB24+oB13-oB3;
oB31 = oB1-oB6+oB33;
oB26 = oB2+v43-v47;
oB15 = oB8+v47;
oB21 = oB8+v43;

iC0 = FMMa( oA0, oB0, nmin, peeling, level);
iC1 = FMMa( oA1, oB1, nmin, peeling, level);
iC2 = FMMa( oA2, oB2, nmin, peeling, level);
iC3 = FMMa( oA3, oB3, nmin, peeling, level);
iC4 = FMMa( oA4, oB4, nmin, peeling, level);
iC5 = FMMa( oA5, oB5, nmin, peeling, level);
iC6 = FMMa( oA6, oB6, nmin, peeling, level);
iC7 = FMMa( oA7, oB7, nmin, peeling, level);
iC8 = FMMa( oA8, oB8, nmin, peeling, level);
iC9 = FMMa( oA9, oB9, nmin, peeling, level);
iC10 = FMMa( oA10, oB10, nmin, peeling, level);
iC11 = FMMa( oA11, oB11, nmin, peeling, level);
iC12 = FMMa( oA12, oB12, nmin, peeling, level);
iC13 = FMMa( oA13, oB13, nmin, peeling, level);
iC14 = FMMa( oA14, oB14, nmin, peeling, level);
iC15 = FMMa( oA15, oB15, nmin, peeling, level);
iC16 = FMMa( oA16, oB16, nmin, peeling, level);
iC17 = FMMa( oA17, oB17, nmin, peeling, level);
iC18 = FMMa( oA18, oB18, nmin, peeling, level);
iC19 = FMMa( oA19, oB19, nmin, peeling, level);
iC20 = FMMa( oA20, oB20, nmin, peeling, level);
iC21 = FMMa( oA21, oB21, nmin, peeling, level);
iC22 = FMMa( oA22, oB22, nmin, peeling, level);
iC23 = FMMa( oA23, oB23, nmin, peeling, level);
iC24 = FMMa( oA24, oB24, nmin, peeling, level);
iC25 = FMMa( oA25, oB25, nmin, peeling, level);
iC26 = FMMa( oA26, oB26, nmin, peeling, level);
iC27 = FMMa( oA27, oB27, nmin, peeling, level);
iC28 = FMMa( oA28, oB28, nmin, peeling, level);
iC29 = FMMa( oA29, oB29, nmin, peeling, level);
iC30 = FMMa( oA30, oB30, nmin, peeling, level);
iC31 = FMMa( oA31, oB31, nmin, peeling, level);
iC32 = FMMa( oA32, oB32, nmin, peeling, level);
iC33 = FMMa( oA33, oB33, nmin, peeling, level);
iC34 = FMMa( oA34, oB34, nmin, peeling, level);
iC35 = FMMa( oA35, oB35, nmin, peeling, level);
iC36 = FMMa( oA36, oB36, nmin, peeling, level);
iC37 = FMMa( oA37, oB37, nmin, peeling, level);
iC38 = FMMa( oA38, oB38, nmin, peeling, level);
iC39 = FMMa( oA39, oB39, nmin, peeling, level);

v52 = iC11-iC38;
b6 = iC19+iC0;
b8 = iC12-iC29;
v49 = iC37+iC4;
v54 = iC32+iC7;
b13 = iC27-iC1+iC14-iC20;
b15 = iC27-iC39-iC36+iC26+iC25;
b16 = iC23-iC35-iC21+iC20-iC33;
b20 = iC15-iC18+v54;
b23 = iC34-iC10-v52;
b25 = iC23+iC6-b6;
b29 = iC36-iC35-iC3+iC14+b25;
b31 = v49-iC12-iC16+b13;
b34 = iC2-iC17+b8-b23;
z31 = iC31+b23-iC15+b6-b13;
z30 = iC30+b31-iC2-v54;
z28 = iC28+b20-iC0-v49+b34;
oC3 = iC38-b15-iC37-iC24;
oC1 = iC22+iC32+iC34+b16;
z13 = iC13+iC3-b16+b15+b34;
z9 = iC9+iC26+v52+b25+b31;
z8 = iC8+b20-iC21+v52+b29;
z5 = iC5+iC25+b8-b29;
b45 = z31+z30;
b48 = z30-z31-z9;
t13 = z13+z8;
b51 = oC1+z13-z8-z5;
oC8 = z28+b45;
oC2 = t13-z9-z5;
b54 = oC3-z28+z5+t13;
oC0 = b51-z9;
oC6 = oC3+z9+b51;
oC5 = z28-b45-oC2;
oC4 = b54-b48;
oC7 = oC1+b48+b54;

C = [ oC0 oC1 oC2 ; oC3 oC4 oC5 ; oC6 oC7 oC8 ] ;
  end
end
end
