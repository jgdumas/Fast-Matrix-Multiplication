function C = FMMa_3_6_3(A, B, nmin, peeling, level)
if nargin < 3, nmin = 4; end    % Threshold to conventional
if nargin < 4, peeling = 1; end % Static (1) or Dynamic (2) peeling
if nargin < 5, level = 8; end   % Verbose level
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
    mu=nmin*3^l; ku=nmin*6^l; nu=nmin*3^l;
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
t18 = A(r1,c5)-A(r2,c5);
t19 = A(r0,c0)-A(r2,c3);
t20 = A(r0,c1)-A(r1,c3);
t21 = A(r1,c0)-A(r2,c1);
t22 = A(r1,c0)-A(r1,c1);
t23 = A(r0,c4)+t18;
t24 = A(r1,c2)-A(r1,c4);
t26 = A(r0,c3)+A(r2,c2);
t27 = A(r0,c5)+A(r1,c2)-A(r2,c4);
t28 = A(r0,c0)+A(r0,c2);
t29 = t19+t23;
t30 = A(r0,c1)-A(r0,c4);
t32 = A(r0,c3)-A(r2,c0);
t33 = A(r2,c1)-A(r2,c4);
t34 = A(r2,c0)-A(r2,c2);
t35 = A(r1,c5)+A(r2,c5);
t36 = t20+A(r0,c2)-A(r2,c3);
r56 = t18/2;
r57 = t30/2;
r58 = (t26-t28)/2;
r59 = (-t21+t27)/2;
r60 = (t22+t24)/2;
r61 = (A(r2,c4)-t22)/2;
r62 = (A(r1,c4)+A(r1,c2)+A(r2,c1)-t32)/2;
r63 = (-t22+t24)/2;
r64 = (t21+A(r0,c5)-A(r1,c4))/2;
r65 = (t18+t36)/2;
r66 = (t30+t35)/2;
r67 = (A(r1,c3)+t29)/2;
r68 = (-t19+t20)/2;
r69 = (A(r1,c3)-t29)/2;
r70 = (t21+t27)/2;
r71 = (t19+t20)/2;
r72 = (t33-t34)/2;
r73 = (t26+t28)/2;
oA1 = r58-r61-r56;
oA2 = r61-r73-r56;
oA8 = r62+r66;
oA10 = (A(r0,c1)+A(r2,c1)+t23-t24+t32)/2;
oA11 = r62-r66;
oA12 = r59-r71;
oA13 = r71+r59;
oA14 = r70+r68;
oA15 = r68-r70;
oA16 = r65-r72;
oA17 = r65+r72;
oA20 = r57+r58+r64;
oA21 = r57-r64-r73;
oA26 = r60+r67;
oA27 = r67-r60;
oA29 = (t33+t34-t35-t36)/2;
oA37 = r69+r63;
oA39 = r63-r69;
v40 = oA10-oA15;
oA5 = oA21-oA17-oA2;
v43 = v40-oA27;
v44 = v43-oA26;
oA4 = oA20-oA16-oA1;
v47 = oA29-oA5;
oA9 = v44-oA14;
v49 = oA8-oA13;
v50 = oA11+oA4;
v51 = oA37-oA12;
oA32 = v51-oA2;
oA23 = oA12+oA16+oA9;
oA28 = oA8+v47-v50;
oA3 = v49-oA5;
oA35 = oA39+oA5+oA9;
oA34 = oA10+oA37-oA4;
oA31 = oA10-oA27-oA1;
oA7 = oA2+v40;
oA0 = v50-oA12;
oA6 = v44-oA1;
oA18 = oA8+oA15-oA21;
oA19 = oA11-oA14+oA20;
oA36 = oA39+v44+v49;
oA30 = oA14-v43-oA2;
oA38 = oA11+v40+v51;
oA25 = v47-oA12;
oA33 = oA1-oA13+oA39;
oA24 = v47+v49-oA11;
oA22 = oA17-oA10-oA13;

[m,n] = size(B);
m0 = 0; m1 = 1*m/6; m2 = 2*m/6; m3 = 3*m/6; m4 = 4*m/6; m5 = 5*m/6; m6 = m;
 r0 = m0+1:m1; r1 = m1+1:m2; r2 = m2+1:m3; r3 = m3+1:m4; r4 = m4+1:m5; r5 = m5+1:m6;
n0 = 0; n1 = 1*n/3; n2 = 2*n/3; n3 = n;
 c0 = n0+1:n1; c1 = n1+1:n2; c2 = n2+1:n3;
t18 = B(r2,c2)+B(r4,c0);
t19 = B(r0,c0)+B(r1,c2);
t20 = B(r2,c1)-B(r4,c1);
t21 = B(r3,c1)+t18;
t22 = B(r0,c1)+B(r1,c1);
t23 = B(r5,c0)+B(r5,c2);
t24 = B(r0,c1)-B(r3,c2);
t25 = B(r2,c1)+B(r3,c0);
t26 = B(r5,c0)-B(r5,c2);
t27 = B(r0,c2)+B(r1,c0);
t28 = t19+t21;
t31 = B(r0,c0)+B(r2,c2);
t32 = B(r1,c0)+B(r4,c2)-B(r5,c1);
t33 = B(r2,c0)+B(r3,c1)+B(r4,c2);
t34 = B(r1,c1)-t23;
t35 = B(r4,c1)+t26;
r53 = (t25+t31)/4;
r54 = (t20-t22)/4;
r55 = (B(r5,c1)+t28)/4;
r56 = (t25-t31)/4;
r57 = (t24-t32)/4;
r58 = (-B(r0,c1)+t34)/4;
r59 = (t24+t32)/4;
r60 = (t22-t26)/4;
r61 = (B(r2,c1)+t35)/4;
r62 = (t20+t22)/4;
r63 = (t21-t27)/4;
r64 = (t19+t33)/4;
r65 = (t19-t33)/4;
r66 = (B(r3,c1)-t18-t27)/4;
r67 = (B(r5,c1)-t28)/4;
r68 = (t20+t23)/4;
oB4 = r58-r65;
oB5 = -r58-r65;
oB6 = r60-r64;
oB7 = r60+r64;
oB8 = r56+r57;
oB9 = r59-r53;
oB10 = r56-r57;
oB11 = r53+r59;
oB12 = r63+r68;
oB13 = r68-r63;
oB14 = r66-r61;
oB15 = r61+r66;
oB16 = (B(r1,c2)+B(r2,c2)+B(r3,c2)-B(r4,c2)+t25-t34-B(r0,c2))/4;
oB22 = r54-r55;
oB23 = r54+r55;
oB32 = r62-r67;
oB33 = r62+r67;
oB38 = (B(r0,c0)+B(r1,c0)+B(r2,c0)+B(r3,c0)+B(r4,c0)+t24-t35)/4;
v40 = oB4+oB13;
v41 = oB6+oB15;
v42 = oB12-oB5;
v43 = oB7-oB14;
oB1 = oB10-v40;
oB0 = v41-oB8;
v46 = oB38+oB15;
v47 = oB32-v42;
v48 = oB16+oB4;
v49 = v43+oB22;
oB18 = v46-oB7;
oB30 = v47-oB7-oB9;
oB2 = oB9+v42;
oB28 = oB8+oB13-oB33;
oB20 = v40-v49;
oB27 = oB7-oB10+oB22;
oB35 = v43-v47;
oB21 = oB23+v41-v42;
oB19 = oB11-oB38-oB0;
oB25 = oB12-oB23-oB0;
oB24 = v49-oB11-oB13;
oB17 = oB5-oB9+oB10-v48;
oB37 = oB12-v48;
oB31 = oB33+oB1-oB6;
oB39 = oB9+oB16-oB1;
oB34 = v40+v41-oB33;
oB29 = oB11+oB12-oB32;
oB26 = oB6+oB9+oB23;
oB36 = v46-oB8-oB11-oB14;
oB3 = v43-oB11;

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
b1 = iC31-iC16;
b3 = iC34-iC15;
v48 = iC37+iC26;
v46 = iC27-iC39;
b9 = iC3-iC36;
v51 = iC0+iC19;
b13 = iC31-iC10-iC27+iC1;
v54 = iC22-iC7;
v55 = iC30-iC17;
v47 = iC35-iC14;
b15 = iC6+iC23;
b16 = iC9-iC10;
b21 = iC8-iC18-b9;
b23 = iC38-iC11+v51;
b25 = b15-iC12;
b34 = iC25-iC29+b9+v51-b25;
b37 = b23-b15-iC32-iC7;
z33 = iC33+iC6-v54+b13+b21;
oC8 = iC28+iC19-iC18-iC29+b1+v55;
oC3 = iC36+iC38-v48-v46-iC25-iC24;
z21 = iC21+v47-iC15-b21+b37;
z20 = iC20+b3-iC14+b23+b13;
z13 = iC13+v55+v54-b16+v46+b34;
z5 = iC5+v47+b34;
z4 = iC4+b1+b3+b16+v48+b25;
z2 = iC2+iC26-iC30+iC9-b37;
oC1 = z20-z33-z21;
b46 = z13-oC8;
b48 = z4-z20;
t15 = z5+z2;
b50 = z2-z5;
t11 = z21+b50;
oC2 = z33+z13-z4+b50;
b54 = z33+z4+b46;
t9 = b48-t11;
oC6 = oC3+z13+b48+t11;
oC5 = t15-b54;
oC4 = oC3+t15+b54;
oC0 = z13-t9;
oC7 = oC3+b46+t9;

C = [ oC0 oC1 oC2 ; oC3 oC4 oC5 ; oC6 oC7 oC8 ] ;
  end
end
end
