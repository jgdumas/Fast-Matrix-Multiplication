function C = FMM_6_3_3(A, B, nmin, peeling, level)
if nargin < 3, nmin = 3; end    % Threshold to conventional
if nargin < 4, peeling = 1; end % Static (1) or Dynamic (2) peeling
if nargin < 5, level = 3; end   % Verbose level
[m,k] = size(A); [k2,n] = size(B);
if (k2 ~= k), error('Incompatible matrix dimensions.'); end
% Recursively cuts into nmin*6^l x nmin*3^l x nmin*3^l blocks, with decreasing maximal l
if (m<=nmin)||(k<=nmin)||(n<=nmin)||(m<6)||(k<3)||(n<3)
  % fprintf("# MM Direct: %d x %d x %d\n",m,k,n)
  C = A*B;
else
  C = zeros(m,n);
  mu=m-rem(m,6);ku=k-rem(k,3);nu=n-rem(n,3);
  l=ceil(min([log(mu)/log(6),log(ku)/log(3),log(nu)/log(3)]));
  if (peeling == 1)
    mu=nmin;ku=nmin;nu=nmin;l=0;
    while (mu <= m) && (ku <= k) && (nu <= n)
      l=l+1; mu=mu*6; ku=ku*3; nu=nu*3;
    end
    l=l-1;mu=nmin*6^l; ku=nmin*3^l; nu=nmin*3^l;
  end
  if (mu < m) || (ku < k) || (nu < n)
    % fprintf("# Core SubMatrix[%d]: %d x %d x %d\n",l,mu,ku,nu)
    C(1:mu,1:nu)=FMM_6_3_3(A(1:mu,1:ku),B(1:ku,1:nu),nmin, peeling, level);
    if (m > mu)
      % fprintf("# MM peel m: %d x %d x %d\n",m-mu,k,n)
      C(mu+1:m,1:n)=C(mu+1:m,1:n)+FMM(A(mu+1:m,1:k),B,nmin, peeling, level);
    end
    if (k > ku) && (mu > 0) && (nu > 0)
      % fprintf("# MM peel k: %d x %d x %d\n",mu,k-ku,nu)
      C(1:mu,1:nu)=C(1:mu,1:nu)+FMM(A(1:mu,ku+1:k),B(ku+1:k,1:nu),nmin, peeling, level);
    end
    if (n > nu) && (mu > 0)
      % fprintf("# MM peel n: %d x %d x %d\n",mu,k,n-nu)
      C(1:mu,nu+1:n)=C(1:mu,nu+1:n)+FMM(A(1:mu,1:k),B(1:k,nu+1:n),nmin, peeling, level);
    end
  else
    if l>=level, fprintf("# Core<6;3;3>[%d]: %d x %d x %d\n",l,m,k,n); end

[m,n] = size(A);
m0 = 0; m1 = 1*m/6; m2 = 2*m/6; m3 = 3*m/6; m4 = 4*m/6; m5 = 5*m/6; m6 = m;
r0 = m0+1:m1; r1 = m1+1:m2; r2 = m2+1:m3; r3 = m3+1:m4; r4 = m4+1:m5; r5 = m5+1:m6; 
n0 = 0; n1 = 1*n/3; n2 = 2*n/3; n3 = n;
c0 = n0+1:n1; c1 = n1+1:n2; c2 = n2+1:n3; 
t18 = A(r2,c0)-A(r2,c2);
t19 = A(r2,c0)+A(r2,c2);
t20 = A(r3,c1)-A(r1,c0);
t21 = A(r5,c1)+t18;
r22 = A(r4,c2)/8;
t22 = A(r5,c2)-r22;
r23 = A(r0,c2)/8;
t23 = A(r3,c2)-r23;
r24 = A(r0,c1)/8;
t24 = A(r1,c2)-r24;
r25 = A(r0,c0)/8;
t25 = A(r3,c0)+r25;
r26 = A(r4,c0)/8;
t26 = A(r5,c0)+r26;
r27 = A(r4,c1)/8;
t27 = t19-r27;
t28 = t20-t21;
t29 = t26-t27;
t31 = A(r5,c0)+A(r2,c1)-r22;
t32 = A(r3,c2)-r25;
t34 = t20+t32;
t35 = t18+r27;
t36 = A(r1,c1)-A(r3,c1);
t38 = -t22-t28;
t42 = -t24-t29;
t43 = A(r1,c1)+A(r5,c2)+r22-t25;
t44 = t22-t28;
t45 = t24+t31;
t46 = A(r1,c1)+r25-t22-A(r3,c0);
t47 = A(r1,c0)+t25;
t48 = A(r1,c1)+A(r3,c1);
t49 = t23+A(r1,c2);
t51 = t21-r24;
t52 = t23+A(r5,c0)-r26;
t53 = A(r5,c1)+t19+r24;
t54 = A(r1,c2)-t23;
t55 = A(r3,c2)+r23;
t56 = t26+t27;
oA0 = t46-t51;
oA1 = -t46-t51;
oA2 = t43-t53;
oA3 = t43+t53;
oA4 = t48+t52-t35;
oA5 = t55-t56-t36;
oA6 = t35+t36+t52;
oA7 = t29+t48-t55;
oA12 = t34+t45;
oA13 = t45-t34;
oA15 = t24-t31+t32-t20;
oA16 = t38-t49;
oA19 = t54-t44;
oA22 = t42-t47;
oA28 = t44+t54;
oA31 = t49+t38;
oA32 = t47+t42;
oA35 = t24+t25-t56-A(r1,c0);
v40 = oA0+oA12;
v42 = v40+oA16-oA19;
v43 = oA7+oA15;
v44 = oA13+oA5;
oA14 = oA1+v42;
oA10 = v43-oA2;
v47 = oA4-oA32;
v48 = oA35+oA14;
oA8 = oA3+v44;
oA9 = oA6-v42;
v51 = oA13-oA28;
v52 = oA31-oA15;
oA21 = oA22+v43+v44;
oA23 = oA6+oA19-oA0;
oA34 = oA12+v43-v47;
oA37 = oA2+oA12+oA32;
oA11 = v40-oA4;
oA24 = oA4-v51;
oA29 = v40-v51-oA3;
oA38 = oA7+v40-v47;
oA30 = v42-v52-oA2;
oA17 = oA13+oA22+oA10;
oA25 = oA0+oA28-oA8;
oA36 = oA3+v48;
oA20 = oA1+oA4+oA16;
oA27 = oA10-oA1-oA31;
oA39 = oA35-oA9-oA5;
oA33 = v48-oA6-v44;
oA26 = v52-oA6;
oA18 = oA3-oA7-oA22;

[m,n] = size(B);
m0 = 0; m1 = 1*m/3; m2 = 2*m/3; m3 = m;
r0 = m0+1:m1; r1 = m1+1:m2; r2 = m2+1:m3; 
n0 = 0; n1 = 1*n/3; n2 = 2*n/3; n3 = n;
c0 = n0+1:n1; c1 = n1+1:n2; c2 = n2+1:n3; 
t9 = -B(r2,c1)-B(r1,c1);
t10 = B(r2,c0)-B(r0,c0);
t11 = B(r2,c0)+B(r0,c0);
t12 = B(r0,c2)-t9;
t13 = B(r1,c0)+B(r2,c0);
t14 = B(r0,c1)+B(r2,c1);
t15 = B(r1,c2)+t12;
oB8 = t12-B(r1,c2)-t11;
oB11 = t10+t15;
oB12 = t15-t10;
oB19 = B(r2,c2)-t10-B(r0,c2);
oB21 = t11-t14;
oB22 = t11+t14;
oB23 = B(r0,c1)-B(r2,c1)+t10;
oB24 = t9-t13;
oB25 = -t9-t13;
v40 = oB25-oB24;
oB0 = v40-oB11;
oB33 = oB8+oB21;
oB32 = oB23+oB12;
oB28 = oB19-oB0;
oB7 = oB32-oB22;
oB13 = oB33+oB22;
oB34 = oB23+oB0;
oB3 = v40-oB8;
oB18 = -oB8-oB28;
oB38 = -oB11-oB24;
oB35 = oB3-oB22;
oB29 = -oB11-oB19;
oB6 = oB23+oB33;
oB20 = oB32-oB11;
oB17 = -oB28-oB13;
oB2 = oB21+oB32;
oB26 = oB21+oB23+oB25;
oB14 = oB8-oB0-oB7;
oB1 = oB33+oB20;
oB5 = v40-oB13;
oB9 = oB23+oB35;
oB4 = v40-oB12;
oB15 = -oB21-oB34;
oB10 = oB22-oB34;
oB30 = -oB7-oB18;
oB16 = oB12+oB29;
oB27 = -oB7-oB38;
oB36 = oB8-oB25;
oB39 = -oB24-oB13;
oB31 = oB19+oB6;
oB37 = oB25-oB12;

iC0 = FMM( oA0, oB0, nmin, peeling, level);
iC1 = FMM( oA1, oB1, nmin, peeling, level);
iC2 = FMM( oA2, oB2, nmin, peeling, level);
iC3 = FMM( oA3, oB3, nmin, peeling, level);
iC4 = FMM( oA4, oB4, nmin, peeling, level);
iC5 = FMM( oA5, oB5, nmin, peeling, level);
iC6 = FMM( oA6, oB6, nmin, peeling, level);
iC7 = FMM( oA7, oB7, nmin, peeling, level);
iC8 = FMM( oA8, oB8, nmin, peeling, level);
iC9 = FMM( oA9, oB9, nmin, peeling, level);
iC10 = FMM( oA10, oB10, nmin, peeling, level);
iC11 = FMM( oA11, oB11, nmin, peeling, level);
iC12 = FMM( oA12, oB12, nmin, peeling, level);
iC13 = FMM( oA13, oB13, nmin, peeling, level);
iC14 = FMM( oA14, oB14, nmin, peeling, level);
iC15 = FMM( oA15, oB15, nmin, peeling, level);
iC16 = FMM( oA16, oB16, nmin, peeling, level);
iC17 = FMM( oA17, oB17, nmin, peeling, level);
iC18 = FMM( oA18, oB18, nmin, peeling, level);
iC19 = FMM( oA19, oB19, nmin, peeling, level);
iC20 = FMM( oA20, oB20, nmin, peeling, level);
iC21 = FMM( oA21, oB21, nmin, peeling, level);
iC22 = FMM( oA22, oB22, nmin, peeling, level);
iC23 = FMM( oA23, oB23, nmin, peeling, level);
iC24 = FMM( oA24, oB24, nmin, peeling, level);
iC25 = FMM( oA25, oB25, nmin, peeling, level);
iC26 = FMM( oA26, oB26, nmin, peeling, level);
iC27 = FMM( oA27, oB27, nmin, peeling, level);
iC28 = FMM( oA28, oB28, nmin, peeling, level);
iC29 = FMM( oA29, oB29, nmin, peeling, level);
iC30 = FMM( oA30, oB30, nmin, peeling, level);
iC31 = FMM( oA31, oB31, nmin, peeling, level);
iC32 = FMM( oA32, oB32, nmin, peeling, level);
iC33 = FMM( oA33, oB33, nmin, peeling, level);
iC34 = FMM( oA34, oB34, nmin, peeling, level);
iC35 = FMM( oA35, oB35, nmin, peeling, level);
iC36 = FMM( oA36, oB36, nmin, peeling, level);
iC37 = FMM( oA37, oB37, nmin, peeling, level);
iC38 = FMM( oA38, oB38, nmin, peeling, level);
iC39 = FMM( oA39, oB39, nmin, peeling, level);

b1 = iC24+iC3;
b4 = iC31+iC16+iC1;
b8 = iC2-iC30;
v49 = iC37+iC17;
b13 = iC25-iC0+iC19;
v48 = iC36+iC18;
b22 = iC35-iC21+b8;
b24 = iC35-iC20+b1;
b26 = iC34+iC21-b13;
b28 = iC34+iC20-b4;
z39 = iC39+iC16-v49;
z38 = iC38+v48-iC19;
z33 = iC33+iC31-iC28-iC34;
z32 = iC32+iC30-iC35-iC29;
z23 = iC23+iC26-iC25+iC21;
z22 = iC22+iC24-iC20+iC27;
z15 = iC15+v48+b26;
z14 = iC14-iC36-b24;
z13 = iC13+iC28-iC24+v49+b28;
z12 = iC12+iC25+iC29+iC37+b22;
z11 = iC11+iC29-iC36+iC19-b1;
z10 = iC10+b4-iC37-iC27;
z9 = iC9+iC26-iC16+iC37+b8;
z8 = iC8+iC28-iC36+b13;
z7 = iC7+iC27-iC30-iC18+b24;
z6 = iC6+iC26-iC31+b26;
z5 = iC5+iC17-b22;
z4 = iC4+b28;
t52 = z39+z38;
b35 = z38-z39;
b36 = z33+z23;
t37 = z33-z23;
t56 = z32-z22;
b37 = z32+z22;
b38 = -z15-z14;
t54 = z15-z14;
b39 = -z13-z12;
b40 = z12-z13;
t50 = z11+z9;
t36 = z11-z9;
t55 = z10-z8;
b41 = -z10-z8;
t42 = z7-z6;
t40 = z7+z6;
t58 = z5-z4;
b42 = -z5-z4;
b43 = t54-z38;
t23 = (z39+t58-b39)/8;
t45 = -b40;
t27 = t52+t50;
b46 = b39-b36-b37;
b47 = t42-b42;
t30 = b42+t42;
t32 = t50+t55;
t29 = t56-t37;
b48 = t56+t40+t37;
t33 = (t36-b41)/8;
oC2 = t45+b38;
t24 = (b43-t40)/8;
oC9 = b35/8+t33;
oC16 = (-b37+b36-t32)/8;
oC14 = t32+b47;
t28 = b47/8;
t20 = b40+t29;
oC6 = t52/8+t28;
b53 = t55-t27;
oC3 = (t55+t27+oC2)/8;
b55 = b41+t36+t29;
oC5 = (t30+t29)/8;
oC17 = t23-t24;
oC15 = t24+t23;
oC4 = (b48-t58)/8;
oC1 = t58+b48-b53;
oC11 = b53/8;
oC12 = b35+b38+t20;
oC10 = t28+(t20-b38)/8;
oC13 = b43-z39+b46;
oC7 = t33+(t54-b46)/8;
oC0 = b35+t30+b55;
oC8 = (b38-t45+b55)/8;

C = [ oC0 oC1 oC2 ; oC3 oC4 oC5 ; oC6 oC7 oC8 ; oC9 oC10 oC11 ; oC12 oC13 oC14 ; oC15 oC16 oC17 ] ;
  end
end
end
