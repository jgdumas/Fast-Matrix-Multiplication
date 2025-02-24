function C = FMM_6_3_3(A, B, nmin, peeling, level)
if nargin < 3, nmin = 4; end    % Threshold to conventional
if nargin < 4, peeling = 1; end % Static (1) or Dynamic (2) peeling
if nargin < 5, level = 8; end   % Verbose level
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
    l = min([floor(log(m/nmin)/log(6)),floor(log(k/nmin)/log(3)),floor(log(n/nmin)/log(3))]);
    mu=nmin*6^l; ku=nmin*3^l; nu=nmin*3^l;
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
t18 = A(r2,c0)+A(r5,c1);
t19 = A(r3,c1)-A(r1,c0);
t20 = A(r2,c2)-A(r5,c2);
t21 = A(r3,c2)-A(r1,c2);
r22 = A(r4,c2)/8;
t22 = t19+r22;
t25 = t18-t22;
r27 = A(r4,c0)/8;
t27 = A(r3,c2)-r27;
r28 = (A(r0,c0)-A(r0,c1))/8;
t28 = A(r5,c2)-r28;
r29 = (A(r0,c0)+A(r0,c1))/8;
t29 = A(r1,c2)-r29;
t30 = A(r3,c0)-A(r1,c1)-r22;
r31 = A(r4,c1)/8;
t31 = A(r2,c2)-r31;
t32 = A(r5,c0)+A(r2,c1);
t33 = A(r1,c1)+t27;
t34 = A(r2,c0)-A(r3,c1);
t36 = t32-t22;
t38 = A(r3,c2)+t29;
r39 = A(r0,c2)/8;
t39 = r39-t21;
t40 = A(r1,c0)+t29-t31-A(r5,c0);
t43 = t21+r39;
t45 = A(r2,c2)+t28;
t46 = -t18-t30;
t48 = A(r5,c1)+t28+t33;
t49 = t18-t30;
t50 = t20+t25;
t51 = t21-r28;
t52 = t20+r29;
t54 = A(r2,c1)-A(r3,c1)+r31;
t55 = t31+r39-A(r5,c0);
t56 = t20-t25;
t57 = t19+t32-r22;
t58 = A(r3,c0)+A(r2,c0)+r27;
oA0 = t52+t46;
oA1 = t45-t49;
oA2 = t46-t52;
oA3 = t45+t49;
oA6 = t33+t34-t55;
oA7 = A(r1,c1)-t27-t34-t55;
oA12 = t38+t57;
oA13 = t36-t51;
oA14 = t36+t51;
oA15 = t38-t57;
oA17 = t50-t43;
oA19 = t56+t39;
oA21 = t40-t58;
oA24 = t48-t54;
oA28 = t39-t56;
oA30 = t43+t50;
oA34 = t40+t58;
oA36 = t48+t54;
v40 = oA13-oA28;
v41 = oA2+oA17;
v42 = oA3-v41;
oA4 = v40+oA24;
v44 = oA0+oA12;
v45 = oA21+v42;
v46 = oA15+oA7;
v47 = oA14-oA1;
v49 = v47-oA2-oA30;
v51 = oA36-v45-oA6;
v53 = oA34+oA4-v46;
v54 = v44-oA19;
oA11 = v44-oA4;
oA9 = oA6-v47;
oA5 = oA21-v41;
oA20 = oA14+oA4-v54;
oA10 = v46-oA2;
oA31 = oA15+v49;
oA32 = v53-oA12;
oA18 = oA13+oA15+v42;
oA27 = oA7-oA14+oA30;
oA33 = v51-oA13;
oA22 = v41-v46-oA13;
oA26 = v49-oA6;
oA39 = v51-oA1;
oA35 = oA36-oA3-oA14;
oA8 = oA13+v45;
oA25 = oA0-v40-v45;
oA16 = v47-v54;
oA38 = oA0-oA15+oA34;
oA29 = v44-oA3-v40;
oA23 = oA6+oA19-oA0;
oA37 = oA2+v53;

[m,n] = size(B);
m0 = 0; m1 = 1*m/3; m2 = 2*m/3; m3 = m;
 r0 = m0+1:m1; r1 = m1+1:m2; r2 = m2+1:m3;
n0 = 0; n1 = 1*n/3; n2 = 2*n/3; n3 = n;
 c0 = n0+1:n1; c1 = n1+1:n2; c2 = n2+1:n3;
t9 = B(r2,c1)+B(r2,c0);
t10 = -B(r1,c1)-B(r1,c2);
t11 = B(r0,c0)+B(r0,c2);
t12 = B(r2,c1)+B(r2,c2);
t14 = -t10-t11;
t15 = t11+B(r1,c1)-B(r1,c2);
t16 = B(r2,c0)-B(r2,c1);
oB1 = t15-t9;
oB5 = t14-t16;
oB9 = t14+t16;
oB13 = t9+t15;
oB18 = B(r0,c0)-B(r0,c2)+B(r2,c0)-B(r2,c2);
oB22 = B(r0,c0)+B(r0,c1)+t9;
oB27 = t9-B(r1,c0)-B(r1,c1);
oB29 = t10-t12;
oB30 = t12+t10;
oB3 = oB18-oB29;
v41 = oB30+oB9;
oB17 = oB5+oB29;
v43 = oB13-oB1;
oB7 = -oB30-oB18;
oB31 = oB1+v41;
oB14 = v43-oB3;
oB39 = -oB27-oB1;
oB10 = -oB17-oB31;
oB35 = oB3-oB22;
oB28 = -oB13-oB17;
v51 = oB13-oB3;
oB4 = v41-oB28;
oB8 = oB5+v51;
oB2 = -oB30-oB17;
oB12 = -v41-oB29;
oB26 = oB9-oB39;
oB16 = -v41;
oB32 = oB22+oB7;
oB20 = oB22-v43;
oB23 = oB9-oB35;
oB24 = oB27-v43;
oB21 = oB35-oB5;
oB11 = v43+oB7;
oB36 = oB14-oB27;
oB19 = oB30-oB14;
oB15 = -oB18-oB31;
oB25 = oB5-oB39;
oB33 = oB13-oB22;
oB38 = -oB7-oB27;
oB37 = oB27-oB10;
oB34 = oB22-oB10;
oB0 = oB1+oB5-oB7;
oB6 = oB9+v51;

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
b1 = iC37-iC16;
v50 = iC38+iC18;
v48 = iC16+iC39;
b16 = iC16+iC31+iC1;
b25 = iC2-iC30+b1;
b31 = iC34+iC4-b16;
b33 = iC35-iC39-iC5+b25;
v44 = iC24+b31;
b34 = iC3+iC35+v44;
b35 = b33-iC25;
b38 = iC34-iC38+iC0+b35;
z33 = iC33+iC31-iC34-iC28;
z32 = iC32+iC30-iC29-iC35;
z23 = iC23+iC26+b35;
z22 = iC22+iC27+v44;
z21 = iC21-b33;
z20 = iC20+b31;
z19 = iC19-iC36-v50;
z17 = iC17+b1-iC39;
z15 = iC15+b38;
z14 = iC14-iC36-b34;
z13 = iC13+iC28-iC24-iC4+v48;
z12 = iC12+iC29+iC5+iC25+v48;
z11 = iC11+iC29-iC24-iC3+v50;
z10 = iC10+b16-iC37-iC27;
z9 = iC9+iC26+b25;
z8 = iC8+iC28-iC0+iC25+v50;
z7 = iC7+iC27-iC18-iC30+b34;
z6 = iC6+iC26-iC18-iC36-iC31+b38;
t39 = z33+z32;
b40 = z32-z33;
b41 = z23+z22;
t40 = z23-z22;
t57 = z21-z20;
b42 = z21+z20;
t47 = z19+z17;
b44 = -z15-z14;
t37 = z15-z14;
b45 = z12-z13;
b46 = -z13-z12;
t44 = z11+z10;
t35 = z11-z10;
t54 = z9-z8;
t45 = z9+z8;
t52 = z7-z6;
t51 = z7+z6;
t56 = (z17-z19)/8;
t46 = -t56;
b47 = b44-b45;
t33 = (z17-b46)/8;
t27 = (t54+t44)/8;
t26 = t44+t47;
b50 = b46-t39-b41;
b51 = t51+t39-b41;
oC6 = t52/8;
t29 = t57+t37;
t24 = (z19-t51+t37)/8;
b54 = t35-t45;
t30 = t45+t35;
t28 = t40+b40;
oC2 = t47+b47;
t25 = b47/8;
oC11 = t56-t30/8;
t20 = t52+t28;
oC12 = b45+b44-b42+t28;
oC16 = (t40+b42-b40)/8-t27;
oC14 = t54+t52+t26;
oC9 = (t26-t54)/8;
oC3 = t27+t25;
oC17 = t33-t24;
oC15 = t33+t24;
oC13 = t29+b50;
oC1 = t57+t30+b51;
oC4 = (b51-t57-t47)/8;
oC7 = (t29-b50)/8+oC9;
t19 = t20-b42;
oC10 = (t20+b42)/8-t25;
oC0 = b54+t19;
oC5 = t46+t19/8;
oC8 = t46+(b54+oC12)/8;

C = [ oC0 oC1 oC2 ; oC3 oC4 oC5 ; oC6 oC7 oC8 ; oC9 oC10 oC11 ; oC12 oC13 oC14 ; oC15 oC16 oC17 ] ;
  end
end
end
