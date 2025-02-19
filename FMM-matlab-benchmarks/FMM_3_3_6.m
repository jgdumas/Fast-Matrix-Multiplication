function C = FMM_3_3_6(A, B, nmin, peeling, level)
if nargin < 3, nmin = 3; end    % Threshold to conventional
if nargin < 4, peeling = 1; end % Static (1) or Dynamic (2) peeling
if nargin < 5, level = 3; end   % Verbose level
[m,k] = size(A); [k2,n] = size(B);
if (k2 ~= k), error('Incompatible matrix dimensions.'); end
% Recursively cuts into nmin*3^l x nmin*3^l x nmin*6^l blocks, with decreasing maximal l
if (m<=nmin)||(k<=nmin)||(n<=nmin)||(m<3)||(k<3)||(n<6)
  % fprintf("# MM Direct: %d x %d x %d\n",m,k,n)
  C = A*B;
else
  C = zeros(m,n);
  mu=m-rem(m,3);ku=k-rem(k,3);nu=n-rem(n,6);
  l=ceil(min([log(mu)/log(3),log(ku)/log(3),log(nu)/log(6)]));
  if (peeling == 1)
    l = min([floor(log(m/nmin)/log(3)),floor(log(k/nmin)/log(3)),floor(log(n/nmin)/log(6))]);
    mu=nmin*2^l; ku=nmin*2^l; nu=nmin*2^l;
  end
  if (mu < m) || (ku < k) || (nu < n)
    % fprintf("# Core SubMatrix[%d]: %d x %d x %d\n",l,mu,ku,nu)
    C(1:mu,1:nu)=FMM_3_3_6(A(1:mu,1:ku),B(1:ku,1:nu),nmin, peeling, level);
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
    if l>=level, fprintf("# Core<3;3;6>[%d]: %d x %d x %d\n",l,m,k,n); end

[m,n] = size(A);
m0 = 0; m1 = 1*m/3; m2 = 2*m/3; m3 = m;
r0 = m0+1:m1; r1 = m1+1:m2; r2 = m2+1:m3; 
n0 = 0; n1 = 1*n/3; n2 = 2*n/3; n3 = n;
c0 = n0+1:n1; c1 = n1+1:n2; c2 = n2+1:n3; 
t9 = A(r0,c0)-A(r2,c0);
t10 = A(r2,c1)+t9;
t11 = A(r0,c2)-A(r1,c2);
t12 = A(r0,c0)+A(r2,c0);
t13 = A(r0,c1)+A(r2,c1);
t14 = A(r0,c2)-A(r2,c2);
t15 = t10-A(r1,c1);
oA6 = t11-t15;
oA9 = -t11-t15;
oA16 = t9+t14;
oA18 = t12-A(r0,c2)-A(r2,c2);
oA19 = t9-t14;
oA21 = t12-t13;
oA22 = t12+t13;
oA23 = A(r0,c1)-t10;
oA39 = A(r1,c0)-t11-A(r0,c0);
v40 = oA18-oA22;
oA31 = oA6+oA19;
oA35 = oA9-oA23;
oA33 = oA6-oA23;
v44 = oA16+oA9;
oA17 = v40-oA21;
oA15 = -oA31-oA18;
oA29 = v40-oA35;
oA8 = oA33-oA21;
oA26 = oA9-oA39;
oA34 = -oA15-oA21;
oA13 = oA22+oA33;
oA5 = oA35-oA21;
v53 = -oA16-oA39;
oA28 = -oA18-oA8;
oA3 = oA22+oA35;
oA24 = -oA39-oA13;
oA1 = oA16+oA31;
oA20 = oA16+oA19+oA23;
oA2 = v44-oA17;
oA7 = v44-oA18;
oA4 = -oA28-oA16;
oA14 = -v44-oA19;
oA25 = oA5-oA39;
oA10 = -oA17-oA31;
oA36 = oA6-oA26;
oA27 = v53-oA31;
oA30 = -v44;
oA0 = oA34-oA23;
oA12 = oA16-oA29;
oA37 = v53+oA17;
oA11 = -oA29-oA19;
oA32 = v44-v40;
oA38 = -oA15-oA26;

[m,n] = size(B);
m0 = 0; m1 = 1*m/3; m2 = 2*m/3; m3 = m;
r0 = m0+1:m1; r1 = m1+1:m2; r2 = m2+1:m3; 
n0 = 0; n1 = 1*n/6; n2 = 2*n/6; n3 = 3*n/6; n4 = 4*n/6; n5 = 5*n/6; n6 = n;
c0 = n0+1:n1; c1 = n1+1:n2; c2 = n2+1:n3; c3 = n3+1:n4; c4 = n4+1:n5; c5 = n5+1:n6; 
t18 = B(r0,c0)+B(r2,c2)/8;
t19 = B(r1,c0)-B(r2,c3)/8;
t20 = B(r0,c3)+B(r1,c2);
r21 = B(r0,c1)/8;
t21 = t19+r21;
t22 = B(r0,c4)+B(r2,c1)/8;
r23 = B(r1,c1)/8;
t23 = B(r2,c0)-r23;
t25 = B(r2,c2)-B(r2,c5);
t27 = t22+t18+B(r1,c3)/8;
t28 = B(r1,c4)-B(r1,c2)/8;
t29 = B(r2,c1)-B(r2,c3);
t30 = B(r0,c0)+B(r0,c3)/8;
r31 = B(r1,c5)/8;
t31 = B(r2,c4)-r31;
t33 = B(r0,c4)-B(r0,c5)/8;
t34 = t21+B(r2,c5)/8-B(r1,c4);
t35 = t23+(B(r0,c5)+t20)/8;
t36 = B(r1,c0)+r23;
t38 = t27-t28;
t39 = t31+t21;
r40 = t20/8;
t40 = r40-t18;
t41 = t33-t30;
t42 = -t27-t28;
r43 = B(r0,c2)/8;
t43 = t34-r43;
t45 = t36-r31;
t46 = t34+r43;
t47 = t36+r31;
t48 = t35-B(r2,c4);
t49 = t18+r40;
t50 = B(r2,c4)+t35;
t53 = B(r1,c4)+(B(r0,c3)+B(r2,c3))/8;
t54 = t30+t33;
t56 = t22-r31-r43;
t57 = t31+r21-t19;
r58 = (t25-t29)/8;
r59 = (t25+t29)/8;
oB1 = t53+t56-r23-B(r2,c0);
oB2 = t56-t23-t53;
oB8 = t40-t57;
oB9 = t39-t49;
oB10 = t57+t40;
oB11 = t39+t49;
oB17 = t50-r59;
oB19 = t50+r59;
oB22 = t42-t47;
oB23 = t38-t45;
oB24 = t41-t43;
oB26 = t43+t41;
oB28 = t48-r58;
oB30 = r58+t48;
oB32 = t38+t45;
oB33 = t47+t42;
oB38 = t46+t54;
oB39 = t46-t54;
v40 = oB8-oB33;
v41 = v40-oB28;
v42 = oB32-oB2;
v43 = oB9+oB23;
v44 = oB11-oB19;
v45 = v43+v44;
v46 = oB24-v41;
oB7 = v42-oB30;
oB3 = v46-oB22;
oB5 = oB17+oB39+v41;
v51 = oB26-v45;
v52 = oB1-oB10;
v54 = oB9-oB5;
v55 = oB38+oB2-v45;
oB0 = v44-oB38;
v57 = -oB7-v51;
v58 = oB11+oB3;
oB6 = oB26-v43;
oB29 = oB11-v42-v54;
oB25 = oB5+v55;
oB12 = oB2-v54;
oB15 = oB8-oB38-v51;
oB34 = v40-v52+oB0;
oB18 = oB8+v57;
oB36 = v57+oB3;
oB31 = oB1+oB33-oB6;
oB4 = v41-v52;
oB16 = oB1-oB9+oB39;
oB21 = oB8-v55;
oB13 = -v41;
oB20 = -v46-v52-oB11;
oB27 = oB22+oB7-oB10;
oB14 = oB7-v58;
oB35 = v58-oB9-v42;
oB37 = oB2-oB10+oB17;

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

b2 = iC31-iC27;
v54 = iC29+iC24;
v50 = iC21+iC17;
b4 = iC13+iC17;
b6 = iC4-iC16;
b9 = iC9-iC30;
v49 = iC23+iC19;
v51 = -v54;
b20 = v51-iC12-b6;
b27 = b9-iC14+v49;
v43 = iC34+b20;
b28 = b27-iC39;
b29 = iC15+b2+v43;
b30 = b28-b4-iC21;
z38 = iC38+iC32+v43;
z37 = iC37-b20;
z35 = iC35+b27;
z33 = iC33+iC36-b28;
z26 = iC26+b2-iC30;
z25 = iC25+iC28-v51;
z22 = iC22+v50-iC18;
z20 = iC20+iC16+v49;
z11 = iC11+iC29-iC32+iC12+v49;
z10 = iC10+iC31+iC34+iC15+v50;
z8 = iC8+iC36+iC28-b30;
z7 = iC7-iC32-iC18-b29;
z6 = iC6+iC36+iC14-iC19+b2;
z5 = iC5+v54-iC39-b4;
z3 = iC3+iC24+iC18+b30;
z2 = iC2+iC21-iC30+b29;
z1 = iC1+b9-iC27-iC16-iC39;
z0 = iC0+b6-iC23-iC34-iC28;
t49 = z38+z37;
b35 = z37-z38;
b36 = z25-z26;
b37 = -z26-z25;
b38 = z22-z33;
b39 = -z33-z22;
b40 = z20-z35;
b41 = -z35-z20;
t55 = z10-z8;
t39 = z10+z8;
b42 = -z7-z6;
b44 = z6-z7-z5;
t48 = z3-z1;
t38 = z3+z1;
b45 = z0-z2;
b46 = -z2-z0;
b47 = -z11-t55;
b48 = t55-z11;
t46 = -b46;
b50 = t49-b36;
b51 = b38+b41+t48;
b52 = b35-z7;
b53 = b36+t49+t46;
oC16 = (t46+t48)/8;
t28 = b35+b37;
b54 = z6+z5+b37;
b55 = -z5-b42;
t29 = b40-b39;
t24 = b40+b39;
b57 = t39-z11;
oC12 = (-b48-b55)/8;
b62 = t38-b45+t28;
oC13 = z11+t39-t29;
b63 = b45+t38+t29;
oC9 = b50-b44;
oC8 = b50-b48;
oC15 = b52+b54;
oC5 = b55+t24;
oC1 = b57+t24;
oC11 = t48+b57-b53;
oC7 = t48+b55+b53;
oC3 = b48-t46-b51;
oC2 = b46+b44+b51;
oC10 = (b54-b52+oC13)/8;
oC4 = (b41-b38-b47-oC15)/8;
oC17 = b47+b62;
oC0 = (-t29-b62)/8;
oC14 = b42-z5+b63;
oC6 = (b63-t28)/8;

C = [ oC0 oC1 oC2 oC3 oC4 oC5 ; oC6 oC7 oC8 oC9 oC10 oC11 ; oC12 oC13 oC14 oC15 oC16 oC17 ] ;
  end
end
end
