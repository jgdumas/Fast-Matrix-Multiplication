function C = FMM_3_6_3(A, B, nmin, peeling, level)
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
    C(1:mu,1:nu)=FMM_3_6_3(A(1:mu,1:ku),B(1:ku,1:nu),nmin, peeling, level);
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
    if l>=level, fprintf("# Core<3;6;3>[%d]: %d x %d x %d\n",l,m,k,n); end

[m,n] = size(A);
m0 = 0; m1 = 1*m/3; m2 = 2*m/3; m3 = m;
 r0 = m0+1:m1; r1 = m1+1:m2; r2 = m2+1:m3;
n0 = 0; n1 = 1*n/6; n2 = 2*n/6; n3 = 3*n/6; n4 = 4*n/6; n5 = 5*n/6; n6 = n;
 c0 = n0+1:n1; c1 = n1+1:n2; c2 = n2+1:n3; c3 = n3+1:n4; c4 = n4+1:n5; c5 = n5+1:n6;
t18 = A(r1,c5)+A(r0,c1);
t19 = A(r1,c3)-A(r0,c2);
t20 = A(r2,c5)+A(r2,c1);
r21 = A(r2,c0)/8;
t21 = t18+r21;
t22 = A(r2,c3)+A(r2,c2);
r24 = A(r0,c0)/8;
t24 = A(r2,c5)-r24;
t25 = A(r1,c2)-A(r0,c3);
t26 = A(r1,c1)+t24;
r27 = A(r1,c0)/8;
t27 = A(r1,c5)+r27;
t28 = A(r2,c3)+A(r1,c2);
t30 = t19-t21;
r32 = A(r2,c4)/8;
t32 = A(r0,c3)-r32;
t33 = A(r1,c1)+A(r0,c5)-r21;
t35 = t26+t28;
t37 = t20-t22;
t43 = t26-t32;
r44 = (A(r0,c4)-A(r1,c4))/8;
t44 = t20-r44;
t45 = A(r1,c3)-t27;
t46 = t21+t25;
t47 = -t22-t20;
t48 = A(r0,c1)-A(r2,c1)+r27-A(r1,c3);
t49 = A(r0,c2)+A(r2,c2)+t27;
t50 = t18-t25-r21;
t51 = A(r0,c5)-r24+t28-r32;
r52 = (A(r0,c4)+A(r1,c4))/8;
t52 = A(r2,c1)-A(r2,c5)+r52;
t53 = t30-r32;
t54 = t22-r52;
t55 = t19+t33;
t56 = t30+r32;
oA1 = r27-A(r1,c1)-A(r1,c5)+A(r2,c2)+t24+t32-A(r0,c2);
oA2 = t43-t49;
oA3 = t49+t43;
oA4 = t54+t55;
oA5 = A(r2,c3)+t19-t33-r44-A(r2,c2);
oA7 = t55-t54;
oA8 = t44+t46;
oA9 = t44-t46;
oA10 = t50+t52;
oA11 = t52-t50;
oA12 = t51-t48;
oA14 = t48+t51;
oA16 = t47-t53;
oA18 = -t47-t53;
oA28 = t37-t56;
oA30 = -t37-t56;
oA36 = t35-t45-r44;
oA37 = t35+t45-r52;
v40 = oA4+oA11;
v41 = oA14-oA1;
v42 = oA8-oA5;
v43 = oA7-oA10;
oA0 = v40-oA12;
oA13 = v42-oA3;
oA6 = oA9+v41;
v49 = v43+oA18;
oA32 = oA37-oA12-oA2;
oA35 = oA36-oA3-oA14;
oA34 = oA10+oA37-oA4;
oA21 = oA2+oA8-v49;
oA24 = oA4+oA28-oA13;
oA15 = oA2-v43;
oA23 = oA9+oA12+oA16;
oA22 = oA3-oA7-oA18;
oA19 = oA16+v40-v41;
oA38 = oA7+oA11+oA32;
oA17 = v42-v49;
oA27 = oA7-oA14+oA30;
oA39 = oA35-oA5-oA9;
oA29 = oA28+v40-v42;
oA20 = oA1+oA4+oA16;
oA31 = v41-v43-oA30;
oA25 = oA28+oA0-oA8;
oA26 = -oA2-oA9-oA30;
oA33 = oA36-oA6-oA8;

[m,n] = size(B);
m0 = 0; m1 = 1*m/6; m2 = 2*m/6; m3 = 3*m/6; m4 = 4*m/6; m5 = 5*m/6; m6 = m;
 r0 = m0+1:m1; r1 = m1+1:m2; r2 = m2+1:m3; r3 = m3+1:m4; r4 = m4+1:m5; r5 = m5+1:m6;
n0 = 0; n1 = 1*n/3; n2 = 2*n/3; n3 = n;
 c0 = n0+1:n1; c1 = n1+1:n2; c2 = n2+1:n3;
r20 = B(r2,c2)/8;
t20 = B(r0,c0)+r20;
r21 = B(r3,c2)/8;
t21 = B(r0,c1)-r21;
r23 = B(r5,c1)/8;
r24 = B(r1,c2)/8;
t24 = B(r4,c2)-r24;
r25 = B(r1,c0)/8;
t25 = B(r4,c2)-r23+r25;
r26 = (B(r2,c1)+B(r3,c0))/8;
r27 = (B(r1,c1)-B(r5,c0))/8;
t27 = B(r0,c2)-r26+r27;
t29 = B(r4,c1)-B(r5,c0)/8;
t30 = B(r4,c0)-r25;
r31 = B(r2,c0)/8;
t31 = B(r0,c0)-r31;
r32 = B(r5,c2)/8;
t32 = B(r0,c1)+r32;
t33 = t29-r32;
t34 = t20+t25;
r35 = B(r3,c1)/8;
t36 = t32+r27;
t37 = t24+r35-t31;
t40 = t25-t20;
r41 = (B(r1,c2)+B(r3,c1))/8;
r43 = B(r2,c1)/8;
t47 = t21-r26;
t48 = B(r0,c2)-t30-r35;
t50 = t30-t31;
t51 = t21+r26;
t52 = t21-t33+B(r3,c0)/8;
t53 = -B(r4,c2)-r21-r24;
t54 = t32-r27;
t55 = -t20-B(r4,c0)-r41-r23;
t56 = B(r0,c1)+B(r1,c1)/8+B(r4,c1)-r43;
r57 = (B(r2,c2)+B(r5,c2))/8;
t57 = r57-t27;
t58 = t24-r21;
t59 = B(r4,c2)+r31+B(r0,c0)+r41;
t60 = t27+r57;
oB4 = t37-t54;
oB5 = t54+t37;
oB6 = t36-t59;
oB7 = t36+t59;
oB8 = t51-t34;
oB9 = t47+t40;
oB10 = t40-t47;
oB11 = t34+t51;
oB12 = r57-t29-t48+r43;
oB13 = t48+r43-r20-t33;
oB16 = t57-t58;
oB18 = t53-t60;
oB22 = t55-t56;
oB23 = -t55-t56;
oB24 = t50-t52;
oB25 = t50+t52;
oB29 = t58+t57;
oB31 = -t53-t60;
v40 = oB10-oB4;
v41 = oB8-oB6;
oB0 = oB12-oB23-oB25;
oB3 = oB24+oB13-oB22;
oB15 = oB0+v41;
v47 = oB7-oB11;
v48 = oB5-oB9;
oB1 = v40-oB13;
oB21 = oB5+oB8-oB25;
oB33 = oB6+oB31-oB1;
oB20 = oB4-oB11-oB24;
oB35 = oB29+oB3-oB5;
oB2 = oB12-v48;
oB36 = oB18+oB3-oB8;
oB27 = oB7-oB10+oB22;
oB37 = oB12-oB16-oB4;
oB26 = oB6+oB9+oB23;
oB17 = v48-oB16+v40;
oB39 = oB9+oB16-oB1;
oB32 = oB11+oB12-oB29;
oB30 = v48-oB29-v47;
oB14 = v47-oB3;
oB19 = v41-v47-oB18;
oB28 = v40+v41-oB31;
oB38 = oB7+oB18-oB15;
oB34 = oB10-oB31+oB15;

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
b3 = iC22-iC35;
v51 = iC28-iC4;
b6 = iC24-iC38;
b8 = iC24-iC13-iC11;
v49 = iC22+iC3;
b11 = iC6+iC36+iC33+iC8;
b15 = iC17-iC10;
v54 = iC1-iC27;
b16 = iC27-iC39;
b20 = iC7-iC14-iC30;
b25 = v51-iC34-iC0;
b33 = iC29-v49+b15+b8;
b34 = iC13-iC9+b3-b20;
b36 = iC10-iC31-v54-b11;
b41 = iC8-iC21-iC2-b25+b34;
z26 = iC36+b41-iC26-b6;
z25 = iC25+iC37+b16+b41;
oC1 = iC23+iC34-iC21+iC20-iC33+iC32+b3;
z19 = iC19+iC20-iC14+v51+b8-b11;
z18 = iC18+v49-iC7+b36;
z16 = iC4-v54-b20-b33-iC16-iC20;
z15 = iC15+b6-iC13+b25+b36;
z12 = iC12+iC2-iC37+iC32-b33;
z5 = iC5+b15-b16-b34;
oC3 = z26-z25;
b49 = z18-z19;
b50 = z19+z18-z15;
t15 = z12+z5;
b51 = z12-z5-z16;
oC5 = t15+z15;
oC8 = z16-b49;
oC4 = t15-z15-oC3;
oC7 = z25+z26-oC1+oC5;
oC2 = b51-b50;
oC0 = b51-oC1+b50;
oC6 = b49+z16-oC7;

C = [ oC0 oC1 oC2 ; oC3 oC4 oC5 ; oC6 oC7 oC8 ] ;
  end
end
end
