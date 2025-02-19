function C = FMM_3_6_3(A, B, nmin, peeling, level)
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
t18 = A(r0,c1)+A(r1,c5);
t19 = A(r1,c3)+A(r2,c3);
r20 = A(r0,c4)/8;
t20 = A(r2,c5)-r20;
t21 = A(r1,c3)-A(r2,c3);
r22 = A(r2,c4)/8;
t22 = A(r2,c1)+r22;
r23 = A(r1,c4)/8;
t23 = A(r1,c1)-r23;
r24 = A(r1,c0)/8;
t24 = A(r1,c2)-r24;
r25 = A(r0,c0)/8;
t25 = t19-r25;
t26 = t18+A(r0,c2);
t27 = A(r2,c1)+r23;
r28 = A(r2,c0)/8;
t28 = A(r2,c2)+r28;
t29 = t21-t26;
t31 = A(r1,c2)-A(r0,c3)+r28;
t32 = A(r2,c2)-r28;
t33 = t20+t25;
t34 = A(r2,c1)-r22;
t36 = t24-A(r1,c5);
t37 = t19+t32;
t38 = A(r2,c5)+t28;
t39 = A(r0,c5)+t23;
t41 = t20+t31;
t42 = t20-t31;
t43 = t23+t33;
t45 = t28-A(r2,c5);
t46 = t34+A(r0,c5)+t24;
t48 = A(r0,c1)-t25;
t49 = t22-t29;
t50 = t22-A(r1,c2)-A(r0,c5)-r24;
t51 = t18+t27;
t52 = t18-t27;
t53 = A(r0,c2)+r20;
t54 = A(r0,c1)-t21-r25;
t55 = t22+t29;
oA4 = t37+t39-t53;
oA5 = t19-t32-t39-t53;
oA8 = t41+t51;
oA9 = t42-t52;
oA10 = t51-t41;
oA11 = -t42-t52;
oA12 = t46-t48;
oA13 = t46+t48;
oA14 = t54-t50;
oA15 = t50+t54;
oA16 = t26-t34-t37-A(r2,c5);
oA18 = t38+t49;
oA19 = t55+t45;
oA24 = t43-t36;
oA30 = t45-t55;
oA31 = t49-t38;
oA37 = t36+t43;
oA39 = A(r1,c5)+t23+t24-t33;
v40 = oA11-oA19;
v42 = oA16+oA4+v40;
v43 = oA31-v42;
v44 = v43+oA30;
v45 = oA8-oA13;
oA3 = v45-oA5;
v47 = oA12-oA11;
v48 = oA15-oA37;
v49 = oA18-oA8;
oA7 = oA10-v44;
v51 = oA9+oA39;
oA1 = oA14-v42;
oA20 = oA14-v40;
oA23 = oA9+oA12+oA16;
oA6 = oA9+v42;
oA33 = oA39+oA1-oA13;
oA2 = oA15-v44;
oA27 = oA10-oA14-v43;
oA28 = oA13+oA24-oA4;
oA36 = oA14+v45+v51;
oA17 = v44-v49-oA5;
oA32 = v44-v48-oA12;
oA29 = oA11+oA24-oA3;
oA38 = oA10-v47-v48;
oA34 = oA10+oA37-oA4;
oA22 = oA3-oA7-oA18;
oA25 = oA24-v45-v47;
oA35 = oA5+v51;
oA26 = v43-oA9-oA15;
oA21 = oA15-v49;
oA0 = oA4-v47;

[m,n] = size(B);
m0 = 0; m1 = 1*m/6; m2 = 2*m/6; m3 = 3*m/6; m4 = 4*m/6; m5 = 5*m/6; m6 = m;
r0 = m0+1:m1; r1 = m1+1:m2; r2 = m2+1:m3; r3 = m3+1:m4; r4 = m4+1:m5; r5 = m5+1:m6; 
n0 = 0; n1 = 1*n/3; n2 = 2*n/3; n3 = n;
c0 = n0+1:n1; c1 = n1+1:n2; c2 = n2+1:n3; 
t18 = B(r4,c0)+B(r2,c2)/8;
t20 = B(r0,c0)-B(r0,c1);
r21 = B(r3,c0)/8;
t21 = B(r4,c2)+r21;
t22 = B(r2,c2)-B(r3,c2);
t23 = B(r1,c2)-B(r5,c2);
t24 = B(r4,c1)+t18;
t26 = B(r2,c1)+B(r5,c1);
t27 = t21-B(r1,c1)/8;
r28 = (B(r2,c1)+B(r5,c0))/8;
t28 = B(r0,c2)-r28;
t30 = B(r0,c2)+r28;
t31 = B(r2,c1)-B(r5,c1);
t32 = B(r4,c1)-t18;
r33 = B(r1,c2)/8;
t33 = t24+r33;
r34 = (B(r1,c1)-B(r3,c1))/8;
t34 = r34-t20;
t35 = B(r0,c0)+B(r0,c1);
t41 = t27-t28;
r43 = (B(r1,c0)-B(r3,c1))/8;
t43 = t30+r43;
t44 = -t20-(B(r2,c2)+B(r3,c2))/8;
t46 = -t28-r43;
r47 = B(r5,c2)/8;
t47 = t32-r47;
t48 = B(r4,c2)+B(r5,c0)/8+B(r2,c0)/8;
t49 = B(r4,c2)-r21+(B(r1,c0)-t26)/8;
t50 = t35+t22/8;
t51 = t24-r47;
t52 = t21+t31/8+B(r1,c0)/8;
t53 = t27+t30;
t54 = t20+r34;
r55 = t23/8;
r56 = t26/8;
r57 = (t22-t23)/8;
r58 = (t22+t23)/8;
oB5 = t48-t54-r55;
oB6 = t34-t48-r55;
oB8 = t44-t49;
oB9 = t49+t44;
oB10 = t52-t50;
oB11 = t50+t52;
oB12 = t46-t47;
oB13 = t43-t51;
oB14 = -t43-t51;
oB15 = t47+t46;
oB17 = t53-r58;
oB19 = t53+r58;
oB22 = (t31-B(r1,c1)-B(r3,c1))/8-t33-t35;
oB23 = r33-t32-t34+r56;
oB29 = t41+r57;
oB31 = t41-r57;
oB33 = r56-t33+t34;
oB35 = t33+t54+r56;
oB3 = oB35+oB5-oB29;
v43 = oB6+oB15;
oB1 = oB31+oB6-oB33;
v45 = oB3+oB22;
v46 = oB5-oB17;
oB0 = v43-oB8;
oB39 = oB13+v46;
v51 = oB12-oB23;
v52 = oB10-oB11-v45;
v53 = oB9+oB12;
oB36 = -oB19-oB14-oB6;
oB27 = oB14-v52;
oB7 = oB11+oB14+oB3;
oB4 = oB10-oB13-oB1;
oB21 = oB5+v43-v51;
oB25 = v51-oB0;
oB18 = oB8-oB3+oB36;
oB28 = oB8+oB13-oB33;
oB20 = v52-oB1;
oB24 = v45-oB13;
oB16 = oB1+oB39-oB9;
oB2 = v53-oB5;
oB30 = -oB9-oB14-oB35;
oB37 = v53-oB10-v46;
oB38 = oB11-oB19-oB0;
oB32 = oB11+oB12-oB29;
oB34 = oB10+oB15-oB31;
oB26 = oB6+oB9+oB23;

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

b3 = iC11-iC29+iC32+iC12;
b8 = iC33+iC13+iC8+iC6;
b17 = iC17-iC31+iC18;
b19 = iC8-iC21+iC25+iC5;
v51 = iC29-iC19;
v50 = iC24-iC36;
b24 = iC14-iC35-iC5-iC3+iC36;
v49 = iC28-iC4;
b32 = b17-iC1-b8;
b35 = iC0+iC34-v49+b19;
v44 = iC7+b3;
b40 = iC2+v44-b32;
z27 = b40-iC27-iC37-v50;
z26 = iC26+iC25-iC39-iC38+b40;
z23 = iC23+iC6-iC0-iC12-iC25+v51+b24;
z22 = iC22+iC13-iC24+iC3+b3+b35;
z20 = iC20+iC19-iC11-iC14+v50+v49-b8;
oC8 = iC28+iC30-b17-v51-iC16;
z15 = iC15+b19-iC38-iC18+v44;
z10 = iC10+iC2-iC37-iC17-b35;
z9 = iC9+iC39-iC30-b24+b32;
oC3 = z27-z26;
t15 = z23+z20;
t13 = z15+z10;
b47 = z22-z15+z10;
b48 = z23-z20+oC8+z9;
oC1 = z22+t15;
b49 = t13-z9;
oC5 = z9+t13;
oC2 = b49-oC8;
oC4 = z26+z27-b49;
oC0 = b47-b48;
oC6 = b47-oC3+b48;
oC7 = z22-t15-oC4;

C = [ oC0 oC1 oC2 ; oC3 oC4 oC5 ; oC6 oC7 oC8 ] ;
  end
end
end
