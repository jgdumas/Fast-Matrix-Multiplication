function C = FMM633_mul_6_3_3(A, B, nmin, peeling, level)
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
    C(1:mu,1:nu)=FMM633_mul_6_3_3(A(1:mu,1:ku),B(1:ku,1:nu),nmin, peeling, level);
    if (m > mu)
      % fprintf("# MM peel m: %d x %d x %d\n",m-mu,k,n)
      C(mu+1:m,1:n)=C(mu+1:m,1:n)+FMM633_mul(A(mu+1:m,1:k),B,nmin, peeling, level);
    end
    if (k > ku) && (mu > 0) && (nu > 0)
      % fprintf("# MM peel k: %d x %d x %d\n",mu,k-ku,nu)
      C(1:mu,1:nu)=C(1:mu,1:nu)+FMM633_mul(A(1:mu,ku+1:k),B(ku+1:k,1:nu),nmin, peeling, level);
    end
    if (n > nu) && (mu > 0)
      % fprintf("# MM peel n: %d x %d x %d\n",mu,k,n-nu)
      C(1:mu,nu+1:n)=C(1:mu,nu+1:n)+FMM633_mul(A(1:mu,1:k),B(1:k,nu+1:n),nmin, peeling, level);
    end
  else
    if l>=level, fprintf("# Core<6;3;3>[%d]: %d x %d x %d\n",l,m,k,n); end

[m,n] = size(A);
m0 = 0; m1 = 1*m/6; m2 = 2*m/6; m3 = 3*m/6; m4 = 4*m/6; m5 = 5*m/6; m6 = m;
 r0 = m0+1:m1; r1 = m1+1:m2; r2 = m2+1:m3; r3 = m3+1:m4; r4 = m4+1:m5; r5 = m5+1:m6;
n0 = 0; n1 = 1*n/3; n2 = 2*n/3; n3 = n;
 c0 = n0+1:n1; c1 = n1+1:n2; c2 = n2+1:n3;
t18 = A(r2,c0)-A(r4,c0);
t19 = A(r0,c1)+A(r2,c2);
oA2 = t18-A(r3,c1);
t21 = A(r1,c1)+A(r5,c1);
t22 = A(r2,c1)-A(r1,c0);
oA8 = A(r4,c2)+t19;
t24 = A(r1,c2)-A(r3,c2);
t26 = A(r0,c2)+t21;
t27 = A(r4,c1)-A(r3,c0);
t28 = A(r0,c0)+oA2-A(r0,c1);
t29 = A(r4,c0)+A(r5,c2);
oA1 = t21-A(r4,c1);
oA4 = A(r4,c0)-t24;
oA7 = oA2-t22;
oA16 = A(r3,c0)-A(r4,c0)-A(r1,c1);
oA17 = A(r2,c1)+A(r4,c2)-t28;
oA19 = A(r3,c2)-t27;
oA20 = A(r5,c1)-t24-t27;
oA21 = A(r2,c1)+oA8-A(r0,c0);
oA22 = t22-t28;
oA24 = t19+t29-A(r1,c2);
oA27 = -t21-t22-A(r5,c0);
oA28 = A(r5,c2)+oA8-A(r3,c2);
oA29 = A(r2,c2)+t29;
oA30 = -oA2-A(r1,c1)-A(r5,c0);
oA31 = A(r2,c1)+A(r4,c1)+A(r5,c0);
oA33 = t26-A(r4,c1)-A(r4,c2);
oA34 = A(r1,c0)+t18+t24;
oA35 = A(r0,c2)+A(r1,c1)+A(r2,c2);
oA36 = t19+t26;
oA38 = A(r1,c2)+t18-t22;
oA0 = A(r3,c2);
oA3 = A(r0,c1);
oA5 = A(r2,c2);
oA6 = A(r4,c1);
oA9 = A(r1,c1);
oA10 = A(r1,c0);
oA11 = A(r1,c2);
oA12 = A(r4,c0);
oA13 = A(r4,c2);
oA14 = A(r5,c1);
oA15 = A(r2,c1);
oA18 = A(r0,c0);
oA23 = A(r3,c0);
oA25 = A(r5,c2);
oA26 = A(r5,c0);
oA32 = A(r3,c1);
oA37 = A(r2,c0);
oA39 = A(r0,c2);

[m,n] = size(B);
m0 = 0; m1 = 1*m/3; m2 = 2*m/3; m3 = m;
 r0 = m0+1:m1; r1 = m1+1:m2; r2 = m2+1:m3;
n0 = 0; n1 = 1*n/3; n2 = 2*n/3; n3 = n;
 c0 = n0+1:n1; c1 = n1+1:n2; c2 = n2+1:n3;
t9 = B(r0,c2)-B(r2,c1);
t10 = B(r0,c1)+B(r1,c2);
oB22 = t9+B(r2,c0);
t12 = B(r0,c0)+B(r0,c1);
oB16 = B(r0,c2)+B(r1,c0);
oB9 = B(r1,c1)-t9;
oB13 = t10+B(r0,c0);
t16 = -B(r1,c1)-B(r2,c1);
oB14 = t10-B(r1,c1);
oB37 = B(r2,c2)-B(r0,c1);
oB39 = -B(r2,c2)-B(r0,c0);
oB11 = B(r2,c1)+t10;
t21 = B(r0,c2)-B(r1,c1);
oB0 = t9-B(r0,c1);
oB2 = t12+oB9;
oB5 = B(r0,c2)-t12;
oB6 = oB13-t9;
oB8 = B(r1,c2)+t21;
oB15 = t21-B(r0,c0);
oB17 = oB16-t12;
oB18 = B(r1,c0)+B(r1,c1);
oB19 = -B(r1,c0)-oB11;
oB20 = oB22-t10;
oB21 = t12-t16-B(r0,c2)*2-B(r2,c0);
oB24 = B(r2,c2)-t10;
oB25 = B(r0,c2)+oB37;
oB26 = oB9-oB39;
oB28 = -B(r1,c2)-oB16;
oB30 = t16-B(r1,c0);
oB31 = B(r0,c0)-oB16;
oB32 = B(r0,c2)+B(r2,c0);
oB33 = oB13-oB22;
oB34 = oB22-B(r0,c1);
oB35 = B(r1,c1)-oB22;
oB36 = oB14-B(r2,c2);
oB38 = -B(r2,c1)-B(r2,c2);
oB1 = B(r0,c0);
oB3 = B(r1,c1);
oB4 = B(r1,c2);
oB7 = B(r2,c1);
oB10 = B(r0,c1);
oB12 = B(r0,c2);
oB23 = B(r2,c0);
oB27 = B(r2,c2);
oB29 = B(r1,c0);

iC0 = FMM633_mul( oA0, oB0, nmin, peeling, level);
iC1 = FMM633_mul( oA1, oB1, nmin, peeling, level);
iC2 = FMM633_mul( oA2, oB2, nmin, peeling, level);
iC3 = FMM633_mul( oA3, oB3, nmin, peeling, level);
iC4 = FMM633_mul( oA4, oB4, nmin, peeling, level);
iC5 = FMM633_mul( oA5, oB5, nmin, peeling, level);
iC6 = FMM633_mul( oA6, oB6, nmin, peeling, level);
iC7 = FMM633_mul( oA7, oB7, nmin, peeling, level);
iC8 = FMM633_mul( oA8, oB8, nmin, peeling, level);
iC9 = FMM633_mul( oA9, oB9, nmin, peeling, level);
iC10 = FMM633_mul( oA10, oB10, nmin, peeling, level);
iC11 = FMM633_mul( oA11, oB11, nmin, peeling, level);
iC12 = FMM633_mul( oA12, oB12, nmin, peeling, level);
iC13 = FMM633_mul( oA13, oB13, nmin, peeling, level);
iC14 = FMM633_mul( oA14, oB14, nmin, peeling, level);
iC15 = FMM633_mul( oA15, oB15, nmin, peeling, level);
iC16 = FMM633_mul( oA16, oB16, nmin, peeling, level);
iC17 = FMM633_mul( oA17, oB17, nmin, peeling, level);
iC18 = FMM633_mul( oA18, oB18, nmin, peeling, level);
iC19 = FMM633_mul( oA19, oB19, nmin, peeling, level);
iC20 = FMM633_mul( oA20, oB20, nmin, peeling, level);
iC21 = FMM633_mul( oA21, oB21, nmin, peeling, level);
iC22 = FMM633_mul( oA22, oB22, nmin, peeling, level);
iC23 = FMM633_mul( oA23, oB23, nmin, peeling, level);
iC24 = FMM633_mul( oA24, oB24, nmin, peeling, level);
iC25 = FMM633_mul( oA25, oB25, nmin, peeling, level);
iC26 = FMM633_mul( oA26, oB26, nmin, peeling, level);
iC27 = FMM633_mul( oA27, oB27, nmin, peeling, level);
iC28 = FMM633_mul( oA28, oB28, nmin, peeling, level);
iC29 = FMM633_mul( oA29, oB29, nmin, peeling, level);
iC30 = FMM633_mul( oA30, oB30, nmin, peeling, level);
iC31 = FMM633_mul( oA31, oB31, nmin, peeling, level);
iC32 = FMM633_mul( oA32, oB32, nmin, peeling, level);
iC33 = FMM633_mul( oA33, oB33, nmin, peeling, level);
iC34 = FMM633_mul( oA34, oB34, nmin, peeling, level);
iC35 = FMM633_mul( oA35, oB35, nmin, peeling, level);
iC36 = FMM633_mul( oA36, oB36, nmin, peeling, level);
iC37 = FMM633_mul( oA37, oB37, nmin, peeling, level);
iC38 = FMM633_mul( oA38, oB38, nmin, peeling, level);
iC39 = FMM633_mul( oA39, oB39, nmin, peeling, level);
t27 = iC34-iC31;
oC0 = iC20-iC27-iC24-iC22;
t26 = iC38+iC18;
b7 = iC17-iC39;
b14 = iC34+iC22+iC10;
b17 = iC28+iC18+iC8;
t28 = iC17+iC5;
t30 = iC22-iC3;
oC11 = iC16-iC37-b7;
oC4 = iC31-iC27+iC10+iC1+b7;
b27 = iC21-iC2+t28;
oC17 = iC33-iC28-t27;
oC7 = iC6-iC23-iC19+iC0+t27;
b29 = t26-iC19;
oC6 = iC29-iC24+iC11-iC3+t26;
oC10 = iC28+iC17+iC13+b14;
oC3 = iC24-iC37+iC4+b14;
b30 = iC38-iC0+b17;
oC13 = iC34+iC21+iC15+b17;
oC5 = iC36+b29;
oC15 = iC27-iC35+iC14+t30+b29;
oC8 = iC30-iC35+b27;
oC9 = iC32-iC29-b27;
oC2 = iC7-iC18-t30+b27;
t19 = -b30;
oC12 = iC25+b30;
b31 = t19-iC23;
oC14 = iC37+iC29+iC12+t28+t19;
oC1 = iC26+iC21-b31;
oC16 = iC39-iC35+iC9+iC5+b31;

C = [ oC0 oC1 oC2 ; oC3 oC4 oC5 ; oC6 oC7 oC8 ; oC9 oC10 oC11 ; oC12 oC13 oC14 ; oC15 oC16 oC17 ] ;
  end
end
end
