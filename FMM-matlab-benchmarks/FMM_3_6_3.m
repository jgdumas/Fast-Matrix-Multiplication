function C = FMM_3_6_3(A, B, nmin)
[m,k] = size(A);
[k,n] = size(B);
if (m<=nmin)||(k<=nmin)||(n<=nmin)||(m<=3)||(k<=6)||(n<=3)
  fprintf("# MM Direct: %d x %d x %d\n",m,k,n)
  C = A*B;
else
  C = zeros(m,n);
  mu=nmin;ku=nmin;nu=nmin;l=0;
  while (mu < m) && (ku < k) && (nu < n)
    l=l+1; mu=mu*3; ku=ku*6; nu=nu*3;
  end
  l=l-1;mu=m-mod(m,nmin*3^l); ku=k-mod(k,nmin*6^l); nu=n-mod(n,nmin*3^l);
  if (mu < m) || (ku < k) || (nu < n)
    fprintf("# Core SubMatrix: %d x %d x %d\n",mu,ku,nu)
    C(1:mu,1:nu)=FMM_3_6_3(A(1:mu,1:ku),B(1:ku,1:nu),nmin);
    if (m > mu)
      % fprintf("# MM peel m: %d x %d x %d\n",m-mu,k,n)
      C(mu+1:m,1:n)=C(mu+1:m,1:n)+FMM(A(mu+1:m,1:k),B);
    end
    if (k > ku) && (mu > 0) && (nu > 0)
      % fprintf("# MM peel k: %d x %d x %d\n",mu,k-ku,nu)
      C(1:mu,1:nu)=C(1:mu,1:nu)+FMM(A(1:mu,ku+1:k),B(ku+1:k,1:nu));
    end
    if (n > nu) && (mu > 0)
      % fprintf("# MM peel n: %d x %d x %d\n",mu,k,n-nu)
      C(1:mu,nu+1:n)=C(1:mu,nu+1:n)+FMM(A(1:mu,1:k),B(1:k,nu+1:n));
    end
  else
    fprintf("# Core<3;6;3>: %d x %d x %d\n",m,k,n)
[m,n] = size(A);
m0 = 0; m1 = 1*m/3; m2 = 2*m/3; m3 = m;
r0 = m0+1:m1; r1 = m1+1:m2; r2 = m2+1:m3; 
n0 = 0; n1 = 1*n/6; n2 = 2*n/6; n3 = 3*n/6; n4 = 4*n/6; n5 = 5*n/6; n6 = n;
c0 = n0+1:n1; c1 = n1+1:n2; c2 = n2+1:n3; c3 = n3+1:n4; c4 = n4+1:n5; c5 = n5+1:n6; 
oA0 = A(r0,c0)+A(r0,c3);
oA1 = A(r0,c1)-A(r0,c3);
oA2 = A(r0,c2)-A(r0,c0);
oA3 = A(r0,c0)+A(r0,c1);
oA6 = A(r0,c2)+A(r0,c3);
oA4 = A(r0,c0);
oA5 = A(r0,c3);

[m,n] = size(B);
m0 = 0; m1 = 1*m/6; m2 = 2*m/6; m3 = 3*m/6; m4 = 4*m/6; m5 = 5*m/6; m6 = m;
r0 = m0+1:m1; r1 = m1+1:m2; r2 = m2+1:m3; r3 = m3+1:m4; r4 = m4+1:m5; r5 = m5+1:m6; 
n0 = 0; n1 = 1*n/3; n2 = 2*n/3; n3 = n;
c0 = n0+1:n1; c1 = n1+1:n2; c2 = n2+1:n3; 
oB0 = B(r0,c0)+B(r1,c0);
oB1 = B(r0,c2)+B(r1,c0);
oB2 = B(r0,c0)+B(r0,c1);
oB4 = B(r0,c1)-B(r1,c0);
oB5 = B(r0,c2)-B(r0,c0);
oB3 = B(r1,c0);
oB6 = B(r0,c0);

iC0 = FMM( oA0, oB0, nmin);
iC1 = FMM( oA1, oB1, nmin);
iC2 = FMM( oA2, oB2, nmin);
iC3 = FMM( oA3, oB3, nmin);
iC4 = FMM( oA4, oB4, nmin);
iC5 = FMM( oA5, oB5, nmin);
iC6 = FMM( oA6, oB6, nmin);
iC7 = FMM( oA7, oB7, nmin);
iC8 = FMM( oA8, oB8, nmin);
iC9 = FMM( oA9, oB9, nmin);
iC10 = FMM( oA10, oB10, nmin);
iC11 = FMM( oA11, oB11, nmin);
iC12 = FMM( oA12, oB12, nmin);
iC13 = FMM( oA13, oB13, nmin);
iC14 = FMM( oA14, oB14, nmin);
iC15 = FMM( oA15, oB15, nmin);
iC16 = FMM( oA16, oB16, nmin);
iC17 = FMM( oA17, oB17, nmin);
iC18 = FMM( oA18, oB18, nmin);
iC19 = FMM( oA19, oB19, nmin);
iC20 = FMM( oA20, oB20, nmin);
iC21 = FMM( oA21, oB21, nmin);
iC22 = FMM( oA22, oB22, nmin);
iC23 = FMM( oA23, oB23, nmin);
iC24 = FMM( oA24, oB24, nmin);
iC25 = FMM( oA25, oB25, nmin);
iC26 = FMM( oA26, oB26, nmin);
iC27 = FMM( oA27, oB27, nmin);
iC28 = FMM( oA28, oB28, nmin);
iC29 = FMM( oA29, oB29, nmin);
iC30 = FMM( oA30, oB30, nmin);
iC31 = FMM( oA31, oB31, nmin);
iC32 = FMM( oA32, oB32, nmin);
iC33 = FMM( oA33, oB33, nmin);
iC34 = FMM( oA34, oB34, nmin);
iC35 = FMM( oA35, oB35, nmin);
iC36 = FMM( oA36, oB36, nmin);
iC37 = FMM( oA37, oB37, nmin);
iC38 = FMM( oA38, oB38, nmin);
iC39 = FMM( oA39, oB39, nmin);
oC2 = iC6+iC5;
oC1 = iC4+iC3;
oC3 = iC4-iC6+iC2+iC0;
oC0 = iC5-iC3+iC1+iC0;

C = [ oC0 oC1 oC2 ; oC3 oC4 oC5 ; oC6 oC7 oC8 ] ;
  end
end
end
