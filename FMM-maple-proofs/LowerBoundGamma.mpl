
# Proof of the Holder lower bound on the growth factor of Fast Matrix Multiplication in Strassen's orbit
restart;
with(LinearAlgebra);
with(Student[VectorCalculus]);
with(PolynomialTools);
HolderP := proc(A, p) local i, s; s := 0; for i to RowDimension(A) do s := s + MatrixNorm(A[i], Frobenius, conjugate = false)^p; end do; s^(1/p); end proc;
# L represents the left linear pre-additions performed by the original Strassen's algorithm on left-hand side of A*B
# K represents all the transformation of any L, R or P matrix within Strassen's orbit
L := Matrix(7, 4, [[1, 0, 0, 1], [0, 1, 0, -1], [-1, 0, 1, 0], [1, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 1]]);
W := <<r | x>, <0 | 1/r>>;
V := <<s | y>, <0 | 1/s>>;
K := KroneckerProduct(W, V);
L, W, V, K;
# It is easier to study the Holder norm to the power -1/z, and use a nicer form
holderexp := -z;
HNlk := HolderP(L . K, 1/holderexp);
G := ((r^2 + x^2)*(s^2 + y^2) + (2*x*y + 1/(r*s))/(r*s))^(-1/(2*z)) + (r*s)^(1/z) + (s^2 + (y + 1/s)^2)^(-1/(2*z))*((r^2 + x^2)^(-1/(2*z)) + r^(1/z)) + (r^2 + (x - 1/r)^2)^(-1/(2*z))*((s^2 + y^2)^(-1/(2*z)) + s^(1/z)) + ((r^2 + x^2)*(s^2 + y^2))^(-1/(2*z));
E := simplify(G^holderexp, symbolic);
simplify(E - HNlk, symbolic);
# We compute the gradient of the Holder norm, and check its root
fx := diff(E, x);
fy := diff(E, y);
fr := diff(E, r);
fs := diff(E, s);
gradE := [fx, fy, fr, fs];
explminpoint := simplify(subs({r = root[4](3/4)}, subs({s = r}, subs({x = (2*r^3)/3, y = -(2*s^3)/3}, [r, s, x, y]))));
subminpoint := solve([r, s, x, y] - explminpoint);
map(simplify, subs(subminpoint, gradE));
# We now compute the Hessian at that point, then the associated characteristic polynomial and eigenvalues
H := map(x -> simplify(radnormal(x, rationalized)), Hessian(E, [r, s, x, y] = explminpoint));
charP := FromCoefficientVector(map(x -> simplify(x, symbolic), CoefficientVector(simplify(expand(CharacteristicPolynomial(H, X)), radical), X)), X);
solutions := map(x -> simplify(x, symbolic), [solve(charP)]);
# We now look for a nicer form of the eigenvalues for z>0
tau := 3^(1/(2*z));
nu := 3^((z + 1)/(2*z));
lambda := 2^(1/(2*z));
d1 := 24*(1 - 16*z)*tau*z*lambda + (1344*z^2 - 384*z + 39)*tau^2 + 48*lambda^2*z^2;
d2 := 72*(272*z - 237)*tau*z*lambda + (12096*z^2 - 20736*z + 8991)*tau^2 + 8112*lambda^2*z^2;
b1 := (32*z - 11)*nu + 4*sqrt(3)*z*2^(1/(2*z));
b2 := (96*z - 63)*nu + 52*sqrt(3)*z*lambda;
n := (2^(1/(2*z)) + 6*3^(1/(2*z)))^(-z - 1)/(6*z);
eigs := [n*(b1 + sqrt(d1)), n*(b1 - sqrt(d1)), n/3*(b2 + sqrt(d2)), n/3*(b2 - sqrt(d2))];
sols := [subs(solutions[1], X), subs(solutions[2], X), subs(solutions[3], X), subs(solutions[4], X)];
map(y -> simplify(y, symbolic), eigs - sols);
# We check that the eigenvalues are real on the positive plane: the minimal values of d1 and d2 are positive
dd1 := diff(d1, z);
dd2 := diff(d2, z);
min1 := solve(dd1);
min2 := solve(dd2);
evalf([min1, min2]);
[[evalf(subs(z = min1 - 1/10, dd1)), evalf(subs(z = min1 + 1/10, dd1))], [evalf(subs(z = min2 - 1/10, dd2)), evalf(subs(z = min2 + 1/10, dd2))]];
evalf([subs(z = min1, d1), subs(z = min2, d2)]);
# We end by checking that the eigenvalues are positive after 0.5171, hence the Hessian is definite positive and the extremal point is a local minimum
sbd1 := b1^2 - d1;
sbd2 := b2^2 - d2;
bd1 := diff(sbd1, z);
bd2 := diff(sbd2, z);
bin1 := solve(bd1);
bin2 := solve(bd2);
evalf([bin1, bin2]);
[[evalf(subs(z = bin1 - 1/10, bd1)), evalf(subs(z = bin1 + 1/10, bd1))], [evalf(subs(z = bin2 - 1/10, bd2)), evalf(subs(z = bin2 + 1/10, bd2))]];
evalf([subs(z = bin1, sbd1), subs(z = bin2, sbd2)]);
evalf([[solve(b1^2 - d1)], [solve(b2^2 - d2)]]);
sin1 := solve(diff(b1, z));
sin2 := solve(diff(b2, z));
evalf([sin1, sin2]);
evalf([subs(z = 1, diff(b1, z)), subs(z = 1, diff(b2, z))]);
evalf(subs(z = 0.5171, [b1, b2]));
;
evalf(subs(z = 0.5171 + rand(), eigs));
# The extremal point is a minimum, we search for the limit of the combined Holder norms at infinity
MinPoint := simplify(expand(simplify(subs(subminpoint, E))));
LowerBound1 := limit(MinPoint^3*7^(1 + 3*z), z = infinity);
LowerBound2 := limit(MinPoint^2*subs(z = -2*z - 1, MinPoint), z = infinity);
evalf([LowerBound1, LowerBound2]);
