
# Proof of the Frobenius norm minimal point of Fast Matrix Multiplication in Strassen's orbit
restart;
with(LinearAlgebra);
with(Student[VectorCalculus]);
with(PolynomialTools);
# L represents the left linear pre-additions performed by the original Strassen's algorithm on left-hand side of A*B
# K represents all the transformation of any L, R or P matrix within Strassen's orbit
L := Matrix(7, 4, [[1, 0, 0, 1], [0, 1, 0, -1], [-1, 0, 1, 0], [1, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 1]]);
W := <<r | x>, <0 | 1/r>>;
V := <<s | y>, <0 | 1/s>>;
K := KroneckerProduct(W, V);
L, W, V, K;
# It is easier to study the square of the Frobenius norm, first its gradient
Nlk := MatrixNorm(L . K, Frobenius, conjugate = false);
E := Nlk^2;
fx := diff(E, x);
fy := diff(E, y);
fr := diff(E, r);
fs := diff(E, s);
gradE := [fx, fy, fr, fs];
# We have found a real root of the gradient and now check that this point is indeed an extremal point
explminpoint := simplify(subs({r = root[4](3/4)}, subs({s = r}, subs({x = (2*r^3)/3, y = -(2*s^3)/3}, [r, s, x, y]))));
subminpoint := solve([r, s, x, y] - explminpoint);
map(simplify, subs(subminpoint, gradE));
# We now compute the Hessian at that point  (multiplying 9/(4) to simplify the following computations) , then the associated characteristic polynomial and eigenvalues
H := map(x -> simplify(radnormal(x, rationalized)), 9/4*Hessian(E, [r, s, x, y] = explminpoint));
charP := FromCoefficientVector(map(x -> simplify(x, symbolic), CoefficientVector(simplify(expand(CharacteristicPolynomial(H, X)), radical), X)), X);
eigs := map(x -> simplify(x, symbolic), [solve(charP)]);
# We end by checking that these eigenvalues are all positive, hence the Hessian is definite positive and the extremal point is a local minimum
map(x -> 0 < x, eigs);
evalf(%);
map(evalb, %);
