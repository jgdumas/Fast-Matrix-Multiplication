function C = recursive(A, B, M, K, N)
%          Recursive matrix multiplication algorithm,
%          with MxK blocks times KxN blocks
[m,k] = size(A);
[k2,n] = size(B);
if (k2 ~= k)
   error('Incompatible matrix dimensions.')
end
if (m <= M) || (k <= K) || (n <= N)
  % fprintf("# MM Direct: %d x %d x %d\n",m,k,n)
  C = A*B;
else
  C = zeros(m,n);
  mu=M*floor(m/M);ku=K*floor(k/K);nu=N*floor(n/N);
  if (mu < m) || (ku < k) || (nu < n)
    % fprintf("# Core SubMatrix: %d x %d x %d\n",mu,ku,nu)
    C(1:mu,1:nu)=recursive(A(1:mu,1:ku),B(1:ku,1:nu), M, K, N);
    if (m > mu)
      % fprintf("# MM peel m: %d x %d x %d\n",m-mu,k,n)
      C(mu+1:m,1:n)=C(mu+1:m,1:n)+recursive(A(mu+1:m,1:k),B, M, K, N);
    end
    if (k > ku) && (mu > 0) && (nu > 0)
      % fprintf("# MM peel k: %d x %d x %d\n",mu,k-ku,nu)
      C(1:mu,1:nu)=C(1:mu,1:nu)+recursive(A(1:mu,ku+1:k),B(ku+1:k,1:nu), M, K, N);
    end
    if (n > nu) && (mu > 0)
      % fprintf("# MM peel n: %d x %d x %d\n",mu,k,n-nu)
      C(1:mu,nu+1:n)=C(1:mu,nu+1:n)+recursive(A(1:mu,1:k),B(1:k,nu+1:n), M, K, N);
    end
  else
    m=m/M;k=k/K;n=n/N;
    fprintf("# %d Block multiplications: %d x %d x %d\n",M*K*N,m,k,n)
    for i = 1:M
      for j = 1:K
        for l = 1:N
           C((i-1)*m+1:i*m,(l-1)*n+1:l*n) = C((i-1)*m+1:i*m,(l-1)*n+1:l*n) + recursive(A((i-1)*m+1:i*m,(j-1)*k+1:j*k), B((j-1)*k+1:j*k,(l-1)*n+1:l*n), M, K, N);
        end
      end
    end
  end
end
end
