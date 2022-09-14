function N = comNorm(A,B)
    N = max(max(abs(A),[],'all'), max(abs(B),[],'all'));
end