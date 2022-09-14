function Y = ReLu_mat(X)
    [m,n] = size(X);
    Y = X;
    for i=1:m
        for j=1:n
            Y(i,j) = relu(X(i,j));
        end
    end
end