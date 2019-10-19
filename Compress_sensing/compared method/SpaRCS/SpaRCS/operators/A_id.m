%"Identity" operator for matrices.  Just convert the matrix into a vector

function b = A_id(X)

b = reshape(X,[],1);