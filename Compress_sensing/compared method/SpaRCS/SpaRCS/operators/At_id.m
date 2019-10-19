%"Identity" transpose operator for matrices.  Just reshape a vector to the
%right matrix dimensions

function X = At_id(b,rows,cols)

X = reshape(b,rows,cols);