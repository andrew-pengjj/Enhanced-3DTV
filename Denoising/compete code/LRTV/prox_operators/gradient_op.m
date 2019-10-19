function [dx, dy] = gradient_op(I, weights_dx, weights_dy)

dx = [I(2:end, :)-I(1:end-1, :) ; zeros(1, size(I, 2))];
dy = [I(:, 2:end)-I(:, 1:end-1) , zeros(size(I, 1), 1)];

if nargin>1
    dx = dx .* weights_dx;
    dy = dy .* weights_dy;
end

end