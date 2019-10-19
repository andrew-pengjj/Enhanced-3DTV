function y = TV_norm(u)

    
        [dx, dy] = gradient_op(u);
    
temp = sqrt(abs(dx).^2 + abs(dy).^2);
y = sum(temp(:));

end