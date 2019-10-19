function cube = Matrices2Cube(X, data, par)

    [row col bands] = size(data);
    win = par.win;
    step = par.step;
    cube = zeros(size(data));
    W = zeros(size(data));
    num = 0;
    
    for i = 1:step:row-win+1
        for j = 1:step:col-win+1
            num = num+1;
            for k = 1:bands
                cube(i:i+win-1, j:j+win-1, k) = cube(i:i+win-1, j:j+win-1, k)+reshape(X(:, k, num), [win win]);
                W(i:i+win-1, j:j+win-1, k) = W(i:i+win-1, j:j+win-1, k)+1;
            end
        end
    end
    
    cube = cube./W;
    
end