function  [X num]  =  Cube2Matrices( data, par )

    [row col bands] = size(data);
    win = par.win;
    step = par.step;
    num = 0;
    X = zeros(win^2, bands, floor((row-win)/step+1)*floor((row-win)/step+1));
    
    for i = 1:step:row-win+1
        for j = 1:step:col-win+1
            num = num+1;
            for k = 1:bands
                patch = data(i:i+win-1, j:j+win-1, k);
                X(:,k,num) = patch(:)';
            end
        end
    end

end