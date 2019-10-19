function mse=MSE(sig0,sig1)

mse = mean( ( sig0(:) - sig1(:) ).^2 );


end