function   [x, svp]  =  solve_Lp_w( y, lambda, p )

% Modified by Dr. xie yuan
% lambda here presents the weights vector
J     =   2;  %2
% tau is generalized thresholding vector
tau   =  (2*lambda.*(1-p)).^(1/(2-p)) + p*lambda.*(2*(1-p)*lambda).^((p-1)/(2-p));
x     =   zeros( size(y) );
% i0表示 thresholding 后非0的个数
i0    =   find( abs(y)>tau );
svp   =   length(i0);

if length(i0)>=1
    % lambda  =   lambda(i0);
    y0    =   y(i0);
    t      =   abs(y0);
    lambda0 = lambda(i0);
    for  j  =  1 : J
        t    =  abs(y0) - p*lambda0.*(t).^(p-1);
    end
    x(i0)   =  sign(y0).*t;
end