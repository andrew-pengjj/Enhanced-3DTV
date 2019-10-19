function x = hardthre(a, tau)
x         = a;
x(find(abs(a)<=tau)) =0;
end