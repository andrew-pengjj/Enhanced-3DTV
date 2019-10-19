function [fmax,S,thmax] = maxFMeasure(E,S0)
fmax = 0;
Vth = linspace(min(abs(E(:))),max(abs(E(:))),10000);
for idx = 1:length(Vth)
    th = Vth(idx);
    Stmp = abs(E) > th;
    rec = sum(Stmp(:)&S0(:))/(sum(S0(:))); 
    pre = (sum(Stmp(:)&S0(:))+1)/(sum(Stmp(:))+1);
    ftmp = 2*rec*pre/(pre+rec);
    if ftmp > fmax
        fmax = ftmp;
        S = Stmp;
        thmax =th;
    end
end

end