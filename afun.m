
function Av = afun(A,At,v,flag)
if strcmp(flag,'notransp')
    Av = A(v);
elseif strcmp(flag,'transp')
    Av = At(v);
end
end