function F = f_ReturnFn(aprime,a,z,K,alpha,delta,sigma)
% Action space: (a',a,z)

[r,w] = f_prices(K,alpha,delta);

cons = (1+r)*a + w*z - aprime;

F = -inf;

if cons>0
    F = cons^(1-sigma)/(1-sigma);
end

end %end function