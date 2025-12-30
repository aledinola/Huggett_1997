function cons = f_consumption(aprime,a,z,K,alpha,delta)
% Action space: (a',a,z)

[r,w] = f_prices(K,alpha,delta);

cons = (1+r)*a + w*z - aprime;

end %end function