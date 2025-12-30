function [r,w,Y] = f_prices(K,alpha,delta)

r = alpha*K^(alpha-1) - delta;
w = (1-alpha)*K^alpha;
Y = K^alpha;

end