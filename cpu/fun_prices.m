function [r,wage] = fun_prices(KLbar,alpha,delta)


r     = alpha*KLbar.^(alpha-1)-delta;
wage  = (1-alpha)*KLbar.^alpha;


end % end function "fun_prices"

