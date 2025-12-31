function [F,cons] = ReturnFn(aprime_val,a_val,z_val,par)
% Computes return function F(a',a,z) 

wage  = par.wage;
r     = par.r;
sigma = par.sigma;

cons = wage*z_val+(1+r)*a_val-aprime_val;

if isscalar(cons)
    F = -Inf;
    if cons>0
        F = util(cons,sigma);
    end
else
    F = -Inf(size(cons));
    ind_pos = cons>0;
    F(ind_pos) = util(cons(ind_pos),sigma);
end

end %end function "ReturnFn"

