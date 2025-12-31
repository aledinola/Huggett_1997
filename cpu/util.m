function U = util(c,sigma)
        
    if sigma==1
        U = log(c);
    else
        U = c.^(1-sigma)/(1-sigma);
    end

end %END function