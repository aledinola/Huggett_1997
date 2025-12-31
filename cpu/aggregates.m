function [KL_implied] = aggregates(a_grid,mu,Lbar,par)

na = par.na;
nz = par.nz;

K_implied = 0;
for z_c = 1:nz 
    for a_c=1:na 
        K_implied = K_implied +a_grid(a_c)*mu(a_c,z_c);
    end
end

KL_implied = K_implied/Lbar;

end % end function "aggregates"

