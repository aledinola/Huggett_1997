function [KLbar,V,Policy,Dist,residual] = compute_steady_state(par)

%{
Inputs:
par: structure with grids and model parameters

Outputs:
K_L: capital-labor ratio in equilibrium
V: value function, dim: (na,nz)
Policy: Structure with policy functions apol and cpol, each with dim:
(na,nz)
Dist: stationary distribution, dim: (na,nz)
%} 

% Solve general equilibrium conditions:
%x0 is the initial interval that should bracket a root of the market
%clearing condition
x0 = par.KL_0;

options = optimset('Display','iter','TolX',1e-6);
[KLbar,residual] = fzero(@(KL) excess_demand(KL,par),x0,options);
fprintf("K/L      = %f \n", KLbar)
fprintf("residual = %f \n", residual)

%Now that the equilibrium prices (actually, the equilibrium K/L) have been find, 
%compute V and Policy
[~,~,V,apol,cpol,Dist] = excess_demand(KLbar,par);

Policy.apol = apol;
Policy.cpol = cpol;

end %end function "compute_steady_state"

