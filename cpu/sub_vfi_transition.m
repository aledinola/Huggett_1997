function [PolicyPath] = sub_vfi_transition(V_final,PricePathOld,ParamPath,T,par)

% Unpack stuff
alpha   = par.alpha;
delta   = par.delta;
na      = par.na;
nz      = par.nz;
pi_z    = par.pi_z;
a_grid  = par.a_grid;
z_grid  = par.z_grid;
verbose = par.verbose;

% Initialize output
% We need policies only up to T-1 since we never update
% price(T)=price_final by construction
PolicyPath = zeros(na,nz,T-1);

Vnext = V_final; % Terminal condition for backward iterations

% Precompute arrays for calculation of ReturnMatrix
Aprime = reshape(a_grid, na, 1, 1);   % na × 1 × 1
A      = reshape(a_grid, 1, na, 1);   % 1 × na × 1
Z      = reshape(z_grid, 1, 1, nz);   % 1 × 1 × nz

disp("Start backward iteration from T-1 to 1")
for t = T-1:-1:1
    if verbose==2
        fprintf('t = %d \n',t);
    end
    % Parameters that change in transition
    %par.tau_k = ParamPath.tau_k(t);
    
    % Prices implied by K/L
    KLbar_t = PricePathOld(t);
    [par.r,par.wage] = fun_prices(KLbar_t,alpha,delta);
    % Static return changes during each t of the transition since r,w
    % change, so it has to be recomputed for each t
    
    ReturnMatrix = ReturnFn(Aprime, A, Z, par);

    [V,~,apol_val] = sub_vfi_onestep(Vnext,pi_z,ReturnMatrix,par);

    Vnext = V;
    PolicyPath(:,:,t) = apol_val;

end


end %end function "sub_vfi_transition"

