function [V1,apol_ind,apol] = sub_vfi_onestep(V0,pi_z,ReturnMatrix,par)

% DESCRIPTION:
% Performs one step of the Bellman operator. Called to find fixed point for
% steady-state or to iterate backwaReturnMatrixrds for transition
% INPUTS:
%   "V0"    Current-period value function, dim: (na,nz)
%   "pi_z"  Markov chain for exogenous state, dim: (nz,nz)
%   "ReturnMatrix" Payoff matrix, dim: (na',na,nz)
%   "par"   Struct with model parameters and flags
% OUTPUTS:
%   "V1"    Updated value function, dim: (na,nz)

[na,nz] = size(V0);

beta   = par.beta;
a_grid = par.a_grid;

%V1       = zeros(na,nz);
%apol     = zeros(na,nz);
%apol_ind = ones(na,nz);

%Compute expected value EV(a',z)=sum_z_prime V(a',z')*Ptran(z,z')'
EV = V0*pi_z';

% for z_c = 1:nz % Exogenous shock
%     for a_c = 1:na % Endogenous state
%         RHS_vec             = ReturnMatrix(:,a_c,z_c)+beta*EV(:,z_c);
%         [V1(a_c,z_c),ap_ind] = max(RHS_vec);
%         apol_ind(a_c,z_c)   = ap_ind;
%         apol(a_c,z_c)       = a_grid(ap_ind);
%     end
% end    

% RetMat is (a',a,z), EV is (a',z)-->EV(a',1,z)
RHS_mat     = ReturnMatrix+beta*reshape(EV,[na,1,nz]);
[V1,ap_ind] = max(RHS_mat,[],1);
apol_ind    = reshape(ap_ind,[na,nz]);
apol        = a_grid(apol_ind);
V1          = reshape(V1,[na,nz]);

end %end function "sub_vfi_onestep"

