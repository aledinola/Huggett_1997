function [Ptran] = fun_Ptran(Policy,a_grid)

% DESCRIPTION:
% Creates the big transition matrix Ptran which gives the transition (a,z)
% to (a',z). This is for first step of Tan improvement.
% INPUTS:
%   "Policy"  Current-period policy function, dim: (na,nz)
%   "a_grid"  Grid for endogenous state
% OUTPUTS:
%   "Ptran"   Transition matrix with rows summing up to one, dim:(na*nz,na*nz)
% AUXILIARY
%   "find_loc_vec"  Used to find interp indeces and weights (based on histc)

[na,nz] = size(Policy);
NA = (1:na)';

% Find interpolating indexes and weights
% aInd is the location of the left point
[aInd,omega] = find_loc_vec2(a_grid,Policy);

% Build matrix G
G = cell(nz,1);

% For each z, build an na*na matrix transition matrix, called G{z}, 
% that represents the policy function a'=g(a,z).
for iz=1:nz
    % each iz-cell is [na,na] matrix
    G{iz} = sparse(NA,aInd(:,iz),omega(:,iz),na,na)+...
       sparse(NA,aInd(:,iz)+1,1-omega(:,iz),na,na); 
end %close iz
    
% Build transition matrix Q, dim: [na*nz,na*nz]. Rows sum to 1
% In the first step of Tan, the shock z does not change: hence the big
% transition matrix is block diagonal
Ptran = blkdiag(G{:});

end % end function "fun_Ptran"

