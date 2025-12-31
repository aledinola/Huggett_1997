function mu1 = MyDist_onestep(mu0,Ptran,pi_z,na,nz)

% INPUTS:
% mu0: Current-period distribution, dim: (na,nz)
% Ptran: Big transition matrix, dim: [na*nz,na*nz]. Rows sum to 1
% OUTPUTS:
% mu1: Next-period distribution, dim: (na,nz)

check_in = sum(mu0(:));
if abs(check_in-1)>1e-6
    error("MU0 does not sum to one")
end

mu0_vec = mu0(:);   % column vector

mu1_vec = Ptran'*mu0_vec; % Column vector of size [na*nz,1]

mu1 = reshape(mu1_vec,na,nz); % Matrix of size [na,nz]

mu1 = mu1*pi_z;

%Check
check_out = sum(mu1(:));
if abs(check_out-1)>1e-6
    error("MU1 does not sum to one")
end

end %end function "MyDist_onestep"

