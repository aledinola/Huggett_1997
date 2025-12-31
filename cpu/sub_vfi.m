function [V1,apol,cpol] = sub_vfi(a_grid,z_grid,pi_z,par)
% Value function iteration for a model written in SLP notation

na        = par.na;
nz        = par.nz;
beta      = par.beta;
disp_vfi  = par.disp_vfi;
do_howard = par.do_howard;
n_howard  = par.n_howard;

%% Precompute return matrix F(a',a,z) SLP notation
ReturnMatrix = zeros(na,na,nz);
for z_c = 1:nz
    for a_c = 1:na
        z_val = z_grid(z_c);
        a_val = a_grid(a_c);
        aprime_val = a_grid;
        ReturnMatrix(:,a_c,z_c)=ReturnFn(aprime_val,a_val,z_val,par);
    end
end

%% Value function iteration
V0 = zeros(na,nz);

maxit_v = 1000;
dist = 10;
tol_v = par.tol_v;
iter = 0;


disp("Start VFI..")
while (dist>tol_v && iter<=maxit_v)
    iter = iter+1;
    
    [V1,apol_ind,apol] = sub_vfi_onestep(V0,pi_z,ReturnMatrix,par);
    
    if do_howard==1
        V1 = sub_howard(V1,apol_ind,ReturnMatrix,pi_z,beta,n_howard);
    end
  
    %Compute error
    dist = max(abs(V0(:)-V1(:)));
    % Update V
    V0 = V1;
    
    if disp_vfi==1
        fprintf('iter = %d, dist = %f\n', iter,dist);
    end
    
end %end while

if dist>tol_v
    warning('VFI failed!')
end

%% Compute other policy functions

cpol = zeros(na,nz);
for z_c = 1:nz 
    for a_c=1:na 
        z_val = z_grid(z_c);
        a_val = a_grid(a_c);
        aprime_val = apol(a_c,z_c);
        [~,cpol(a_c,z_c)] = ReturnFn(aprime_val,a_val,z_val,par);
    end
end

end %end function "sub_vfi"

