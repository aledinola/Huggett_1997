function [mu1,flag_mu] = MyDist(apol_val,pi_z,a_grid,tol_dist,maxiter_dist,par)

[na,nz] = size(apol_val);

disp_mu = par.disp_mu;


mu0  = ones(na,nz)/(na*nz);

dist = 10;
iter = 0;

flag_mu = 0;
disp('Start distribution..')

%% Precomputation
[Qmat] = fun_Ptran(apol_val,a_grid);

%% Iterations to find fixed point

while dist>tol_dist && iter<=maxiter_dist
    iter = iter+1;
   
    mu1 = MyDist_onestep(mu0,Qmat,pi_z,na,nz);
    
    %Check
    check = sum(mu1(:));
    
    %Compute error
    dist = max(abs(mu0-mu1),[],"all");
    % Update mu
    mu0 = mu1/check;
    
    if disp_mu==1
        fprintf('sum(mu1) = %f\n', check);
        fprintf('iter = %d, dist = %f\n', iter,dist);
    end
    
end %end while

if dist>tol_dist
    flag_mu = -1;
end


end %end function "MyDist"

