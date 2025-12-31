function [PricePathNew,r_path,wage_path,err_tran_vec] = compute_transition(PricePath0,...
    ParamPath, T, V_final, AgentDist_init,par)

% DESCRIPTION
% Computes the transition to MIT shock. The parameters that change between
% steady state are stored in ParamPath.
% INPUTS
% PricePath0:     Struct with initial guess for price sequence (e.g.
%                 interest rate or K/L ratio)
% ParamPath:      Struct with parameters that change along the transition
%                 (e.g. tax rates)
% T:              Length of transition
% V_final:        Value function final steady-state
% AgentDist_init: Distribution initial steady-state



%Unpack structure
%beta  = par.beta;
%tau_k = ParamPath.tau_k(1); %parameter that changed at t=1 (and stays there)
delta = par.delta;
alpha = par.alpha;
%sigma = par.sigma;

na = par.na;
nz = par.nz;
a_grid = par.a_grid;
%z_grid = par.z_grid;
pi_z   = par.pi_z;
Lbar   = par.Lbar;
verbose = par.verbose;

%amin = min(a_grid);
%amax = max(a_grid);

%Transition parameters
oldpathweight = 0.9*ones(1,T);

%speed of updating the interest rate
%xi = 1*(exp(-0.05*(1:T)) - exp(-0.05*T));
%oldpathweight = 1-xi;

tol = 1e-3;
pathcounter = 0;
maxiterations = 1000;
err_tran = 10;
%PolicyPath = zeros(na,nz,T);

err_tran_vec = zeros(maxiterations,1);
PricePathOld = PricePath0.KL;

tic

while err_tran>tol && pathcounter<maxiterations
    
    pathcounter = pathcounter+1;
    
    %% Backward iteration on Value Function
    %First, go from T-1 to 1 calculating the Value function and Optimal
    %policy function at each step. Since we won't need to keep the value
    %functions for anything later we just store the next period one in
    %Vnext, and the current period one to be calculated in V
    
    [PolicyPath] = sub_vfi_transition(V_final,PricePathOld,ParamPath,T,par);
    % PolicyPath has dim: (na,nz,T-1)
     
    %% Forward iteration on the Distribution
    % Now we have the full PolicyPath, we go forward in time from 1
    % to T using the policies to update the agents distribution generating a
    % new price path
    fprintf(' \n')
    disp("Start forward iteration on the Distribution")
    AgentDist       = AgentDist_init;
    PricePathNew    = zeros(T,1);
    PricePathNew(T) = PricePath0.KL(T); % equal to KL_final
    
    for t = 1:T-1
        if verbose==2
            fprintf('t = %d \n',t);
        end
        
        % Get the current optimal policy
        Policy = PolicyPath(:,:,t);
        Ptran  = fun_Ptran(Policy,a_grid);
        
        % Compute next-period distribution
        AgentDistnext = MyDist_onestep(AgentDist,Ptran,pi_z,na,nz);
        
        KL_implied      = aggregates(a_grid,AgentDist,Lbar,par);
        PricePathNew(t) = KL_implied;
        
        AgentDist = AgentDistnext;
    
    end
    
    %% See how far apart the price paths are and then update
    err_tran=max(abs(reshape(PricePathNew(1:T-1,:)-PricePathOld(1:T-1,:),[numel(PricePathOld(1:T-1,:)),1])));
    err_tran_vec(pathcounter) = err_tran;
    % Notice that the distance is always calculated ignoring the time t=1 &
    % t=T periods, as these needn't ever converges
    
     fprintf("------------------------------------ \n")
     fprintf("iter_tran = %d \n", pathcounter)
     fprintf("err_tran  = %f \n", err_tran)
     fprintf("------------------------------------ \n")
    
    if verbose>=2
        disp('Old, New')
        disp([PricePathOld,PricePathNew])
    end
    
    PricePathOld(1:T-1)=oldpathweight(t)*PricePathOld(1:T-1)+(1-oldpathweight(t))*PricePathNew(1:T-1);
    
    if mod(pathcounter,1)==0
        figure(1);
        %clf;
        %subplot(2,2,1);
        plot(PricePathOld,'linewidth',2);
        hold on;
        title('KL ratio');
%         plot(r_t);
%         subplot(2,2,2);
%         plot(dS);
%         title('dS');
%         subplot(2,2,3);
%         plot(zeros(1,21),'r--');
%         hold on;
%         plot(dS(1:20));
%         title('dS(1:20)');
%         subplot(2,2,4);
%         plot(dS(N-20:N));
%         hold on;
%         title('dS(N-20:N)');
        pause(0.1);
    end
    
end %end outer loop over prices

disp('Time to compute transition:')
toc

%From path for K(t)/L(t), compute r(t) and w(t)
[r_path,wage_path] = fun_prices(PricePathNew,alpha,delta);

end %end function "compute_transition"

