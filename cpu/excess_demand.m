function [ED,KL_implied,V1,apol,cpol,mu] = excess_demand(KLbar,par)
%Compute excess demand as capital demand minus
%capital supply or KL minus KLimplied

%Unpack some parameters:
a_grid = par.a_grid;
z_grid = par.z_grid;
pi_z   = par.pi_z;
Lbar   = par.Lbar;
alpha  = par.alpha;
delta  = par.delta;

%% Prices implied by K/L
[par.r,par.wage] = fun_prices(KLbar,alpha,delta);

%% VFI

[V1,apol,cpol] = sub_vfi(a_grid,z_grid,pi_z,par);


%% Given policy function apol_val and Markov chain for shocks pi_z, compute the stationary distrib mu(a,z)
tol_dist     = par.tol_dist;
maxiter_dist = 10000;
[mu,flag_mu] = MyDist(apol,pi_z,a_grid,tol_dist,maxiter_dist,par);

if flag_mu<0
    warning('Distribution failed!')
end

%% Plot stuff to check
% %assets
% figure(1)
% plot(a_grid,a_grid,'--k',a_grid,apol_val(:,1),'b',a_grid,apol_val(:,nz),'r','linewidth',2)
% legend('45 line','zmin','zmax','location','northwest')
% %consumption
% figure(2)
% plot(a_grid,cpol(:,1),'b',a_grid,cpol(:,nz),'r','linewidth',2)
% %distrib
% figure(3)
% plot(a_grid,sum(mu,2))

%% Compute aggregate K and N implied by household problem

[KL_implied] = aggregates(a_grid,mu,Lbar,par);

ED = KLbar - KL_implied;

fprintf("------------------------------------ \n")
fprintf("errKL   = %f \n", ED)
fprintf("K/L old = %f \n", KLbar)
fprintf("K/L new = %f \n", KL_implied)
fprintf("------------------------------------ \n")

end %END function

%% function nvaluef
% function [val] = nvaluef(kp,np,PARA)
% % negative value function at (kgrid(ik),sgrid(is),hgrid(ih)) 
% % by choosing "kp" as capital next period and 
% % "np" as labor supply this period. 
% % Case 1: np = 0 (no work)
% % Case 2: np = 1 (work)
% 
% delta=PARA.delta; r=PARA.r; w=PARA.w; T=PARA.T; d=PARA.d;
% beta=PARA.beta; chi=PARA.chi; phi=PARA.phi; 
% uform=PARA.uform; gamma=PARA.gamma; EZ=PARA.EZ;
% kgrid=PARA.kgrid; sgrid=PARA.sgrid; hgrid=PARA.hgrid;  
% ik=PARA.ik; is=PARA.is; ih=PARA.ih;  
%   
% % consumption decision given decisions for capital and labor
% cp = (1-delta+r)*kgrid(ik) + sgrid(is)*w*np + T(is,ih) + d - kp;
% % value function (piecewise linear) 
% vfh = utilf(cp,uform,gamma) -chi*np*(1-hgrid(ih))^(phi)+ beta*qinterp1(kgrid,EZ(:,is),kp); 
% val= -vfh; 

%% function utilf
% function [u] = utilf(c,uspec,gamma)
% % [u] = utilf(c,uspec,gamma)
% % 
% % Utility function for "CRRA" or "CARA"
% 
% if strcmp(uspec,'CRRA')==1;
%     % Utility is only defined if c>=0 for gamma<1 and c>0 if gamma>1
%     %u = -inf*ones(size(c)); %start -inf, change only those c positive
%     u = (-1e+100)*ones(size(c)); % this might be better when doing golden section search
%     cpos = find(c>0);       %find positive c
%     if gamma==1; u(cpos)=log(c(cpos));
%     else u(cpos) = (c(cpos).^(1-gamma))./(1-gamma);
%     end
% elseif strcmp(uspec,'CARA')==1;
%     u = -exp(-gamma*c)/gamma;
% end

