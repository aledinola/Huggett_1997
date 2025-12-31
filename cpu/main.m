% Huggett 1997 model, inelastic labor supply (i.e. always equal to one)
% Transition after arbitrary initial distribution, see JME paper for
% details
clear,clc,close all

format long g
FigDir = fullfile('figures'); 
MatDir = fullfile('mat_files');

%diary myDiaryFile

%{
Huggett (JME 1997) model:
u(c) = c^(1-sigma)/(1-sigma), log(c) if sigma=1
c+a'=(1+r*)*a+w*z
z follows a Marko chain with two states, z=[0.8,1.2]
and pi_z = [0.5,0.5;0.5,0.5] so that E(z)=1
Prices: r = alpha*(K/L)^(alpha-1), w = (1-alpha)(K/L)^alpha
but note that L=E(z)=1
Aggregate resource constraint: Y=C+delta*K
%}

%% Set flags and computational parameters
do_plots    = 1;     % flag 0/1 to draw plots
do_save     = 1;     % flag save plots as png
verbose     = 1;     % 0=no display,1=limited info,2=all
disp_vfi    = 0;     % flag 0/1 to display iter in VFI steady-state
disp_mu     = 0;     % flag 0/1 to display iter in distrib. steady-state
do_golden   = 0;     % flag 0/1 for golden algo
do_howard   = 1;     % flag 0/1 for Howard
n_howard    = 50;    % No. iteration in Howard algorithm
do_ss_final = 1;     % flag 0/1 for final s.s.
do_trans    = 1;     % flag 0/1 for transition
tol_v       = 1e-9;  % Tolerence for VFI 
tol_dist    = 1e-9;  % Tolerance for distribution

%% Set economic parameters
beta  = 0.96;
alpha = 0.36;
sigma = 1.5; % CRRA
delta = 0.1;

KL_0_final = 4.31; % Initial condition for K/L for final steady state

%% Income process
nz = 2;
z_grid = [0.8,1.2]';
pi_z = [0.5,0.5;0.5,0.5];
Lbar = 1;

%% Set grid for assets

r_ss=1/beta-1;
K_ss_norisk=((r_ss+delta)/alpha)^(1/(alpha-1)); %The steady state capital in the absence of aggregate uncertainty.

na     = 2000;    % No. grid points for assets
amin   = 0.0;     % Minimum assets
amax   = 10*K_ss_norisk; % Maximum assets
curv   = 3;       % Spacing of grid (=1 is linspace)

a_grid = amin + (amax-amin)*(linspace(0,1,na).^curv)';

%% Pack parameters into a structure

par = pack_into_struct(Lbar,beta,sigma,alpha,delta,nz,na,a_grid,z_grid,pi_z,...
    verbose,do_golden,disp_vfi,disp_mu,do_howard,n_howard,tol_v,tol_dist);

%% Define initial distribution
% Huggett: The initial distribution puts 20% of the agents exactly at
%zero asset holdings and equal numbers of agents at all capital levels between
%0 and 10.8104. He does not report how he constructed the asset grid though
a_cut_init = find(a_grid>10.81, 1 );
pdf_assets = zeros(na,1);
pdf_assets(1) = 0.2;
pdf_assets(2:a_cut_init) = 0.8/numel(2:a_cut_init);

AgentDist_init = zeros(na,nz);
AgentDist_init(:,1) = 0.5*pdf_assets;
AgentDist_init(:,2) = 0.5*pdf_assets; 
AgentDist_init = AgentDist_init/sum(AgentDist_init,"all");


%% Compute the final steady-state

if do_ss_final==1

par.KL_0  = KL_0_final;  % improved guess for equilibrium K/L
tic
disp("Final steady-state..")
[KLbar_final,V_final,Policy_final,Dist_final,err_GE] = compute_steady_state(par);
time_ge=toc;

%Check market clearing condition
% Y = C+delta*K
apol = Policy_final.apol;
cpol = Policy_final.cpol;
mu   = Dist_final;

Kss = KLbar_final*Lbar;
Css = sum(cpol.*mu,'all');
Yss = Kss^par.alpha*Lbar^(1-par.alpha);
[rss,wss] = fun_prices(KLbar_final,alpha,delta);

fprintf(' \n')
err_walras = Yss - (Css+delta*Kss);
        
% Open text file for writing
fid = fopen('results_huggett_1997_cpu.txt','w');

fprintf('FINAL STEADY-STATE, OWN CPU CODE: \n');
fprintf(fid,'FINAL STEADY-STATE, OWN CPU CODE: \n');

fprintf('No. grid points assets: %d \n',na);
fprintf(fid,'No. grid points assets: %d \n',na);

fprintf('CapitalMarket residual: %f \n',err_GE);
fprintf(fid,'CapitalMarket residual: %f \n',err_GE);

fprintf('Goods market residual:  %f \n',err_walras);
fprintf(fid,'Goods market residual:  %f \n',err_walras);

fprintf('Capital stock:           %f \n',Kss);
fprintf(fid,'Capital stock:           %f \n',Kss);

fprintf('Capital-to-labor ratio:  %f \n',KLbar_final);
fprintf(fid,'Capital-to-labor ratio:  %f \n',KLbar_final);

fprintf('Capital-to-output ratio: %f \n',Kss/Yss);
fprintf(fid,'Capital-to-output ratio: %f \n',Kss/Yss);

fprintf('Consumption:             %f \n',Css);
fprintf(fid,'Consumption:             %f \n',Css);

fprintf('Interest rate:           %f \n',rss);
fprintf(fid,'Interest rate:           %f \n',rss);

fprintf('Wage:                    %f \n',wss);
fprintf(fid,'Wage:                    %f \n',wss);

fprintf('Run time General Equil: %f \n',time_ge);
fprintf(fid,'Run time General Equil: %f \n',time_ge);

% Close the file
fclose(fid);

save(fullfile(MatDir,'ss_final.mat'),'KLbar_final','V_final','Policy_final','Dist_final');

% Plot stationary distribution
figure
plot(a_grid,sum(Dist_final,2))
xlabel('Assets')
title('Distribution of assets in steady state')
print(fullfile(FigDir,'stadist.png'),'-dpng')

% Replicate Figure 2 of Huggett JME paper
a_cut = find(a_grid>10, 1 );
figure
plot(a_grid(1:a_cut),apol(1:a_cut,1),':','linewidth',2)
hold on
plot(a_grid(1:a_cut),apol(1:a_cut,2),'-.','linewidth',2)
hold on 
plot(a_grid(1:a_cut),a_grid(1:a_cut),'k-','linewidth',2)
legend('a(k,e_1)','a(k,e_2)','45 degree line','Location','southoutside','NumColumns', 3)
xlabel('CAPITAL')
ylabel('CAPITAL NEXT PERIOD')
title('OPTIMAL DECISION RULE')
axis tight
print(fullfile(FigDir,'fig2_huggett.png'),'-dpng')

end

%% Compute the transition path
% For the transition path we will need the initial agents distribution
% For the transition path we will need the final value function

% Compute the transition path
% For this we need the following extra objects: PricePathOld, T, V_final, StationaryDist_init
% (already calculated V_final & StationaryDist_init above)

% Number of time periods to allow for the transition (if you set T too low
% it will cause problems, too high just means run-time will be longer).
if do_trans==1
T = 100;

% We want to look at a one off unanticipated change of alpha. ParamPath & PathParamNames are thus given by
ParamPath.alpha=alpha*ones(T,1); % For each parameter that changes value, ParamPath is matrix of size T-by-1
% (the way ParamPath is set is designed to allow for a series of changes in the parameters)

% We need to give an initial guess for the price path on interest rates
% (or,alternatively, the K/L ratio)
PricePath0.KL = [linspace(KLbar_final, KLbar_final, floor(T/2))'; KLbar_final*ones(T-floor(T/2),1)]; % For each price, PricePath0 is matrix of size T-by-1

fprintf(' \n')
fprintf(' \n')
disp('START TRANSITION')
fprintf(' \n')

[PricePathNew,r_path,wage_path,err_tran_vec] = compute_transition(PricePath0, ParamPath, T, V_final, AgentDist_init,par);

err_tran_vec(err_tran_vec==0) = [];

save(fullfile(MatDir,'temp_transition.mat'))

figure
plot(err_tran_vec,'linewidth',2)
xlabel('Iteration')
ylabel('Distance b/w old path and new path')
title('Error at each iteration of transition loop')
if do_save==1; print(fullfile(FigDir,'err_path'),'-dpng'); end

figure(1)
plot(1:T,PricePathNew,'linewidth',2)
xlabel("Periods")
ylabel("K_t/L_t")
title("The time series K_t/L_t")
if do_save==1; print(fullfile(FigDir,'KL_path'),'-dpng'); end

figure(2)
plot(1:T,r_path,'linewidth',2)
xlabel("Periods")
ylabel("r_t")
title("Interest rate along the transition")
if do_save==1; print(fullfile(FigDir,'r_path'),'-dpng'); end

figure(3)
plot(1:T,wage_path,'linewidth',2)
xlabel("Periods")
ylabel("w_t")
title("Wage along the transition")
if do_save==1; print(fullfile(FigDir,'wage_path'),'-dpng'); end

end

diary off

