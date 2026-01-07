% Example: General equilibrium transition path for Huggett (JME 1997)
%
% Author: Alessandro Di Nola
% 
% This script uses the VFI Toolkit to:
%   1. Solve for the stationary general equilibrium of the Huggett (1997) model.
%   2. Construct an initial cross-sectional distribution.
%   3. Compute a general equilibrium transition path for aggregate capital.
%   4. Produce a few diagnostic plots.

clear; clc; close all;

%% Paths and output folders
do_transition = 1; % Set equal to 1 if want to run transition 
toolkit_path  = 'C:\Users\aledi\Documents\GitHub\VFIToolkit-matlab';
addpath(genpath(toolkit_path));

FigDir = 'figures';
if ~exist(FigDir, 'dir')
    mkdir(FigDir);
end

%% Grid sizes

n_a = 2000;  % number of asset grid points
n_z = 2;     % number of idiosyncratic productivity states

%% Parameters

Params.beta  = 0.96;
Params.alpha = 0.36;
Params.delta = 0.10;
Params.sigma = 1.5;  % curvature of utility (CRRA parameter)

% Idiosyncratic productivity shocks
z_grid = [0.8, 1.2]';
pi_z   = [0.5, 0.5;
          0.5, 0.5];

%% Toolkit options

% --- Value function iteration options
vfoptions                     = struct();
vfoptions.lowmemory           = 0;
vfoptions.verbose             = 0;
vfoptions.tolerance           = 1e-9;   % default: 1e-9
vfoptions.maxiter             = 1000;   % default: Inf
vfoptions.howards             = 50;     % default: 150
vfoptions.maxhowards          = 500;    % default: 500
vfoptions.howardsgreedy       = 0;
vfoptions.gridinterplayer     = 1;
vfoptions.ngridinterp         = 20;
vfoptions.divideandconquer    = 0;
% vfoptions.level1n           = 51; % default: 51 if one a variable, 21 if two

% --- Distribution / simulation options
simoptions                    = struct();  % default options for stationary distribution
simoptions.tolerance          = 1e-9;
simoptions.maxit              = 10000;
simoptions.gridinterplayer    = vfoptions.gridinterplayer;
simoptions.ngridinterp        = vfoptions.ngridinterp;

% --- Heterogeneous-agent GE options
heteroagentoptions                          = struct();
heteroagentoptions.verbose                  = 1;     % 1 = print progress
heteroagentoptions.toleranceGEprices        = 1e-6;  % default: 1e-4
heteroagentoptions.toleranceGEcondns        = 1e-6;  % default: 1e-4
heteroagentoptions.fminalgo                 = 0;     % 0=fzero, 1=fminsearch, 8=lsqnonlin
heteroagentoptions.maxiter                  = 1000;

% --- Transition path options
transpathoptions              = struct();
transpathoptions.verbose      = 1;
transpathoptions.weightscheme = 1;

% Subfolder for figures (depending on interpolation setting)
if vfoptions.gridinterplayer == 0
    FigDir1 = fullfile(FigDir, 'discrete');
elseif vfoptions.gridinterplayer == 1
    FigDir1 = fullfile(FigDir, 'interp');
end
if ~exist(FigDir1, 'dir')
    mkdir(FigDir1);
end

%% Grids

% Check that E[z] = 1
z_mean = MarkovChainMoments(z_grid, pi_z);
if abs(z_mean - 1) > 1e-12
    warning('Average of productivity shocks is not equal to one.');
end

% --- Asset holdings
r_ss = 1 / Params.beta - 1;
% Steady-state capital in the absence of aggregate uncertainty
K_ss = ((r_ss + Params.delta) / Params.alpha)^(1 / (Params.alpha - 1));

% Asset grid
a_min  = 0;
a_max  = 10 * K_ss;
a_grid = a_min + (a_max - a_min) * (linspace(0, 1, n_a)'.^3);

% No d-variable in this model
d_grid = 0;
n_d    = 0;

%% Return function, variables to evaluate, and GE conditions

DiscountFactorParamNames = {'beta'};

% First inputs: (aprime, a, z), followed by parameters
ReturnFn = @(aprime, a, z, K, alpha, delta, sigma) ...
    f_ReturnFn(aprime, a, z, K, alpha, delta, sigma);

GEPriceParamNames = {'K'};
Params.K          = 4.31;  % initial guess for aggregate capital in GE

% Functions to evaluate on the cross-sectional distribution
FnsToEvaluate.A = @(aprime, a, z) a;

% General equilibrium condition (capital market clearing)
GeneralEqmEqns.CapitalMarket = @(K, A) K - A;

fprintf('Grid sizes: %d asset points, %d shock states\n', n_a, n_z);

%% Compute the final steady state

fprintf('Calculating price vector corresponding to the stationary general equilibrium...\n');
tic;
[p_eqm_final, ~, GeneralEqmCondn_final] = HeteroAgentStationaryEqm_Case1( ...
    n_d, n_a, n_z, 0, pi_z, d_grid, a_grid, z_grid, ...
    ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, ...
    DiscountFactorParamNames, [], [], [], GEPriceParamNames, ...
    heteroagentoptions, simoptions, vfoptions);
time_ge = toc;

disp(p_eqm_final);  % equilibrium values of the GE prices
% Note: GeneralEqmCondn_final should be essentially zero (GE conditions).

% Update parameters with their equilibrium values
Params.K = p_eqm_final.K;

% Solve value function with final GE capital
[V_final, Policy_final] = ValueFnIter_Case1( ...
    n_d, n_a, n_z, d_grid, a_grid, z_grid, pi_z, ...
    ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);

% Policy functions in levels
PolicyValues_final = PolicyInd2Val_Case1(Policy_final, n_d, n_a, n_z, d_grid, a_grid, vfoptions);
pol_aprime         = reshape(PolicyValues_final, [n_a, n_z]);

FnsToEvaluate.C = @(aprime, a, z, K, alpha, delta) ...
    f_consumption(aprime, a, z, K, alpha, delta);

ValuesOnGrid=EvalFnOnAgentDist_ValuesOnGrid_Case1(Policy_final,FnsToEvaluate,Params,...
    [],n_d, n_a, n_z, d_grid, a_grid, z_grid, simoptions);
pol_cons = ValuesOnGrid.C;

% Stationary distribution at the final steady state
StationaryDist_final = StationaryDist_Case1(Policy_final, n_d, n_a, n_z, pi_z, simoptions);

% Aggregate variables (check)
AggVars_final = EvalFnOnAgentDist_AggVars_Case1( ...
    StationaryDist_final, Policy_final, FnsToEvaluate, Params, [], ...
    n_d, n_a, n_z, d_grid, a_grid, z_grid, simoptions);

% General equilibrium residuals
err_GE           = GeneralEqmEqns.CapitalMarket(Params.K, AggVars_final.A.Mean);
[r, w, Y_agg]    = f_prices(Params.K, Params.alpha, Params.delta);
err_walras       = Y_agg - (AggVars_final.C.Mean + Params.delta * Params.K);

% Open text file for writing results
fid = fopen('results_huggett_1997_toolkit.txt', 'w');

disp('RESULTS FINAL STEADY STATE, VFI TOOLKIT');
fprintf(fid, 'RESULTS FINAL STEADY STATE, VFI TOOLKIT\n');

% Huggett reports steady-state K ≈ 4.3242
fprintf('No. of asset grid points: %d\n', n_a);
fprintf(fid, 'No. of asset grid points: %d\n', n_a);

fprintf('CapitalMarket residual: %f\n', err_GE);
fprintf(fid, 'CapitalMarket residual: %f\n', err_GE);

fprintf('Goods market residual:  %f\n', err_walras);
fprintf(fid, 'Goods market residual:  %f\n', err_walras);

fprintf('Aggregate capital:      %f\n', p_eqm_final.K);
fprintf(fid, 'Aggregate capital:      %f\n', p_eqm_final.K);

fprintf('Capital-to-labor ratio:  %f\n', p_eqm_final.K / z_mean);
fprintf(fid, 'Capital-to-labor ratio:  %f\n', p_eqm_final.K / z_mean);

fprintf('Capital-to-output ratio: %f\n', p_eqm_final.K / Y_agg);
fprintf(fid, 'Capital-to-output ratio: %f\n', p_eqm_final.K / Y_agg);

fprintf('Consumption:             %f\n', AggVars_final.C.Mean);
fprintf(fid, 'Consumption:             %f\n', AggVars_final.C.Mean);

fprintf('Interest rate:           %f\n', r);
fprintf(fid, 'Interest rate:           %f\n', r);

fprintf('Wage:                    %f\n', w);
fprintf(fid, 'Wage:                    %f\n', w);

fprintf('Run time GE:             %f\n', time_ge);
fprintf(fid, 'Run time GE:             %f\n', time_ge);

fclose(fid);

%% Define initial distribution via cutoff on asset grid
% Huggett (1997) initialization:
%   - 20% of agents are exactly at zero assets (a = 0).
%   - The remaining 80% are spread uniformly across the asset grid up to a
%     cutoff c.
%
% We choose c by bisection so that aggregate capital in the initial
% distribution matches the steady-state value p_eqm_final.K.

c_low   = a_grid(1);
c_high  = a_grid(end);
maxiter = 50;
tol     = 1e-6;  % tolerance on aggregate capital

for it = 1:maxiter
    
    % Midpoint cutoff
    c_mid = 0.5 * (c_low + c_high);

    % First grid point strictly above the cutoff
    a_cut_init = find(a_grid > c_mid, 1);

    % Cross-sectional distribution over assets
    pdf_assets    = zeros(n_a, 1);
    pdf_assets(1) = 0.20;  % 20% at zero assets

    if ~isempty(a_cut_init)
        % Spread remaining 80% uniformly across indices 2:a_cut_init
        pdf_assets(2:a_cut_init) = 0.80 / numel(2:a_cut_init);
    end

    % Joint distribution over (a,z)
    StationaryDist_init      = zeros(n_a, n_z);
    StationaryDist_init(:,1) = 0.5 * pdf_assets;
    StationaryDist_init(:,2) = 0.5 * pdf_assets;

    % Normalize (should already sum to 1, but do it for safety)
    StationaryDist_init = StationaryDist_init ./ sum(StationaryDist_init, "all");

    % Aggregate capital implied by current cutoff
    K_mid = sum(a_grid .* StationaryDist_init, "all");

    fprintf('Iter %2d: cutoff c = %.6f, K = %.6f (target %.6f)\n', ...
        it, c_mid, K_mid, p_eqm_final.K);

    % Stopping rule
    if abs(K_mid - p_eqm_final.K) < tol
        fprintf('Tolerance reached. Stopping bisection.\n');
        break;
    end

    % Bisection update
    if K_mid > p_eqm_final.K
        % Too much capital -> reduce cutoff
        c_high = c_mid;
    else
        % Too little capital -> increase cutoff
        c_low = c_mid;
    end

end % end bisection iterations

%% Compute transition path
if do_transition == 1
disp('Start transition computation...');

% V_final: value function in the final steady state
% StationaryDist_init: initial cross-sectional distribution

T = 500;

% Initial guess for the path of aggregate capital
PricePath0.K = 4.325 * ones(T, 1);

% No change in parameters along the path; this is just a constant path
ParamPath.alpha = Params.alpha * ones(T, 1);

TransPathGeneralEqmEqns.CapitalMarket = @(K, A) K - A;

transpathoptions.GEnewprice = 3;
transpathoptions.tolerance  = 1e-4;

% Specify how to update GE prices from GE conditions (for price scheme 3).
% A row has the form: {GEcondnName, priceName, add, factor}
% Here: CapitalMarket > 0 means K is too big, so we want to reduce K.
transpathoptions.GEnewprice3.howtoupdate = ...
    {'CapitalMarket', 'K', 0, 0.1};
% Note: update formula is roughly:
%   new_price = price + factor * add       * GEcondn_value ...
%                       - factor * (1-add) * GEcondn_value

PricePath = TransitionPath_Case1( ...
    PricePath0, ParamPath, T, V_final, StationaryDist_init, ...
    n_d, n_a, n_z, pi_z, d_grid, a_grid, z_grid, ...
    ReturnFn, FnsToEvaluate, TransPathGeneralEqmEqns, Params, ...
    DiscountFactorParamNames, transpathoptions, vfoptions, simoptions);

%% Plots

% Capital over the transition
figure;
plot(1:T, PricePath.K, 'LineWidth', 1.5);
title('Capital over the transition');
xlabel('Time periods');
ylabel('Capital');
print(fullfile(FigDir1, 'Kt_tran.png'), '-dpng');

end %end do transition

% Stationary distribution of assets
figure;
plot(a_grid, sum(StationaryDist_final, 2), 'LineWidth', 1.5);
xlabel('Assets');
ylabel('Density');
title('Distribution of assets in steady state');
print(fullfile(FigDir1, 'stadist.png'), '-dpng');

% Replicate Figure 2 in Huggett (1997): decision rule a'(k,e)
a_cut = find(a_grid > 10, 1);
figure;
plot(a_grid(1:a_cut), pol_aprime(1:a_cut,1), ':',  'LineWidth', 2); hold on;
plot(a_grid(1:a_cut), pol_aprime(1:a_cut,2), '-.', 'LineWidth', 2);
plot(a_grid(1:a_cut), a_grid(1:a_cut), 'k-',       'LineWidth', 2);
legend('a(k,e_1)', 'a(k,e_2)', '45-degree line', ...
       'Location', 'southoutside', 'NumColumns', 3);
xlabel('Capital');
ylabel('Capital next period');
title('Optimal decision rule for capital');
axis tight;
print(fullfile(FigDir1, 'fig2_huggett.png'), '-dpng');

% Policy for consumption 
figure;
plot(a_grid(1:a_cut), pol_cons(1:a_cut,1), ':',  'LineWidth', 2); hold on;
plot(a_grid(1:a_cut), pol_cons(1:a_cut,2), '-.', 'LineWidth', 2);
legend('c(k,e_1)', 'c(k,e_2)', ...
       'Location', 'southoutside', 'NumColumns', 3);
xlabel('Capital');
ylabel('Consumption');
title('Optimal decision rule for consumption');
axis tight;
print(fullfile(FigDir1, 'fig_consumption.png'), '-dpng');

% Initial vs stationary distributions
figure;
plot(a_grid, sum(StationaryDist_init, 2), '-.', 'LineWidth', 2); hold on;
plot(a_grid, sum(StationaryDist_final, 2), 'k-', 'LineWidth', 2);
legend('Initial distribution', 'Stationary distribution', ...
       'Location', 'southoutside', 'NumColumns', 2);
xlabel('Capital');
ylabel('Density');
title('Initial vs stationary distributions');
