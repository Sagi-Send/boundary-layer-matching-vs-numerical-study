function compare_uniform_vs_numerical
%COMPARE_UNIFORM_VS_NUMERICAL  Matched-asymptotic “uniform” solution vs.
%                              direct numerical BVP solution.
%
%   PROBLEM
%   -------
%   We study the singularly-perturbed boundary-layer ODE
%
%       ε y''(x) + (1/x) y'(x) + y(x) = 0,  x ∈ [ε, 1],
%       y(ε) = 0, y(1) = exp(−1/2).
%
%   A composite uniform approximation y_uniform(x) is obtained via
%   matched-asymptotic expansions (outer + inner layer, then matched).
%   This script compares that analytic approximation with the numerical
%   solution produced by `bvp4c`.
%
%   WORKFLOW
%   --------
%   • For every value in EPSILON_VEC
%       1. Evaluate the composite uniform approximation.
%       2. Solve the BVP numerically with first-order reduction + `bvp4c`.
%       3. Compute the L₁-error on [ε, 1].
%       4. Plot both curves on the same axes.
%   • Fit the convergence rate log(L₁) ≈ C ε^p and show it on a log–log
%     plot.
%   • Save both figures (PNG) to the folder `./fig/`.
%
%   FILES WRITTEN
%   -------------
%   fig/solutions_vs_numerical.png   – uniform vs. numerical curves
%   fig/L1_error_loglog.png          – L₁-error convergence
%
%   AUTHOR  : Sagi Senderovich
%   VERSION : 22-Jun-2025
% -------------------------------------------------------------------------

clc;  close all;

%% --------------------------- user parameters --------------------------
epsilon_vec = [1e-1, 1e-2, 1e-3, 1e-4];
x_common    = linspace(min(epsilon_vec), 1, 1000);
save_folder = fullfile(pwd,'fig');           % ./fig   (auto-created)
if ~exist(save_folder,"dir"), mkdir(save_folder); end
% -----------------------------------------------------------------------

n_cases   = numel(epsilon_vec);
color_map = lines(n_cases);
L1_error  = zeros(n_cases,1);

% ---------------------------- Figure 1 ---------------------------------
fig1 = figure('Name','Legacy vs. numerical solutions');
hold on;

for k = 1:n_cases
    epsilon = epsilon_vec(k);

    % ---------- legacy “uniform” solution (kept verbatim) --------------
    legacy_fun = @(x) exp(-x.^2/2) - exp(-(x - epsilon)./epsilon.^2);
    x_legacy   = linspace(epsilon,1,1000);
    y_legacy   = legacy_fun(x_legacy);

    % ---------------- numerical BVP solution ---------------------------
    ode = @(x,Y)[ Y(2) ;
                 -(1./(x*epsilon)).*Y(2) - (1/epsilon).*Y(1) ];
    bc  = @(Ya,Yb)[ Ya(1) ;
                    Yb(1)-exp(-1/2) ];
    guess = @(x)[ ((x-epsilon)/(1-epsilon))*exp(-1/2) ;
                  exp(-1/2)/(1-epsilon) ];
    sol   = bvp4c(ode,bc,bvpinit(linspace(epsilon,1,400),guess));

    % safe evaluation on [epsilon,1]
    valid_idx      = x_common >= epsilon;
    y_num_full     = NaN(1,numel(x_common));
    y_temp         = deval(sol,x_common(valid_idx));
    y_num_full(valid_idx) = y_temp(1,:);
    y_num = y_num_full;

    % L1 error
    L1_error(k) = trapz(x_common(valid_idx), ...
                        abs(legacy_fun(x_common(valid_idx))-y_num(valid_idx)));

    % plotting
    plot(x_legacy,y_legacy,'-','LineWidth',1.4,'Color',color_map(k,:), ...
        'DisplayName',sprintf('Legacy  \\epsilon = %.0e',epsilon));
    plot(x_common(valid_idx),y_num(valid_idx),'--','LineWidth',1.4, ...
        'Color',color_map(k,:), ...
        'DisplayName',sprintf('Numerical \\epsilon = %.0e',epsilon));
end

xlabel('x'); ylabel('y(x)');
title({'Transition-layer solutions';'(legacy analytic vs. numerical BVP)'});
legend('Location','best'); grid on; box on;

% save figure 1
saveas(fig1, fullfile(save_folder,'solutions_vs_numerical.png'));

% ---------------------------- Figure 2 ---------------------------------
fig2 = figure('Name','L1 error (log-log)');
loglog(epsilon_vec,L1_error,'o-','LineWidth',1.5); grid on;
xlabel('\epsilon'); ylabel('L_1 error');
title('Convergence of legacy solution towards numerical BVP');

% fit slope
coeff = polyfit(log(epsilon_vec(:)), log(L1_error(:)),1);
rate  = coeff(1); intercept = coeff(2);
text(epsilon_vec(2),L1_error(2)*1.4, ...
     sprintf('L_1 \\approx %.2f \\epsilon^{%.2f}',exp(intercept),rate), ...
     'BackgroundColor','w','EdgeColor','k');

fprintf('Observed fit:  log(L1) = %.2f·log(epsilon) + %.2f\n',rate,intercept);

% save figure 2
saveas(fig2, fullfile(save_folder,'L1_error_loglog.png'));
end
