function compare_asymptotic_vs_numerical
%COMPARE_ASYMPTOTIC_VS_NUMERICAL_  Convergence test for a matched solution.
%
%   PURPOSE
%   -------
%   Validate a **matched asymptotic approximation** y_uniform(x; ε) obtained
%   via boundary‐layer theory by showing that its L₁-distance to the direct
%   numerical solution of the BVP decays algebraically as ε → 0.
%
%   MODEL PROBLEM
%   -------------
%       ε y'''(x) − y'(x) + x y(x) = 0,  x ∈ [0, 1],
%       y(0) = 1, y'(0) = 1, y(1) = 1.
%
%   The script:
%     1) Evaluates the matched asymptotic formula y_uniform(x; ε).
%     2) Solves the full BVP with `bvp4c` for several ε values.
%     3) Computes   L₁(ε) = ∫₀¹ |y_uniform − y_num| dx.
%     4) Plots both curves for visual comparison.
%     5) Fits   log L₁ ≈ p log ε + C   to expose the observed convergence
%        rate p.
%
%   EXPECTED RESULT
%   ---------------
%   On a log–log plot, L₁(ε) should decrease with slope p < 0,
%   confirming that the asymptotic matching captures the correct interior
%   and boundary‐layer behaviour as ε becomes small.
%
%   AUTHOR  : Sagi Senderovich
%   VERSION : 22-Jun-2025
% -------------------------------------------------------------------------

clc; close all;

%% --------------------------- user parameters ---------------------------
EPSILON_VEC = [1e-3, 1e-4, 1e-6, 1e-8, 1e-10];
X_GRID      = linspace(0, 1, 1000);        % common evaluation grid

% output folder
fig_dir = fullfile(pwd,'fig');
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end
% ------------------------------------------------------------------------

n_cases   = numel(EPSILON_VEC);
color_map = lines(n_cases);
L1_error  = zeros(n_cases, 1);

%% ----------------------- figure 1: solutions ----------------------------
fig1 = figure('Name', 'Asymptotic vs. Numerical');
ylim([1 1.7])
hold on;

for k = 1:n_cases
    eps = EPSILON_VEC(k);

    % asymptotic (asymptotic) approximation  y_uniform(x; ε)
    uniform_fun = @(x) ...
        ((sqrt(eps) + 1) .* exp(x.^2 / 2) ...
        - sqrt(eps) .* exp(-x ./ sqrt(eps)) ...
        + (1 - sqrt(exp(1)) .* (sqrt(eps) + 1)) ...
          .* exp((x - 1) ./ sqrt(eps)));

    y_uniform = uniform_fun(X_GRID);

    % numerical BVP  (first-order system)
    %   Y(1) = y ,  Y(2) = y' ,  Y(3) = y''
    ode = @(x, Y)[ Y(2) ; ...
                   Y(3) ; ...
                   (Y(2) - x .* Y(1)) / eps ];

    bc  = @(Ya, Yb)[ Ya(1) - 1 ;   % y(0)  = 1
                     Ya(2) - 1 ;   % y'(0) = 1
                     Yb(1) - 1 ];  % y(1)  = 1

    init_guess = @(x)[ 1 ; 1 ; 0 ];                 % coarse initial guess
    sol        = bvp4c(ode, bc, bvpinit(X_GRID, init_guess));
    y_num      = deval(sol, X_GRID);

    % L₁-error
    L1_error(k) = trapz(X_GRID, abs(y_uniform - y_num(1 , :)));

    % plotting
    plot(X_GRID, y_uniform, '-',  'LineWidth',1.4, ...
        'Color',color_map(k,:), ...
        'DisplayName',sprintf('Asymptotic \\epsilon = %.0e', eps));

    plot(X_GRID, y_num(1 , :), '--', 'LineWidth',1.4, ...
        'Color',color_map(k,:), ...
        'DisplayName',sprintf('Numerical \\epsilon = %.0e', eps));
end

xlabel('x'); ylabel('y(x)');
title({'Boundary-Layer Solutions'; '(asymptotic vs. numerical BVP)'});
legend('Location','best'); grid off; box off;

saveas(fig1, fullfile(fig_dir,'solutions_vs_numerical.png'));


%% ------------------- figure 2: L1 error convergence ---------------------
fig2 = figure('Name','L1 error (log‒log)');
loglog(EPSILON_VEC, L1_error, 'o-', 'LineWidth',1.5);
xlim([0 0.001]);
grid off; box off; xlabel('\epsilon'); ylabel('L_1 error');
title('Convergence of Asymptotic Solution');

% fit slope  log(L1) = p log ε + C
coeff = polyfit(log(EPSILON_VEC(:)), log(L1_error(:)), 1);
rate  = coeff(1);  intercept = coeff(2);

text('Units','normalized', 'Position',[0.55 0.25], ...
     'String',sprintf('L_1 \\approx %.2f \\epsilon^{%.2f}', ...
                      exp(intercept), rate), ...
     'BackgroundColor','w', 'EdgeColor','k');

fprintf('Observed fit:  log(L1) = %.2f·log(ε) + %.2f\n', rate, intercept);
saveas(fig2, fullfile(fig_dir,'L1_error_loglog.png'));
end
