% Global Senstivity Analysis using-
% SAFE (Sensitivity Analysis For Everybody) Toolbox- http://bristol.ac.uk/cabot/resources/safe-toolbox/
% Pianosi, F., Sarrazin, F., Wagener, T. (2015), 
% A Matlab toolbox for Global Sensitivity Analysis, Environmental Modelling & Software, 70, 80-85. 


n_params =10; % [thetae thetai thetaa thetaa2 
              % taue taui taua taua2 tauadap beta]

params_lower= [0.1 15 -10 0.1 0 0 10 10 450 0];

params_upper= [20 35 -0.1 10 20 10 30 30 550 10];

DistrFun = 'unif'; % TODO change distribution
DistrPar = cell(n_params,1);
for i =1:n_params; DistrPar {i} = [ params_lower(i) params_upper(i) ]; end

params_label = {'\theta_{E}', '\theta_{I}', '\theta_{A_{d}}', '\theta_{A_{h}}', '\tau_{E}', '\tau_{I}', '\tau_{A_{d}}', '\tau_{A_{h}}', '\tau_{a}', '\beta'};
fun = 'EE_findperUP';

r = 100;
SampStrategy = 'lhs';
design_type = 'radial';
X = OAT_sampling (r,n_params, DistrFun, DistrPar, SampStrategy, design_type);
Y= model_evaluation (fun, X);
[mi, sigma] = EET_indices(r, params_lower, params_upper, X, Y, design_type);


Nboot = 100;
[mi, sigma, EE ,mi_sd, sigma_sd, mi_lb, sigma_lb, mi_ub, sigma_ub] = ...
EET_indices(r,params_lower, params_upper, X, Y, design_type, Nboot);

%%%%
%plotting code
colors = distinguishable_colors(n_params);

% plotting mi and std
for i= 1:n_params
    h = fill([mi_lb(i),mi_lb(i),mi_ub(i),mi_ub(i)], ...
        [sigma_lb(i),sigma_ub(i),sigma_ub(i),sigma_lb(i)], colors(i,:), 'EdgeColor', colors(i,:));
    set(h,'facealpha',0.1)
    hold on
end


for i= 1:n_params
    k(i) = plot(mi(i),sigma(i), 'ok', "MarkerFaceColor", colors(i,:), 'MarkerSize', 14);
    hold on
end
hold off

legend(k, params_label, 'Location', 'southeast', 'FontSize', 15);
xlabel('Mean of EEs','FontSize', 18);
ylabel('Standard deviation of EEs','FontSize',18);

ax = gca;
ax.FontSize = 16; 




