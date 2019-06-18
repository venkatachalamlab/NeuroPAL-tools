function stats = computeFunctStructStats(ca_table, varargin)
%COMPUTEFUNCTSTRUCTSTATS Compute correlations between the functional &
%   structural connectomes for a set of data inclusion criteria.
%
% Input:
%   ca_table   = the connectome to use
%   scale_func = a scaling function for weights (e.g., log)
%                default: weight = weight
%
% Output:
%   stats = the statistical correlations between the connectomes

% Are we scaling the weights?
scale = @(x) (x);
if ~isempty(varargin)
    scale = varargin{1}; 
end

% Initialize the weights.
ca_weight = scale(1 - ca_table.calcium_dist);
chemical_weight = scale(ca_table.chemical_weight);
electrical_weight = scale(ca_table.electrical_weight);
total_weight = scale(ca_table.total_weight);

% Compute the chemical weight statistics.
data = [];
i = ca_table.chemical_weight > 0;
[data.pearson.r, data.pearson.p] = ...
    corr(chemical_weight(i), ca_weight(i), 'Type', 'Pearson');
[data.spearman.r, data.spearman.p] = ...
    corr(chemical_weight(i), ca_weight(i), 'Type', 'Spearman');
[data.kendall.r, data.kendall.p] = ...
    corr(chemical_weight(i), ca_weight(i), 'Type', 'Kendall');
[pfit, pgof] = fit(chemical_weight(i), ca_weight(i), 'poly1');
data.poly1.fit = pfit;
data.poly1.gof = pgof;
[pfit, pgof] = fit(chemical_weight(i), ca_weight(i), 'poly2');
data.poly2.fit = pfit;
data.poly2.gof = pgof;
stats.chemical = data;

% Compute the electrical weight statistics.
data = [];
i = ca_table.electrical_weight > 0;
[data.pearson.r, data.pearson.p] = ...
    corr(electrical_weight(i), ca_weight(i), 'Type', 'Pearson');
[data.spearman.r, data.spearman.p] = ...
    corr(electrical_weight(i), ca_weight(i), 'Type', 'Spearman');
[data.kendall.r, data.kendall.p] = ...
    corr(electrical_weight(i), ca_weight(i), 'Type', 'Kendall');
[pfit, pgof] = fit(electrical_weight(i), ca_weight(i), 'poly1');
data.poly1.fit = pfit;
data.poly1.gof = pgof;
[pfit, pgof] = fit(electrical_weight(i), ca_weight(i), 'poly2');
data.poly2.fit = pfit;
data.poly2.gof = pgof;
stats.electrical = data;

% Compute the total weight statistics.
data = [];
i = ca_table.total_weight > 0;
[data.pearson.r, data.pearson.p] = ...
    corr(total_weight(i), ca_weight(i), 'Type', 'Pearson');
[data.spearman.r, data.spearman.p] = ...
    corr(total_weight(i), ca_weight(i), 'Type', 'Spearman');
[data.kendall.r, data.kendall.p] = ...
    corr(total_weight(i), ca_weight(i), 'Type', 'Kendall');
[pfit, pgof] = fit(total_weight(i), ca_weight(i), 'poly1');
data.poly1.fit = pfit;
data.poly1.gof = pgof;
[pfit, pgof] = fit(total_weight(i), ca_weight(i), 'poly2');
data.poly2.fit = pfit;
data.poly2.gof = pgof;
stats.total = data;
end
