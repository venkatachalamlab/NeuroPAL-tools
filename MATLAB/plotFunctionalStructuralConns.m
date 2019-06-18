%% Plot the functional vs. structural connectomes.

% Load the functional connectome.
load('ca_conn_raw.mat');

% Remove low sampling.
N_thresh = 3;
ca_table = ca_table(ca_table.samples > N_thresh,:);

% Determine the chemical/electrical synapses.
chem_i = true(size(ca_table,1),1);
elec_i = true(size(ca_table,1),1);
%chem_i = ca_table.chemical_weight > 0;
%elec_i = ca_table.electrical_weight > 0;

% Determine the L/R neurons.
LR_i = ca_table.neuron1_class == ca_table.neuron2_class;

% Determine the nearby neurons.
near_thresh = 5; % 5um is near.
near_i = ca_table.spatial_dist <= near_thresh;

% Determine the non L/R or near neurons.
other_i = ~LR_i & ~near_i;

% Are we scaling the graph?
scale = @(x) (x);
%scale = @(x) log(x+1);
%scale = @(x) log10(x+1);

%% Compute the correlations.
stats = [];

% Compute the correlation for all.
stats.all = computeFunctStructStats(ca_table);

% Compute the correlation excluding L/R neurons.
stats.no_LR = computeFunctStructStats(ca_table(~LR_i,:));

% Compute the correlation excluding nearby neurons.
stats.no_near = computeFunctStructStats(ca_table(~near_i,:));

% Compute the correlation excluding L/R neurons.
stats.no_LR_near = computeFunctStructStats(ca_table(~LR_i & ~near_i,:));

%% Plot the data.

% Setup the plot markers.
LR_mark = 's';
near_mark = '^';
other_mark = 'o';
mark_alpha = 0.5;
chem_mark_size = 60;
elec_mark_size = 30;
chem_line_width = 0.5;
elec_line_width = 1.5;
chem_edge_color = [1,1,1];
elec_edge_color = [0,0,0];
chem_fill_color = [0,0,0];
elec_fill_color = [1,1,1];
line_width = 2.5;
chem_line = '-';
elec_line = ':';

% Plot the chemical connections.
figure;
hold on;
i = chem_i & LR_i;
scatter(scale(ca_table.chemical_weight(i)), ...
    scale(1 - ca_table.calcium_dist(i)), chem_mark_size, ...
    LR_mark, 'MarkerFaceColor', chem_fill_color, ...
    'MarkerEdgeColor', chem_edge_color, 'MarkerFaceAlpha', mark_alpha, ...
    'LineWidth', chem_line_width);
i = chem_i & near_i;
scatter(scale(ca_table.chemical_weight(i)), ...
    scale(1 - ca_table.calcium_dist(i)), chem_mark_size, ...
    near_mark, 'MarkerFaceColor', chem_fill_color, ...
    'MarkerEdgeColor', chem_edge_color, 'MarkerFaceAlpha', mark_alpha, ...
    'LineWidth', chem_line_width);
i = chem_i & other_i;
scatter(scale(ca_table.chemical_weight(i)), ...
    scale(1 - ca_table.calcium_dist(i)), chem_mark_size, ...
    other_mark, 'MarkerFaceColor', chem_fill_color, ...
    'MarkerEdgeColor', chem_edge_color, 'MarkerFaceAlpha', mark_alpha, ...
    'LineWidth', chem_line_width);

% Plot the electrical connections.
i = elec_i & LR_i;
scatter(scale(ca_table.electrical_weight(i)), ...
    scale(1 - ca_table.calcium_dist(i)), elec_mark_size, ...
    LR_mark, 'MarkerFaceColor', elec_fill_color, ...
    'MarkerEdgeColor', elec_edge_color, 'MarkerFaceAlpha', mark_alpha, ...
    'LineWidth', elec_line_width);
i = elec_i & near_i;
scatter(scale(ca_table.electrical_weight(i)), ...
    scale(1 - ca_table.calcium_dist(i)), elec_mark_size, ...
    near_mark, 'MarkerFaceColor', elec_fill_color, ...
    'MarkerEdgeColor', elec_edge_color, 'MarkerFaceAlpha', mark_alpha, ...
    'LineWidth', elec_line_width);
i = elec_i & other_i;
scatter(scale(ca_table.electrical_weight(i)), ...
    scale(1 - ca_table.calcium_dist(i)), elec_mark_size, ...
    other_mark, 'MarkerFaceColor', elec_fill_color, ...
    'MarkerEdgeColor', elec_edge_color, 'MarkerFaceAlpha', mark_alpha, ...
    'LineWidth', elec_line_width);

% Determine the data limits.
x_lim = [0, max(max(ca_table.chemical_weight), ...
    max(ca_table.electrical_weight))];
y_lim = [0, max(1 - ca_table.calcium_dist)];

% Plot the chemical correlation line.
struc_x = scale(ca_table.chemical_weight(chem_i));
funct_y = scale(1 - ca_table.calcium_dist(chem_i));
%chem_fit = polyfit(struc_x, funct_y);
X = [ones(length(struc_x),1), struc_x];
b = X\funct_y;
funct_pred = X*b;
struc_x = struc_x(1,end);
funct_pred = funct_pred(1,end);
plot(struc_x, funct_pred, 'LineWidth', line_width, 'LineStyle', chem_line);

% Plot the fits.
fplot = plot(stats.all.electrical.poly1.fit);
fplot.LineWidth = line_width;
fplot.LineStyle = elec_line;
fplot = plot(stats.all.chemical.poly1.fit);
fplot.LineWidth = line_width;
fplot.LineStyle = chem_line;

% Label the plot.
xlabel('Structural Weight (Counts)');
ylabel('Functional Weight (Correlation)');
legend off;

% Set the plot limits.
xlim(x_lim);
ylim(y_lim);
box on;

% Set the plot to log scale.
%gca.xscale = 'log';
%gca.yscale = 'log';
