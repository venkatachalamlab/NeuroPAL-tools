%%

% Correlations to find
% "filtered" or "filtered_diff" or "raw"
method = "filtered";
%method = "raw";

%"blacklist" or "whitelist"
neuron_inclusion_filter = "blacklist"; 

% Used for blacklisting
prefixes_to_exclude = {"Q", "STITCH"};

% Used for whitelisting
prefixes_to_include = {"A"};

% Pixel to micron conversion.
xyz_scale = [0.5,0.5,1.5];

% Normalized time scale.
norm_fps = 4.1;

% Filtering to use
CROP_FRONT_SAMPLES = 20;
LOWPASS_SAMPLES = 10;
HIGHPASS_SAMPLES = 500;
%filter_type = 'causal';
%passband = [2.5, 120, 5]; % [low, high, crop] in seconds
%filter_text = '_2-5';
%filter_type = 'high';
%passband = [1, 2]; % highpass at 1Hz
%filter_text = '_high_1';
%filter_type = 'high';
%passband = [0.5, 1]; % highpass at 2Hz
%filter_text = '_high_0-5';
%filter_type = 'low';
%passband = [5, 10]; % lowpass at 0.2Hz
%filter_text = '_low_5';
filter_type = 'low';
passband = [10, 20]; % lowpass at 0.12Hz
filter_text = '_low_10';

datasets = struct(...
    "animal", {1, 22, 23, 55, 56, 65, 66}, ...
    "run", {402, 101, 201, 301, 401, 301, 101});

%root = "D:\data\180503 23pentanedione_butanone_NaCl\";
root = './';
deps = fullfile(root, "deps");
addpath(genpath(deps));

%% Collect calcium correlations in a table.

data = table('Size', [0, 4], ...
    'VariableTypes', ...
        {'categorical', 'double', 'double', 'int32'}, ...
    'VariableNames', ...
        {'neurons', 'calcium_distance', 'spatial_distance', 'samples'});
data_idx = 0;

wb = waitbar(0);

% Clean the traces -- UNUSED.
clean_trace = @(x) ...
        highpass(...
            lowpass(...
                crop_front(...
                    interp_nans(x, 'linear'), ...
                    CROP_FRONT_SAMPLES), ...
                LOWPASS_SAMPLES), ...
             HIGHPASS_SAMPLES);
diff_trace = @(x) [0 sign(diff(x))];

% Raw, filtered, or derivative traces?
if method == "filtered"
    trace_filter = clean_trace;
elseif method == "filtered_diff"
    trace_filter = @(x) diff_trace(clean_trace(x));
elseif method == "raw"
    trace_filter = @(x) x;
    filter_type = [];
    passband = [];
    filter_text = '';
end

% Get the pairwise coeffcient of variation for traces.
% Compute the correlation coefficient & subtract the squared sum from 1.
% 0 = identical traces & 1 = totally uncorrelated
%get_pdist = @(x) pdist(x, 'correlation');
get_pdist = @(x) pdist(x, @r_squared_dist);
    
% Run pairwise correlations.
N_datasets = length(datasets);
for animal_idx = 1:N_datasets
    waitbar((animal_idx-1)/length(datasets), ...
        wb, "Animal " + string(animal_idx));

    % Construct the filename.
    dataset = datasets(animal_idx); 
    animal = dataset.animal;
    part = "head";
    run = dataset.run;

    % Load the data.
    filename = sprintf('%d_head_run%d_traces.mat', animal, run);
    %dataset_file= fullfile(root, 'Clean_Traces', filename);
    dataset_file = fullfile(root, filename);
    S = load(dataset_file);
    all_neuron_names = categorical(S.neuron_names);

    % Which neurons are we keeping?
    if neuron_inclusion_filter == "blacklist"
        names_to_keep = true(size(all_neuron_names));
        for i = 1:length(prefixes_to_exclude)
            rej = startsWith(string(all_neuron_names), ...
                prefixes_to_exclude{i});
            names_to_keep = names_to_keep & ~rej;
        end
    elseif neuron_inclusion_filter == "whitelist"
        names_to_keep = false(size(all_neuron_names));
        for i = 1:length(prefixes_to_exclude)
            acc = startsWith(string(all_neuron_names), ...
                prefixes_to_include{i});
            names_to_keep = names_to_keep | acc;
        end
    end

    % Remove unwanted neurons.
    neuron_names = all_neuron_names(names_to_keep);
    G = S.gcamp(names_to_keep, :);
    P = S.positions(names_to_keep, :, :);

    % Remove ON/OFF from AWC.
    for j = 1:length(neuron_names)
        name = char(neuron_names(j));
        if strncmp('AWC', name, 3) && length(name) > 4
            neuron_names(j) = categorical(cellstr(name(1:4)));
        end
    end
            
    % Normalize the timing.
    data_times = S.times.times;
    [G, times] = normalizeTimedData(G, data_times, norm_fps);
    [P, ~] = normalizeTimedData(P, data_times, norm_fps);
    
    % Take the mean position.
    P = mean(P, 3) .* xyz_scale;
    
    % Clean the traces.
    filter_band = round(passband * norm_fps);
    if isempty(passband)
        [G,~,~] = cleanTraces(G, norm_fps, 10, 3);
    else
        [G,~,~] = cleanTraces(G, norm_fps, 10, 3, [], filter_type, filter_band);
    end
    
    % Center & scale the trace.
    G_filtered = G;
    if 0
        G_filtered = zeros(size(G), 'double');
        for i = 1:size(G, 1)
            x = G(i,:);
            x = x - mean(x);
            x = x / quantile(x, 0.9);
            try
                G_filtered(i,:) = trace_filter(x);
            catch err
                G_filtered(i,:) = nan(1, size(G, 2));
            end
        end
    end

    % Get neuron pairs ranked by correlation distance.
    pairs = ranked_pairs_from_pdist(get_pdist(G_filtered));
    get_name = @(x) neuron_names(x);
    named_pairs = arrayfun(get_name, pairs(:,1:2));
        
    % Organize the cross correlations into a table.
    run_data = table('Size', [size(pairs,1), 4], ...
        'VariableTypes', ...
            {'categorical', 'double', 'double', 'int32'}, ...
        'VariableNames', ...
            {'neurons', 'calcium_distance', 'spatial_distance', 'samples'});
    N_pairs = length(pairs);
    wb_offset = (animal_idx-1)/length(datasets); % progress bar
    for i = 1:length(pairs)
        waitbar(wb_offset +  i/N_pairs/N_datasets, wb); drawnow;
        
        n1 = named_pairs(i,1);
        n2 = named_pairs(i,2);

        % Alphabetize keys to make pair comparison easier
        if string(n1) < string(n2)
            run_data{i,'neurons'} = string(n1) + ":" + string(n2);
        else
            run_data{i,'neurons'} = string(n2) + ":" + string(n1);
        end

        run_data{i, 'calcium_distance'} = pairs(i, 3);
        run_data{i, 'spatial_distance'} = ...
            norm(P(pairs(i,1),:) - P(pairs(i,2),:));

    end
    
    % Organize the data.
    data = [data;run_data];
end

% Compute the mean & SD for the activity & distance correlations.
[G, pair_names] = findgroups(data.neurons);
calcium_distance_mean = splitapply(@mean, data.calcium_distance, G);
spatial_distance_mean = splitapply(@mean, data.spatial_distance, G);
calcium_distance_std = splitapply(@std, data.calcium_distance, G);
spatial_distance_std = splitapply(@std, data.spatial_distance, G);
samples = splitapply(@length, data.calcium_distance, G);

ca_table = table(pair_names, samples, calcium_distance_mean, ...
    calcium_distance_std, spatial_distance_mean, spatial_distance_std);
ca_table.Properties.VariableNames = {'connection', 'samples', ...
    'calcium_dist', 'calcium_std', 'spatial_dist', 'spatial_std'};

%% Add synapses.

% Load the connectome.
S = load("connectome_syn_gj.mat");
edges = S.conn_edges;

% % Organize the structural connectome for comparison.
% wb = waitbar(0);
% weight = zeros(size(ca_table, 1), 1);
% synapse_type = repmat("none", [size(ca_table, 1), 1]);
% for i = 1:size(edges, 1)
%     
%     waitbar(i/size(edges,1), wb);
%     
%     % Get the strauctural connectome edge.
%     row = edges(i,:);
%     n1 = row.Source;
%     n2 = row.Target;
%     
%     % Alphabetize keys to make pair comparison easier.
%     if string(n1) < string(n2)
%         pairname = string(n1) + ":" + string(n2);
%     else
%         pairname = string(n2) + ":" + string(n1);
%     end
%     
%     % Find the functional connectome equivalent edge.
%     % Note 1: since the functional connectome is undirected, the structural
%     % connectome directionality will be degenerate when matching.
%     % Note 2: since we're indexing using the structural connectome, we'll hit
%     % same pair chemical & electrical edges independently & ensure they
%     % both have representation.
%     row_idx = find(ca_table.pair_names == pairname);
%     if row_idx
%         weight(row_idx) = row.Weight;
%         synapse_type(row_idx) = string(row.Type);
%     end    
% end
% delete(wb);
% % Add the synapse types & weights to the joint connectome.
% ca_table.synapse_type = categorical(synapse_type);
% ca_table.synapse_weight = weight;


%% Build the calcium connectome.

% Remove 0-padding from the connectome.
for i = 1:size(edges,1)
    
    % Fix the source.
    name = char(edges.Source(i));
    if length(name) > 3 && name(3) == '0'
        edges.Source(i) = categorical(cellstr(name([1:2,4:end])));
    end
    
    % Fix the target.
    name = char(edges.Target(i));
    if length(name) > 3 && name(3) == '0'
        edges.Target(i) = categorical(cellstr(name([1:2,4:end])));
    end
end

% Initialize the calcium connectome.
chemical = categorical("chemical");
electrical = categorical("electrical");
ca_table.chemical_weight(:) = 0;
ca_table.electrical_weight(:) = 0;
ca_table.total_weight(:) = 0;
ca_table.is_mNTR(:) = false;
ca_table.neuron1(:) = categorical("");
ca_table.neuron1_class(:) = categorical("");
ca_table.neuron2(:) = categorical("");
ca_table.neuron2_class(:) = categorical("");

% Build the calcium connectome.
wb = waitbar(0);
for i = 1:size(ca_table, 1)
    
    % Update the progress bar.
    waitbar(i/size(ca_table,1), wb, 'Building calcium connectome.');
    
    % Determine the neuron pair.
    pair = cellstr(ca_table.connection(i));
    neurons = split(pair, ':');
    neuron1 = categorical(cellstr(neurons{1}));
    neuron2 = categorical(cellstr(neurons{2}));
    ca_table.neuron1(i) = neuron1;
    ca_table.neuron2(i) = neuron2;
    
    % Determine the neuron classes.
    ca_table.neuron1_class(i) = ...
        edges.source_class(find(edges.Source == neuron1, 1));
    ca_table.neuron2_class(i) = ...
        edges.source_class(find(edges.Source == neuron2, 1));
    
    % Find the structural connectome's equivalent edges.
    pair_i = find((edges.Source == neuron1 & edges.Target == neuron2) | ...
        (edges.Source == neuron2 & edges.Target == neuron1));

    % Compute the connection's properties.
    is_mNTR = false;
    chemical_weight = 0;
    electrical_weight = 0;
    for j = 1:length(pair_i)
        k = pair_i(j);
        switch edges.Type(k)
            case chemical
                is_mNTR = is_mNTR | edges.mNTR_edge(k);
                chemical_weight = chemical_weight + edges.Weight(k);
            case electrical
                electrical_weight = electrical_weight + edges.Weight(k);
        end
    end
    ca_table.is_mNTR(i) = is_mNTR;
    ca_table.chemical_weight(i) = chemical_weight;
    ca_table.electrical_weight(i) = electrical_weight;
    ca_table.total_weight(i) = chemical_weight + electrical_weight;
end

%% Find nearest spatial neighbors (controls) in the calcium connectome.

% Initialize the nearest neighbors.
ca_table.nearest_neuron1(:) = categorical("");
ca_table.nearest_neuron1_calcium(:) = 0;
ca_table.nearest_neuron1_spatial(:) = 0;
ca_table.nearest_neuron2(:) = categorical("");
ca_table.nearest_neuron2_calcium(:) = 0;
ca_table.nearest_neuron2_spatial(:) = 0;

% Find the nearest neighbors.
wb = waitbar(0);
for i = 1:size(ca_table, 1)
    
    % Update the progress bar.
    waitbar(i/size(ca_table,1), wb, 'Finding nearest neghbors.');
    
    % Find neuron 1's nearest neighbor.
    neuron1 = ca_table.neuron1(i);
    near_i = find((ca_table.neuron1 == neuron1) | (ca_table.neuron2 == neuron1));
    [min_dist, min_i] = min(ca_table.spatial_dist(near_i));
    j = near_i(min_i);
    nearest_neuron = ca_table.neuron1(j);
    if neuron1 == nearest_neuron
        nearest_neuron = ca_table.neuron2(j);
    end
    ca_table.nearest_neuron1(i) = nearest_neuron;
    ca_table.nearest_neuron1_calcium(i) = ca_table.calcium_dist(j);
    ca_table.nearest_neuron1_spatial(i) = ca_table.spatial_dist(j);
    
    % Find neuron 2's nearest neighbor.
    neuron2 = ca_table.neuron2(i);
    near_i = find((ca_table.neuron1 == neuron2) | (ca_table.neuron2 == neuron2));
    [min_dist, min_i] = min(ca_table.spatial_dist(near_i));
    j = near_i(min_i);
    nearest_neuron = ca_table.neuron1(j);
    if neuron2 == nearest_neuron
        nearest_neuron = ca_table.neuron2(j);
    end
    ca_table.nearest_neuron2(i) = nearest_neuron;
    ca_table.nearest_neuron2_calcium(i) = ca_table.calcium_dist(j);
    ca_table.nearest_neuron2_spatial(i) = ca_table.spatial_dist(j);
end
delete(wb);

% Sort the table.
ca_table = sortrows(ca_table, "calcium_dist");

% Save the table.
save('ca_conn_' + method + filter_text, 'ca_table');

%% Compute the correlation coefficient.

% Get the ranks for each connectome.
ca_filt = ca_table.samples > 3 & ... % & ca_table.chemical_weight > 0
    (ca_table.neuron1_class ~= ca_table.neuron2_class);
ca = 1 - ca_table.calcium_dist(ca_filt);
syn = ca_table.total_weight(ca_filt);

% Compute the correlation between the connectomes.
% TODO: organize by edge type & use Spearmann!!!
corr_matrix = corrcoef(ca, syn);
r = corr_matrix(2)
[r, p] = corr(ca, syn,'Type','Spearman')

% Scatter plot the correlation.
figure(1); clf; plot(ca, syn, '.')

% Compute the regression then plot the predicted synaptic weights.
X = [ones(length(ca),1) ca];
b = X\syn;
syn_pred = X*b;
figure(1); hold on; plot(ca, syn_pred, 'LineWidth', 2);

% Label the plots.
xlabel("Distance between calcium traces (1-r^2)");
ylabel("Synapse weight");
title(["Correlation coefficient: " num2str(r)]);
