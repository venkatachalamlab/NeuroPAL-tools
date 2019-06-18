%%

% Correlations to find
top_N = 10;
bottom_N = 10;
method = "filtered"; % "filtered" or "filtered_diff"

neuron_inclusion_filter = "blacklist"; %"blacklist" or "whitelist"

% Used for blacklisting
prefixes_to_exclude = {"Q", "STITCH"};

% Used for whitelisting
prefixes_to_include = {"A"};

% Filtering to use
LOWPASS_SAMPLES = 30;
HIGHPASS_SAMPLES = 300;

datasets = struct(...
    "animal", {1, 22, 23, 55, 56, 65, 66}, ...
    "run", {402, 101, 201, 301, 401, 301, 101});

animal_idx = 1;

root = "~/Dropbox/180503 23pentanedione_butanone_NaCl/";

deps = fullfile(root, "deps");
addpath(genpath(deps));

%%
dataset = datasets(animal_idx); 
animal = dataset.animal;
part = "head";
run = dataset.run;

animal_name = sprintf("animal_%03d", animal);
dataset_name = sprintf("%s_%s", animal_name, part);
run_name = sprintf("run%03d", run);

dataset_path = fullfile(root, animal_name, dataset_name);
gcamp_file = fullfile(dataset_path, ...
    sprintf("%s_annotations_gcamp.mat", run_name));
filters_path = fullfile(root, "filters");
addpath(filters_path);


S = load(gcamp_file);
annotations = S.annotations_gcamp;

neuron_names = unique(annotations.neuron_id);


if neuron_inclusion_filter == "blacklist"
    names_to_keep = true(size(neuron_names));
    for i = 1:length(prefixes_to_exclude)
        rej = startsWith(string(neuron_names), prefixes_to_exclude{i});
        names_to_keep = names_to_keep & ~rej;
    end
elseif neuron_inclusion_filter == "whitelist"
    names_to_keep = false(size(neuron_names));
    for i = 1:length(prefixes_to_exclude)
        acc = startsWith(string(neuron_names), prefixes_to_include{i});
        names_to_keep = names_to_keep | acc;
    end
end

neuron_names = neuron_names(names_to_keep);

size_N = length(neuron_names);
size_T = max(annotations.t);
G = zeros(size_N, size_T);

for n =1:size_N
    name = neuron_names(n);
    A = annotations(annotations.neuron_id == name, :);
    t_idx = A.t;
    gc_val = A.gcamp;
    G(n, t_idx) = gc_val;
end

clean_trace = @(x) ...
    highpass(...
        lowpass(...
            crop_front(interp_nans(x, 'linear'), 20), ...
            LOWPASS_SAMPLES), ...
         HIGHPASS_SAMPLES);

diff_trace = @(x) [0 sign(diff(x))];


G1 = zeros(size(G), 'double');
G2 = zeros(size(G), 'double');

for i = 1:size(G, 1)
    x = G(i,:);
    x = x - mean(x);
    x = x / quantile(x, 0.9);
    try
        G1(i,:) = clean_trace(x);
    catch err
        G1(i,:) = nan(1, size(G, 2));
    end
    try
        G2(i,:) = diff_trace(clean_trace(x));
    catch
        G2(i,:) = nan(1, size(G, 2));
    end
end

get_pdist = @(x) pdist(x, 'correlation');

% Get pairs


pairs_1 = ranked_pairs_from_pdist(get_pdist(G1));
pairs_2 = ranked_pairs_from_pdist(get_pdist(G2));

% Show pair names.

if method == "filtered"
    pairs = pairs_1;
elseif method == "filtered_diff"
    pairs = pairs_2;
end

get_name = @(x) neuron_names(x);

named_pairs = arrayfun(get_name, pairs);

top_pairs = named_pairs(1:top_N, :)
bottom_pairs = flipud(named_pairs(end-bottom_N+1:end,:))



