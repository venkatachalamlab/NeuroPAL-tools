% if getenv('computername') == "DESKTOP-FG5FH37"
%     root_drive = 'D:\';
% elseif getenv('computername') == "DESKTOP-VGB0UKT"
%     root_drive = 'Z:\';
% end

root_drive = '~/';
np_annotation_file = '~/Documents/Github/NeuroPAL/annotator/annotations.json';

local_root = '~/Dropbox/180503 23pentanedione_butanone_NaCl';
% remote_root = fullfile(root_drive, 'Dropbox', 'data', ...
%     '180503 23pentanedione_butanone_NaCl');
%np_root = fullfile(root_drive, 'workspace', 'NeuroPAL', 'matlab');
%vv_tools = fullfile(root_drive, 'workspace', 'vv_tools', 'MATLAB');
%filters_path = fullfile(local_root, "filters");

%addpath(genpath(np_root));
%addpath(genpath(vv_tools));
%addpath(genpath(filters_path));

neuron_inclusion_filter = "blacklist"; %"blacklist" or "whitelist"

% Used for blacklisting
%prefixes_to_exclude = {"Q", "STITCH"};
prefixes_to_exclude = {};

% Used for whitelisting
prefixes_to_include = {"A"};

% Filtering to use
LOWPASS_SAMPLES = 30;
HIGHPASS_SAMPLES = 300;

datasets = struct(...
    "animal", {1, 22, 23, 55, 56, 65, 66}, ...
    "run", {402, 101, 201, 301, 401, 301, 101}, ...
    "id_volume", {nan, nan, nan, "id301.ome.tiff", "id402.ome.tiff", nan, "id101.ome.tiff"});

load_status;

%% Prepare color volume

dataset_idx = 7;
dataset = datasets(dataset_idx);

animal = dataset.animal;
part = "head";
run = dataset.run;

animal_name = sprintf("animal_%03d", animal);
dataset_name = sprintf("%s_%s", animal_name, part);
run_name = sprintf("run%03d", run);

dataset_path = fullfile(local_root, animal_name, dataset_name);
gcamp_file = fullfile(dataset_path, ...
    sprintf("%s_annotations_gcamp.mat", run_name));
trace_file = fullfile(dataset_path, ...
    sprintf("%s_traces.mat", run_name));

color_data_file = fullfile(dataset_path, dataset.id_volume);
color_data = load_tiff_stack(color_data_file);
S = size(color_data);
color_data = reshape(color_data, [S(1), S(2), 21, 4]);

R = color_data(:,:,:,1);
G = color_data(:,:,:,2);
B = color_data(:,:,:,3);

% Clean up haze in red channel (high pass filter)
R64 = double(R)/2^16;
R64 = R64 - smooth3(R64, 'gaussian', [31 31 11], 15);
R = uint16(2^16*R64);
% done.

threshold = @(x, p) x - quantile(x(:), p);
median_filter = @(x, Nx, Nz) medfilt3(x, [Nx, Nx, Nz]);
smooth = @(x) smooth3(double(x)/double(quantile(x(:), 0.9999)), 'gaussian');
process = @(x) uint8(255*smooth(median_filter(threshold(x, 0.7), 5, 3)));
mipz = @max_intensity_z;

combine = @(R, G, B) cat(3, mipz(R), mipz(G), mipz(B));
show3D = @(x) imshow(mipz(x));
show3DRGB = @(R, G, B) imshow(cat(3, mipz(R), mipz(G), mipz(B)));

colorscale = [1.5, 1.0, 1.5];

R_out = process(R) * colorscale(1);
G_out = process(G) * colorscale(1);
B_out = process(B) * colorscale(1);

figure(6); clf; show3DRGB(R_out, G_out, B_out);

%% Create a mask to select points to animate.

X = R_out + G_out + B_out;
mask = uint8(X > quantile(X(:), 0.96));
R_masked = R_out.*mask;
G_masked = G_out.*mask;
B_masked = B_out.*mask;

figure(7); clf;  show3DRGB(R_masked, G_masked, B_masked);

%%

S = load(gcamp_file);
TF = load(trace_file);
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

if 0
for i = 1:size(G, 1)
    x = G(i,:);
    x = x - mean(x, 'omitnan');
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
else
    G1 = cleanTraces(G, TF.fps);
end

%% Get Ev's annotations
all_annotations_struct = load_json(np_annotation_file);
fns = fieldnames(all_annotations_struct);

all_annotations_struct_out = struct();
for i = 1:length(fns)
    if strcmp(all_annotations_struct.(fns{i}).dataset_id, dataset_name)
        a = all_annotations_struct.(fns{i});
        a.c = 0;
        a.t = 1;
        all_annotations_struct_out.(fns{i}) = a;
    end
end

fns = fieldnames(all_annotations_struct_out);

all_annotations_struct_array = all_annotations_struct_out.(fns{1});
for i = 2:length(fns)
    all_annotations_struct_array(end+1) = all_annotations_struct_out.(fns{i});
end

all_annotations = struct2table(all_annotations_struct_array);

all_annotations.neuron_id = categorical(all_annotations.neuron_id);
all_annotations.dataset_id = categorical(all_annotations.dataset_id);

translate_coord = @(x, max, scale) x/(scale*max)*(max-1)+1;

all_annotations.x = translate_coord(all_annotations.x, 256, 2);
all_annotations.y = translate_coord(all_annotations.y, 128, 2);
all_annotations.z = translate_coord(all_annotations.z, 21, 8);

all_annotations.Properties.RowNames = all_annotations.id;

all_id_annotations = all_annotations(all_annotations.t==1, :);

figure(7); hold on; plot(all_id_annotations.x, all_id_annotations.y, '.');

%% Assign each pixel to the nearest neuron

X = all_id_annotations{:,{'y', 'x', 'z'}};
DT = delaunayTriangulation(X);
T = DT.ConnectivityList;
I_linear = find(mask);
[I, J, K] = ind2sub(size(mask), I_linear);
XI = [I, J, K];
IDX = dsearchn(X, T, XI);

IDX_mask = mask;
IDX_mask(I_linear) = IDX;

figure(8); clf; show3D(IDX_mask);

%% Remove mask points that aren't near neurons.

R_max = 10;
metric = diag([1, 1, 4]);
mask_2 = mask;
for i = 1:numel(mask)
    neuron_idx = IDX_mask(i);
    if neuron_idx == 0
        continue;
    end
    n = X(neuron_idx, :); % neuron coordinate
    [y, x, z] = ind2sub(size(mask), i);
    p = [y, x, z]; % point coordinate
    d = p-n;
    d2 = d*metric*d';
    if d2 > R_max^2
        mask_2(i) = 0;
    end
end

IDX_mask_2 = IDX_mask .* mask_2;

figure(9); clf; show3D(mask_2.*IDX_mask);

%% Modulate color intensity using the neuronal traces.

F = 1.2.^(G1);

size_T = size(G, 2);
A = zeros([size(R_masked), 3, size_T], 'uint8');
A(:,:,:,1,:) = repmat(R_masked, [1 1 1 1 size_T]);
A(:,:,:,2,:) = repmat(G_masked, [1 1 1 1 size_T]);
A(:,:,:,3,:) = repmat(B_masked, [1 1 1 1 size_T]);

wb = waitbar(0, "Generating animated array.");
for i = 1:numel(mask_2)
    neuron_idx = IDX_mask_2(i);
    if neuron_idx == 0
        continue;
    end
    
    [y, x, z] = ind2sub(size(mask_2), i);
    neuron_name = all_id_annotations.neuron_id(neuron_idx);
    neuron_gcamp_idx = find(neuron_names == neuron_name);
    
    if neuron_gcamp_idx
        F1 = F(neuron_gcamp_idx, :);
        F3 = cat(1, F1, F1, F1);
        F3 = reshape(F3, [1, 1, 1, size(F3)]);

        A(y,x,z,:,:) = uint8(double(A(y,x,z,:,:)) .* F3);
    end
    percent_complete = i/numel(mask_2)*100;
    waitbar(percent_complete, wb, "Generating animated array.");
end

%%

Az = max_intensity_z(A);
Rz = squeeze(Az(:,:,1,:));
Gz = squeeze(Az(:,:,2,:));
Bz = squeeze(Az(:,:,3,:));

movie_file = sprintf("%s_animated_with_run_%03d.mp4", ...
    color_data_file, dataset.run);

movie_from_colors(Rz, Gz, Bz, movie_file);

