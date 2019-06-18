function update_with_run_jpegs(animal_number, animal_part, id_run_number, ...
    id_flip, run_numbers, frame_offsets, frame_counts, annotation_file, ...
    original_jpegs, data_location_in, data_location_out)

if strcmp(animal_part, 'tail')
    id_run = sprintf('id_tail_%03d.nd2', id_run_number);
elseif strcmp(animal_part, 'head')
    id_run = sprintf('id%03d.nd2', id_run_number);
else
    error('Unknown part: %s', animal_part);
end   

runs = {};
for i = 1:length(run_numbers)
    if strcmp(animal_part, 'tail')
        runs{i} = sprintf('run_tail_%03d.nd2', run_numbers(i));
    elseif strcmp(animal_part, 'head')
        runs{i} = sprintf('run%03d.nd2', run_numbers(i));
    else
        error('Unknown part: %s', animal_part);
    end     
end

dead_frames = 2;

dataset_id = sprintf('animal_%03d_%s', animal_number, animal_part);

data_root_in = fullfile(data_location_in, sprintf('animal_%03d', animal_number));
data_root_out = fullfile(data_location_out, sprintf('animal_%03d', animal_number));

output_folder = fullfile(data_root_out, dataset_id);
make_directory(output_folder);

annotator_folder = fullfile(data_location_out, 'annotator_public');
make_directory(annotator_folder);

annotator_image_folder = fullfile(annotator_folder, 'images');
make_directory(annotator_image_folder);

id_shape = [256 512 21];

shapes = {};
for i = 1:length(runs)
    shapes{i} = [128 256 23 floor(frame_counts(i)/23)];
end


%% Get frame 1 from the id run

A = get_array_data(BioFormatsArray(fullfile(data_root_in, id_run)));

gcamp = A(:,:,:,4,1);

gcamp_small = imresize(gcamp, [128 256]);

if id_flip
    gcamp_small = flip(gcamp_small, 1);
    gcamp_small = flip(gcamp_small, 3);
end

%% Get annotations from annotator

all_annotations_struct = load_json(annotation_file);
fns = fieldnames(all_annotations_struct);

all_annotations = all_annotations_struct.(fns{1});
for i = 2:length(fns)
    all_annotations(end+1) = all_annotations_struct.(fns{i});
end

annotations = struct2table(all_annotations(...
    strcmp({all_annotations.dataset_id}, dataset_id)));

translate_coord = @(x, max, scale) x/(scale*max)*(max-1)+1;

annotations.x = translate_coord(annotations.x, 256, 4);
annotations.y = translate_coord(annotations.x, 128, 4);
annotations.z = translate_coord(annotations.z, 21, 8);

%%

crop_dead = @(x) x(:,:,(1+dead_frames):end);
get_clean_slice = @(x,t) crop_dead(get_slice(x,t));

% f = fullfile(data_root, runs{run_idx});
% N = ND2Array(f, shapes{run_idx}, frame_offsets(run_idx));
% a = get_clean_slice(N,1);
% figure; 
% imshowpair(max_intensity_y(gcamp), max_intensity_y(imresize(a, 2)));

%% Find a good range to make uint8

convert_to_uint8 = {};

for run_idx = 1:length(runs)
    quantile_target = 0.9;

    f = fullfile(data_root_in, runs{run_idx});
    N = ND2Array(f, shapes{run_idx}, frame_offsets(run_idx));

    size_T = size(N, 4);

    hi = 0;
    t_idx = 0;
    qs = [];
    for i = 1:10

        t = randi(size_T);

        a = get_clean_slice(N, t);
        a_hi = max_all(a);
        qs(end+1) = quantile(a(:), quantile_target);
        if a_hi > hi
            hi = a_hi;
            t_idx = t;
        end
    end

    lo = max(qs);
    hi = double(hi)*0.9;
    convert_to_uint8{run_idx} = @(x) uint8((double(x)-lo)/(hi-lo)*255);

    for i = 1:10

        t = randi(size_T);

        a = get_clean_slice(N, t);

        az = max_intensity_z(a);
        az8 = convert_to_uint8{run_idx}(az);

        figure(2); clf; imshow(imresize(az8, 4));

        pause;

    end

end

%% Write out uint8 volumes (matlab)

for run_idx = 1:length(runs)

    f = fullfile(data_root_in, runs{run_idx});
    N = ND2Array(f, shapes{run_idx}, frame_offsets(run_idx));

    f_new = fullfile(output_folder, [runs{run_idx}(1:end-3) 'mat']);

    file = matfile(f_new);

    file.Properties.Writable = true;
    S = size(N);
    S(3) = S(3) - dead_frames;
    file.data = zeros(S, 'uint8');
    for t = 1:size(N, 4)
        slice = get_clean_slice(N, t);
        slice = convert_to_uint8{run_idx}(slice);
        if id_flip
            slice = flip(slice, 1);
            slice = flip(slice, 3);
        end
        file.data(:,:,:,t) = slice;
    end

end


%%

make_original_filename = @(view, view_idx, channel_idx, t_idx) fullfile(...
    original_jpegs, ...
    sprintf('%s_%d_%d_%d.jpg', view, view_idx, channel_idx, t_idx));

data_idx = load_json(fullfile(original_jpegs, 'index.json'));

% Update shape to include new volumes
data_idx.shape_t = 1 + length(runs);
data_idx.shape_c = 2;


%% Gather color volumes

image_folder = fullfile(output_folder, 'jpegs');
make_directory(image_folder);

make_filename = @(view, view_idx, channel_idx, t_idx) fullfile(...
    image_folder, ...
    sprintf('%s_%d_%d_%d.jpg', view, view_idx, channel_idx, t_idx));

process_image = @(x) uint8((x - quantile(x(:), 0.9))/2^5);
NP_GCaMP = process_image(gcamp);

write = @(x, view, view_idx, channel_idx, t_idx) imwrite(...
    x, ...
    make_filename(view, view_idx, channel_idx, t_idx), ...
    'Quality', 50);

write_x = @(im, idx, c, t) write(im, 'X', idx, c, t);
write_y = @(im, idx, c, t) write(im, 'Y', idx, c, t);
write_z = @(im, idx, c, t) write(im, 'Z', idx, c, t);

pluck_x = @(im, idx) squeeze(im(:,idx,:));
pluck_y = @(im, idx) squeeze(im(idx,:,:))';
pluck_z = @(im, idx) squeeze(im(:,:,idx));

%

t = 1; % NeuroPAL ID volume
c = 1; % NeuroPAL / MIPs

% Z
for i = 1:data_idx.shape_z

    src = make_original_filename('Z',i,0,0);
    dst = make_filename('Z',i,c,t);
    copyfile(src, dst);

end

% Y
for i = 1:data_idx.shape_y

    src = make_original_filename('Y',i,0,0);
    dst = make_filename('Y',i,c,t);
    copyfile(src, dst);

end

% X
for i = 1:data_idx.shape_x

    src = make_original_filename('X',i,0,0);
    dst = make_filename('X',i,c,t);
    copyfile(src, dst);

end

%
t = 1; % NeuroPAL ID volume
c = 2; % GCaMP

% Z
for i = 1:data_idx.shape_z

    im = pluck_z(NP_GCaMP, i);
    write_z(im, i, c, t);

end

% Y
for i = 1:data_idx.shape_y

    im = pluck_y(NP_GCaMP, i);
    write_y(im, i, c, t);

end

% X
for i = 1:data_idx.shape_x

    im = pluck_x(NP_GCaMP, i);
    write_x(im, i, c, t);

end

% run slices

% Times to take out of runs
run_t_idxs = [1];
t = 1;

for run_idx = 1:length(runs)

    filename = fullfile(output_folder, [runs{run_idx}(1:end-3) 'mat']);
    S = load(filename);
    run_data = 2 * S.data(:,:,:,run_t_idxs);
    clear S;
    
    for j = 1:length(run_t_idxs)
        
        t = t + 1;
        
        A = run_data(:,:,:,j);

        % Z
        mipz = max_intensity_z(A);
        for i = 1:data_idx.shape_z

            im = pluck_z(A, i);
            write_z(im, i, 2, t);
            write_z(mipz, i, 1, t);

        end

        % Y
        mipy = max_intensity_y(A);
        for i = 1:data_idx.shape_y

            im = pluck_y(A, ceil(i/2));
            write_y(im, i, 2, t);
            write_y(mipy, i, 1, t);

        end

        % X
        mipx = max_intensity_x(A);
        for i = 1:data_idx.shape_x

            im = pluck_x(A, ceil(i/2));
            write_x(im, i, 2, t);
            write_x(mipx, i, 1, t);

        end
        
    end
    
end

save_json(data_idx, fullfile(image_folder, 'index.json'));

% Update index


copyfile(image_folder, fullfile(annotator_image_folder, dataset_id));

make_dataset_index(annotator_image_folder);

% Update annotations
all_annotations_struct = load_json(annotation_file);
fns = fieldnames(all_annotations_struct);

output_annotations = struct();
for i = 1:length(fns)
    if strcmp(all_annotations_struct.(fns{i}).dataset_id, dataset_id)
        a = all_annotations_struct.(fns{i});
        a.c = 0;
        a.t = 1;
        output_annotations.(fns{i}) = a;
    end
end

output_filename = fullfile(annotator_folder, 'annotations.json');

if exist(output_filename, 'file')
    current_annotations = load_json(output_filename);
else
    current_annotations = struct();
end

current_annotations = merge_struct(current_annotations, output_annotations);
save_json(current_annotations, output_filename);
