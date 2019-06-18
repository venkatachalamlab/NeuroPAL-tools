dataset_id = 'animal_001_tail';

annotation_file = 'C:\Users\vivek\workspace\NeuroPAL\annotator\annotations.json';
data_root = 'G:\My Drive\data\171103 neuropal microfluidics\180503 23pentanedione_butanone_NaCl\animal_001';

id_run = 'id_tail_404.nd2';
flip = false;
runs = {'run_tail_404.nd2', 'run_tail_406.nd2'};

id_shape = [256 512 21];
shape = [128 256 23 floor(23806/23)];

dead_frames = 2;

%% Get frame 1 from the id run

A = get_array_data(BioFormatsArray(fullfile(root, id_run)));

gcamp = A(:,:,:,4,1);

gcamp_small = imresize(gcamp, [128 256]);

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

% [0 4*256] -> [1 256]
% [0 8*21] -> [1 21]
