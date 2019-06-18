root = "./";

animal = 66;
part = "head";
run = 101;

animal_name = sprintf("animal_%03d", animal);
dataset_name = sprintf("%s_%s", animal_name, part);
run_name = sprintf("run%03d", run);

%dataset_path = fullfile(root, animal_name, dataset_name);
%annotation_file = fullfile(dataset_path, sprintf("%s_annotations.mat", run_name));
dataset_path = './';
annotation_file = 'run101_annotations.mat';
mips_file = fullfile(dataset_path, sprintf("%s_mip.mat", run_name));

S = load(annotation_file);
annotations = S.annotations;

S = load(mips_file);
images = S.data;

%%
frame = 200;

frame_annotations = annotations(annotations.t==frame,:);
frame_mip = images(:,:,frame);

figure(1); clf; 
imshow(frame_mip);
hold on;
plot(frame_annotations.x, frame_annotations.y, '.');