function annotations_from_tracks(track_folder, name, time_idx, invert, scale)
% ANNOTATIONS_FROM_TRACKS(track_folder, time_idx)
%
%  This creates an annotations.json in the track folder using labeled spots
%  from the first time frame. This file can then be used in the annotator.

if nargin < 4
    scale = 1;
end

if nargin < 3
    invert = false;
end

all_tracks = tracks_from_mamut(track_folder, [1 1 2]);

frame_tracks = all_tracks(all_tracks.time==time_idx, :);

A = BigDataViewerArray(track_folder);

S = size(A);

x = frame_tracks.y;
y = frame_tracks.x;
z = frame_tracks.z;

if invert
    S = S([2 1 3 4]);
    y = S(1) - y + 1;
    z = S(3) - z + 1;
end

x = x*scale(1);
y = y*scale(2);
z = z*scale(3);

annotations = struct();

for i = 1:size(frame_tracks, 1)
    neuron_name = sprintf('Q%03d', 100+i);
    key = sprintf('%s_%s_0', name, neuron_name);
    annotations.(key) = struct(...
        'id', key, ...
        'dataset_id', name, ...
        'neuron_id', neuron_name, ...
        'x', x(i), ...
        'y', y(i), ...
        'z', z(i));
end

json_annotations = jsonencode(annotations);

output_filename = fullfile(track_folder, 'annotations.json');
fid = fopen(output_filename, 'wt');
fprintf(fid, json_annotations);
fclose(fid);