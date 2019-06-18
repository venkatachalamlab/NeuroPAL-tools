trace_root = 'G:\My Drive\data\171103 neuropal microfluidics\nacl_and_23pentanedione\animal_012';
name_to_trace_file = fullfile(trace_root, 'animal_012_head.xlsx');
annotation_file = 'D:\workspace\NeuroPAL\annotator\annotations.json';
dataset_id = 'animal_012_head';
traces = struct(...
    'Run201', load(fullfile(trace_root, 'run201', 'traces.mat'))  );
    %'Run202', load(fullfile(trace_root, 'run202', 'traces.mat')) ...

normalize = @(x) (x - min(x))/((std(x)));
    
trace_fields = fields(traces);
trace_count = length(trace_fields);
trace_lengths = [];

for i = 1:trace_count
    trace_data = traces.(trace_fields{i});
    trace_lengths(i) = size(trace_data.gcamp, 1);
end
trace_linear_idx = [1 cumsum(trace_lengths)+1];
N_points = trace_linear_idx(end) - 1;

lookup_table = readtable(name_to_trace_file);
annotations = load_json(annotation_file);

annotation_ids = fields(annotations);

for i = 1:length(annotation_ids)
    annotation = annotations.(annotation_ids{i});
    if ~strcmp(annotation.dataset_id, dataset_id)
        continue
    end
    neuron_id = annotation.neuron_id;
    rows = strcmp(lookup_table.Annotator, neuron_id);
    lookups = lookup_table(rows, 2:end);
    if size(lookups, 1) == 0
        fprintf('Neuron %s not found\n', neuron_id);
        continue
    elseif size(lookups, 1) > 1
        fprintf('Neuron %s found in multiple rows. Skipping.\n', neuron_id);
        continue
    end
    trace = NaN(N_points, 1);
    for j = 1:trace_count
        start_idx = trace_linear_idx(j);
        end_idx = trace_linear_idx(j+1) - 1;

        lookup_idx = lookups{1, j} + 1;

        calcium = traces.(trace_fields{j}).gcamp;
        trace(start_idx:end_idx) = calcium(:, lookup_idx);
    end

    trace = normalize(trace);
    trace = round(trace, 2);

    annotation.trace = trace;
	annotations.(annotation_ids{i}) = annotation;
end

save_json(annotations, annotation_file);