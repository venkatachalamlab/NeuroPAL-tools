% Put the timing info into the "*_traces.mat" files.

% Get the files.
%files = natsort(split(ls('*/*/*_traces.mat')));
files = natsort(split(ls('*_traces.mat')));
files = files(~cellfun('isempty', files));

% Fix the time in each file.
for i = 1:length(files)
    
    % Construct the timing info filename.
    trace_file = files{i};
    time_file = strrep(trace_file, '_traces.mat', 'times.mat');
    time_file = strrep(time_file, '_tail_', '');
    
    % Find the timing file.
    if ~exist(time_file, 'file')
        warning('Cannot find "%s"!\n', time_file);
        continue;
    end

    % Load the data.
    data = load(time_file);
    stimulus = data.stimulus;
    times = data.times;

    % Compute the mean frames/second.
    fps = 1/nanmedian(diff(times.times));

    % Fix the time in each file.
    disp(['Fixing (' num2str(i) '): "' trace_file '" ...']);
    save(trace_file, 'stimulus', 'fps', 'times', '-append');
end
