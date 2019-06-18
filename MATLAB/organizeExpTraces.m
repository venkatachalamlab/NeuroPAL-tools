%% Are we smoothing the data?
%is_smooth = true;
is_smooth = false;

%% Get the file names.
files = natsort(split(ls('*_traces.mat')));
files = files(~cellfun('isempty', files));

% Set the fps.
fps = 4.1;

% Set the passband.
if is_smooth
    passband = round([2, 240, 1] * fps); % ~(8,960,4)
end

% Initialize the stimuli.
% #1 = 2-Butanone
% #2 = 2,3-Pentanedione
% #3 = NaCl
load('stimuli.mat');
stim_list = stimulus;
stim_list_keywords = {stim_list.keyword};

% Initialize the light off time (adjustment period).
%light_off_time = 1;
light_off_time = 0.5;

%% Load the experiments, normalize their timing, & clean the traces.
if 1
exps = [];
for i = 1:length(files)
    
    % Load the file.
    disp(['Cleaning traces: ' files{i}]);
    data = load(files{i});

    % Initialize the data variables.
    neuron_names = cellstr(data.neuron_names);
    confs = ones(length(neuron_names),1); %%% FIX ME!!!
    data_times = data.times.times;
    traces = data.gcamp;
    pos = data.positions;
    
    % Remove unlabeled neurons.
    unlabeled = ...
        startsWith(neuron_names, 'Q') | startsWith(neuron_names, 'STITCH');
    neuron_names(unlabeled) = [];
    confs(unlabeled) = [];
    traces(unlabeled,:) = [];
    pos(unlabeled,:,:) = [];
    
    % Normalize the traces & positions.
    [traces, times] = normalizeTimedData(traces, data_times, fps);
    [pos, ~] = normalizeTimedData(pos, data_times, fps);
    
    % Clean up the traces.
    if is_smooth
        [traces, scales, offs] = cleanTraces(traces, fps, 10, 3, [], ...
            'causal', passband);
    else
        [traces, scales, offs] = cleanTraces(traces, fps, 10, 3);
    end
    
    % Store the stimuli.
    stims = data.stimulus;
    for j = 1:length(stims)
        stim_name = lower(stims{j,1});
        iStim = find(cellfun(@(x) contains(stim_name,x), stim_list_keywords));
        switch iStim
            case 1 % 2-Butanone
                exps(i).stims.but(1) = stims{j,2};
                exps(i).stims.but(2) = stims{j,3};
            case 2 % 2,3-Pentanedione
                exps(i).stims.pent(1) = stims{j,2};
                exps(i).stims.pent(2) = stims{j,3};
            case 3 % NaCl
                exps(i).stims.nacl(1) = stims{j,2};
                exps(i).stims.nacl(2) = stims{j,3};
            otherwise % Unknown.
                error('Unknown stimulus: %s', stim_name);
        end
    end

    % Store the traces & positions.
    exps(i).name = files{i};
    exps(i).neuron_names = neuron_names;
    exps(i).confs = confs;
    exps(i).traces = traces;
    exps(i).scales = scales;
    exps(i).offs = offs;
    exps(i).pos = pos;
    exps(i).times = times;
    exps(i).response.light = true;
    exps(i).response.but = true;
    exps(i).response.pent = true;
    exps(i).response.nacl = true;
end
end

%% Reorganize the experiments by neurons & stimuli.
if 1
% Organize the list of experiment names.
exp_names = {exps.name};
num_exps = length(exp_names);

% Aggregate the neurons present.
neuron_names = exps(1).neuron_names;
for i = 2:length(exps)
    neuron_names = union(neuron_names, exps(i).neuron_names);
end
neuron_names = sort(neuron_names);

% Are we being strict about stimulus vs. non-stimulus data?
% Note: the frames during which stimuli are turned on or off average both
% stimulus & non-stimulus responses. To be strict, we should discard these
% "averaged" frames.
is_strict_data = true;

% Initialize memory allocation
nan_exps = nan(1,num_exps);
nan_5s = nan(num_exps, round(fps * 5) + 1);
nan_10s = nan(num_exps, round(fps * 10) + 1);
nan_20s = nan(num_exps, round(fps * 20) + 1);
data_5s = 1:size(nan_5s,2);
data_10s = 1:size(nan_10s,2);
data_20s = 1:size(nan_20s,2);

% Initialize stimuli sizes & indices.
%light_nan = nan_5s;
%light_i = data_5s;
%light_data_i = data_5s;
light_nan = nan_10s;
light_i = data_10s;
light_data_i = data_10s;
but_nan = nan_10s;
but_i = data_10s;
but_data_i = data_10s;
pent_nan = nan_10s;
pent_i = data_10s;
pent_data_i = data_10s;
nacl_nan = nan_20s;
nacl_i = data_20s;
nacl_data_i = data_20s;

% Compute start & end data.
but_start_data_i = but_data_i(1:ceil(length(but_data_i) / 2));
but_end_data_i = but_data_i(floor(length(but_data_i) / 2):end);
pent_start_data_i = pent_data_i(1:ceil(length(pent_data_i) / 2));
pent_end_data_i = pent_data_i(floor(length(pent_data_i) / 2):end);
nacl_start_data_i = nacl_data_i(1:ceil(length(nacl_data_i) / 2));
nacl_end_data_i = nacl_data_i(floor(length(nacl_data_i) / 2):end);

 % Remove stimulus transition frames from data computations.
if is_strict_data
    
    % Recompute stimulus data indices.
    light_data_i = light_data_i(2:(end-1));
    but_data_i = but_data_i(2:(end-1));
    pent_data_i = pent_data_i(2:(end-1));
    nacl_data_i = nacl_data_i(2:(end-1));

    % Recompute stimulus start & end data indices.
    but_start_data_i = but_start_data_i(2:end);
    but_end_data_i = but_end_data_i(1:(end-1));
    pent_start_data_i = pent_start_data_i(2:end);
    pent_end_data_i = pent_end_data_i(1:(end-1));
    nacl_start_data_i = nacl_start_data_i(2:end);
    nacl_end_data_i = nacl_end_data_i(1:(end-1));
end


% Aggregate left/right neurons.
neurons = [];
used_neurons = [];
for i = 1:length(neuron_names)
    
    % Strip the L/R from the neuron name.
    nameLR = neuron_names{i};
    name = stripNeuronLR(nameLR);
        
    % Store the neuron name.
    if isempty(find(strcmp(used_neurons, name), 1))
        used_neurons{end+1} = name;
        neurons(end+1).name = name;
        neurons(end).isLR = ~strcmp(nameLR, name);
        
        % NaN all values.
        if ~neurons(end).isLR
            % NaN the light traces.
            neurons(end).light.stim.trace = light_nan;
            neurons(end).light.post.trace = light_nan;
            % NaN the light response.
            neurons(end).light.stim.power = nan_exps;
            neurons(end).light.post.power = nan_exps;
            neurons(end).light.stim.mean = nan_exps;
            neurons(end).light.post.mean = nan_exps;
            neurons(end).light.off.diff = nan_exps;
            neurons(end).light.off.ratio = nan_exps;
            % NaN the Butanone traces.
            neurons(end).but.pre.trace = but_nan;
            neurons(end).but.stim.trace = but_nan;
            neurons(end).but.post.trace = but_nan;
            % NaN the Butanone response.
            neurons(end).but.pre.power = nan_exps;
            neurons(end).but.stim.power = nan_exps;
            neurons(end).but.post.power = nan_exps;
            neurons(end).but.pre.mean = nan_exps;
            neurons(end).but.stim.mean = nan_exps;
            neurons(end).but.post.mean = nan_exps;
            neurons(end).but.stim.end.mean = nan_exps;
            neurons(end).but.post.start.mean = nan_exps;
            neurons(end).but.on.diff = nan_exps;
            neurons(end).but.on.ratio = nan_exps;
            neurons(end).but.off.diff = nan_exps;
            neurons(end).but.off.ratio = nan_exps;
            % NaN the Pentanedione traces.
            neurons(end).pent.pre.trace = pent_nan;
            neurons(end).pent.stim.trace = pent_nan;
            neurons(end).pent.post.trace = pent_nan;
            % NaN the Pentanedione response.
            neurons(end).pent.pre.power = nan_exps;
            neurons(end).pent.stim.power = nan_exps;
            neurons(end).pent.post.power = nan_exps;
            neurons(end).pent.pre.mean = nan_exps;
            neurons(end).pent.stim.mean = nan_exps;
            neurons(end).pent.post.mean = nan_exps;
            neurons(end).pent.stim.end.mean = nan_exps;
            neurons(end).pent.post.start.mean = nan_exps;
            neurons(end).pent.on.diff = nan_exps;
            neurons(end).pent.on.ratio = nan_exps;
            neurons(end).pent.off.diff = nan_exps;
            neurons(end).pent.off.ratio = nan_exps;
            % NaN the NaCl traces.
            neurons(end).nacl.pre.trace = nacl_nan;
            neurons(end).nacl.stim.trace = nacl_nan;
            neurons(end).nacl.post.trace = nacl_nan;
            % NaN the NaCl response.
            neurons(end).nacl.pre.power = nan_exps;
            neurons(end).nacl.stim.power = nan_exps;
            neurons(end).nacl.post.power = nan_exps;
            neurons(end).nacl.pre.mean = nan_exps;
            neurons(end).nacl.stim.mean = nan_exps;
            neurons(end).nacl.post.mean = nan_exps;
            neurons(end).nacl.stim.end.mean = nan_exps;
            neurons(end).nacl.post.start.mean = nan_exps;
            neurons(end).nacl.on.diff = nan_exps;
            neurons(end).nacl.on.ratio = nan_exps;
            neurons(end).nacl.off.diff = nan_exps;
            neurons(end).nacl.off.ratio = nan_exps;
        else
            %% Left side.
            % NaN the light traces.
            neurons(end).left.light.stim.trace = light_nan;
            neurons(end).left.light.post.trace = light_nan;
            % NaN the light response.
            neurons(end).left.light.stim.power = nan_exps;
            neurons(end).left.light.post.power = nan_exps;
            neurons(end).left.light.stim.mean = nan_exps;
            neurons(end).left.light.post.mean = nan_exps;
            neurons(end).left.light.off.diff = nan_exps;
            neurons(end).left.light.off.ratio = nan_exps;
            % NaN the left Butanone traces.
            neurons(end).left.but.pre.trace = but_nan;
            neurons(end).left.but.stim.trace = but_nan;
            neurons(end).left.but.post.trace = but_nan;
            % NaN the left Butanone response.
            neurons(end).left.but.pre.power = nan_exps;
            neurons(end).left.but.stim.power = nan_exps;
            neurons(end).left.but.post.power = nan_exps;
            neurons(end).left.but.pre.mean = nan_exps;
            neurons(end).left.but.stim.mean = nan_exps;
            neurons(end).left.but.post.mean = nan_exps;
            neurons(end).left.but.stim.end.mean = nan_exps;
            neurons(end).left.but.post.start.mean = nan_exps;
            neurons(end).left.but.on.diff = nan_exps;
            neurons(end).left.but.on.ratio = nan_exps;
            neurons(end).left.but.off.diff = nan_exps;
            neurons(end).left.but.off.ratio = nan_exps;
            % NaN the left Pentanedione traces.
            neurons(end).left.pent.pre.trace = pent_nan;
            neurons(end).left.pent.stim.trace = pent_nan;
            neurons(end).left.pent.post.trace = pent_nan;
            % NaN the left Pentanedione response.
            neurons(end).left.pent.pre.power = nan_exps;
            neurons(end).left.pent.stim.power = nan_exps;
            neurons(end).left.pent.post.power = nan_exps;
            neurons(end).left.pent.pre.mean = nan_exps;
            neurons(end).left.pent.stim.mean = nan_exps;
            neurons(end).left.pent.post.mean = nan_exps;
            neurons(end).left.pent.stim.end.mean = nan_exps;
            neurons(end).left.pent.post.start.mean = nan_exps;
            neurons(end).left.pent.on.diff = nan_exps;
            neurons(end).left.pent.on.ratio = nan_exps;
            neurons(end).left.pent.off.diff = nan_exps;
            neurons(end).left.pent.off.ratio = nan_exps;
            % NaN the left NaCl traces.
            neurons(end).left.nacl.pre.trace = nacl_nan;
            neurons(end).left.nacl.stim.trace = nacl_nan;
            neurons(end).left.nacl.post.trace = nacl_nan;
            % NaN the left NaCl response.
            neurons(end).left.nacl.pre.power = nan_exps;
            neurons(end).left.nacl.stim.power = nan_exps;
            neurons(end).left.nacl.post.power = nan_exps;
            neurons(end).left.nacl.pre.mean = nan_exps;
            neurons(end).left.nacl.stim.mean = nan_exps;
            neurons(end).left.nacl.post.mean = nan_exps;
            neurons(end).left.nacl.stim.end.mean = nan_exps;
            neurons(end).left.nacl.post.start.mean = nan_exps;
            neurons(end).left.nacl.on.diff = nan_exps;
            neurons(end).left.nacl.on.ratio = nan_exps;
            neurons(end).left.nacl.off.diff = nan_exps;
            neurons(end).left.nacl.off.ratio = nan_exps;
            
            %% Right side.
            % NaN the light traces.
            neurons(end).right.light.stim.trace = light_nan;
            neurons(end).right.light.post.trace = light_nan;
            % NaN the light response.
            neurons(end).right.light.stim.power = nan_exps;
            neurons(end).right.light.post.power = nan_exps;
            neurons(end).right.light.stim.mean = nan_exps;
            neurons(end).right.light.post.mean = nan_exps;
            neurons(end).right.light.off.diff = nan_exps;
            neurons(end).right.light.off.ratio = nan_exps;
            % NaN the right Butanone traces.
            neurons(end).right.but.pre.trace = but_nan;
            neurons(end).right.but.stim.trace = but_nan;
            neurons(end).right.but.post.trace = but_nan;
            % NaN the right Butanone response.
            neurons(end).right.but.pre.power = nan_exps;
            neurons(end).right.but.stim.power = nan_exps;
            neurons(end).right.but.post.power = nan_exps;
            neurons(end).right.but.pre.mean = nan_exps;
            neurons(end).right.but.stim.mean = nan_exps;
            neurons(end).right.but.post.mean = nan_exps;
            neurons(end).right.but.stim.end.mean = nan_exps;
            neurons(end).right.but.post.start.mean = nan_exps;
            neurons(end).right.but.on.diff = nan_exps;
            neurons(end).right.but.on.ratio = nan_exps;
            neurons(end).right.but.off.diff = nan_exps;
            neurons(end).right.but.off.ratio = nan_exps;
            % NaN the right Pentanedione traces.
            neurons(end).right.pent.pre.trace = pent_nan;
            neurons(end).right.pent.stim.trace = pent_nan;
            neurons(end).right.pent.post.trace = pent_nan;
            % NaN the right Pentanedione response.
            neurons(end).right.pent.pre.power = nan_exps;
            neurons(end).right.pent.stim.power = nan_exps;
            neurons(end).right.pent.post.power = nan_exps;
            neurons(end).right.pent.pre.mean = nan_exps;
            neurons(end).right.pent.stim.mean = nan_exps;
            neurons(end).right.pent.post.mean = nan_exps;
            neurons(end).right.pent.stim.end.mean = nan_exps;
            neurons(end).right.pent.post.start.mean = nan_exps;
            neurons(end).right.pent.on.diff = nan_exps;
            neurons(end).right.pent.on.ratio = nan_exps;
            neurons(end).right.pent.off.diff = nan_exps;
            neurons(end).right.pent.off.ratio = nan_exps;
            % NaN the right NaCl traces.
            neurons(end).right.nacl.pre.trace = nacl_nan;
            neurons(end).right.nacl.stim.trace = nacl_nan;
            neurons(end).right.nacl.post.trace = nacl_nan;
            % NaN the right NaCl response.
            neurons(end).right.nacl.pre.power = nan_exps;
            neurons(end).right.nacl.stim.power = nan_exps;
            neurons(end).right.nacl.post.power = nan_exps;
            neurons(end).right.nacl.pre.mean = nan_exps;
            neurons(end).right.nacl.stim.mean = nan_exps;
            neurons(end).right.nacl.post.mean = nan_exps;
            neurons(end).right.nacl.stim.end.mean = nan_exps;
            neurons(end).right.nacl.post.start.mean = nan_exps;
            neurons(end).right.nacl.on.diff = nan_exps;
            neurons(end).right.nacl.on.ratio = nan_exps;
            neurons(end).right.nacl.off.diff = nan_exps;
            neurons(end).right.nacl.off.ratio = nan_exps;
        end
    end
end

% Aggregate the data per neuron, stimulus, & experiment.
neuron_names = {neurons.name};
for i = 1:length(exps)
    disp(['Aggregating neurons: ' exps(i).name]);
    names = exps(i).neuron_names;
    for j = 1:length(names)
        
        % Find the neuron.
        disp(['Aggregating neuron: ' exps(i).name ' : ' names{j}]);
        nameLR = names{j};
        name =  stripNeuronLR(nameLR);
        iNeuron = find(strcmp(neuron_names, name));        
        
        % Determine the stimuli bounds.
        light_off = round(light_off_time * fps) + 1;
        but_on = round(exps(i).stims.but(1) * fps) + 1;
        %but_off = round(exps(i).stims.but(2) * fps) + 1;
        pent_on = round(exps(i).stims.pent(1) * fps) + 1;
        %pent_off = round(exps(i).stims.pent(2) * fps) + 1;
        nacl_on = round(exps(i).stims.nacl(1) * fps) + 1;
        %nacl_off = round(exps(i).stims.nacl(2) * fps) + 1;
        
        % Determine the stimuli frames.
        light_stim_i = light_i + light_off - 1; % special case
        light_post_i = light_stim_i(end) + light_i - 1; % special case
        but_pre_i = but_i - but_i(end) + but_on;
        but_stim_i = but_i + but_pre_i(end) - 1;
        but_post_i = but_i + but_stim_i(end) - 1;
        pent_pre_i =  pent_i - pent_i(end) + pent_on;
        pent_stim_i = pent_i + pent_pre_i(end) - 1;
        pent_post_i = pent_i + pent_stim_i(end) - 1;
        nacl_pre_i = nacl_i - nacl_i(end) + nacl_on;
        nacl_stim_i = nacl_i + nacl_pre_i(end) - 1;
        nacl_post_i = nacl_i + nacl_stim_i(end) - 1;
        
        % Determine the stimuli traces.
        trace = squeeze(exps(i).traces(j,:));
        light_trace = trace(light_stim_i);
        light_post_trace = trace(light_post_i);
        but_pre_trace = trace(but_pre_i);
        but_trace = trace(but_stim_i);
        but_post_trace = trace(but_post_i);
        pent_pre_trace = trace(pent_pre_i);
        pent_trace = trace(pent_stim_i);
        pent_post_trace = trace(pent_post_i);
        nacl_pre_trace = trace(nacl_pre_i);
        nacl_trace = trace(nacl_stim_i);
        nacl_post_trace = trace(nacl_post_i);
        
        % Compute the light data.
        light_power = nanmean(light_trace(light_data_i).^2);
        post_light_power = nanmean(light_post_trace(light_data_i).^2);
        light_mean = nanmean(light_trace(light_data_i));
        post_light_mean = nanmean(light_post_trace(light_data_i));
        
        % Compute the butanone data.
        pre_but_power = nanmean(but_pre_trace(but_data_i).^2);
        but_power = nanmean(but_trace(but_data_i).^2);
        post_but_power = nanmean(but_post_trace(but_data_i).^2);
        pre_but_mean = nanmean(but_pre_trace(but_data_i));
        but_mean = nanmean(but_trace(but_data_i));
        post_but_mean = nanmean(but_post_trace(but_data_i));
        but_end_mean = nanmean(but_trace(but_end_data_i));
        post_but_start_mean = nanmean(but_post_trace(but_start_data_i));
        
        % Compute the pentanedione data.
        pre_pent_power = nanmean(pent_pre_trace(pent_data_i).^2);
        pent_power = nanmean(pent_trace(pent_data_i).^2);
        post_pent_power = nanmean(pent_post_trace(pent_data_i).^2);
        pre_pent_mean = nanmean(pent_pre_trace(pent_data_i));
        pent_mean = nanmean(pent_trace(pent_data_i));
        post_pent_mean = nanmean(pent_post_trace(pent_data_i));
        pent_end_mean = nanmean(pent_trace(pent_end_data_i));
        post_pent_start_mean = nanmean(pent_post_trace(pent_start_data_i));
                    
        % Compute the NaCl data.
        pre_nacl_power = nanmean(nacl_pre_trace(nacl_data_i).^2);
        nacl_power = nanmean(nacl_trace(nacl_data_i).^2);
        post_nacl_power = nanmean(nacl_post_trace(nacl_data_i).^2);
        pre_nacl_mean = nanmean(nacl_pre_trace(nacl_data_i));
        nacl_mean = nanmean(nacl_trace(nacl_data_i));
        post_nacl_mean = nanmean(nacl_post_trace(nacl_data_i));
        nacl_end_mean = nanmean(nacl_trace(nacl_end_data_i));
        post_nacl_start_mean = nanmean(nacl_post_trace(nacl_start_data_i));
                    
        % Store the left/right neuron's stimuli traces.
        if neurons(iNeuron).isLR
            
            % Are we dealing with AWC?
            LR = nameLR(end);
            if strncmp('AWC', name, 3) && length(nameLR) > 3
                LR = nameLR(4);
            end
            
            % Store the left/right neuron's stimuli traces.
            switch LR
                case 'L' % Left neuron.
                    neurons(iNeuron).left.conf(i) = exps(i).confs(j);
                    neurons(iNeuron).left.scale(i) = exps(i).scales(j);
                    neurons(iNeuron).left.offset(i) = exps(i).offs(j);
                    % Light.
                    neurons(iNeuron).left.light.stim.trace(i,:) = light_trace;
                    neurons(iNeuron).left.light.post.trace(i,:) = light_post_trace;
                    % Light response.
                    neurons(iNeuron).left.light.stim.power(i) = light_power;
                    neurons(iNeuron).left.light.post.power(i) = post_light_power;
                    neurons(iNeuron).left.light.stim.mean(i) = light_mean;
                    neurons(iNeuron).left.light.post.mean(i) = post_light_mean;
                    neurons(iNeuron).left.light.off.diff(i) = post_light_mean - light_mean;
                    neurons(iNeuron).left.light.off.ratio(i) = post_light_mean / light_mean;
                    % Butanone trace.
                    neurons(iNeuron).left.but.pre.trace(i,:) = but_pre_trace;
                    neurons(iNeuron).left.but.stim.trace(i,:) = but_trace;
                    neurons(iNeuron).left.but.post.trace(i,:) = but_post_trace;
                    % Butanone response.
                    neurons(iNeuron).left.but.pre.power(i) = pre_but_power;
                    neurons(iNeuron).left.but.stim.power(i) = but_power;
                    neurons(iNeuron).left.but.post.power(i) = post_but_power;
                    neurons(iNeuron).left.but.pre.mean(i) = pre_but_mean;
                    neurons(iNeuron).left.but.stim.mean(i) = but_mean;
                    neurons(iNeuron).left.but.post.mean(i) = post_but_mean;
                    neurons(iNeuron).left.but.stim.end.mean(i) = but_end_mean;
                    neurons(iNeuron).left.but.post.start.mean(i) = post_but_start_mean;
                    neurons(iNeuron).left.but.on.diff(i) = but_mean - pre_but_mean;
                    neurons(iNeuron).left.but.on.ratio(i) = but_mean / pre_but_mean;
                    neurons(iNeuron).left.but.off.diff(i) = post_but_mean - pre_but_mean;
                    neurons(iNeuron).left.but.off.ratio(i) = post_but_mean / pre_but_mean;
                    % Pentanedione.
                    neurons(iNeuron).left.pent.pre.trace(i,:) = pent_pre_trace;
                    neurons(iNeuron).left.pent.stim.trace(i,:) = pent_trace;
                    neurons(iNeuron).left.pent.post.trace(i,:) = pent_post_trace;
                    % Pentanedione response.
                    neurons(iNeuron).left.pent.pre.power(i) = pre_pent_power;
                    neurons(iNeuron).left.pent.stim.power(i) = pent_power;
                    neurons(iNeuron).left.pent.post.power(i) = post_pent_power;
                    neurons(iNeuron).left.pent.pre.mean(i) = pre_pent_mean;
                    neurons(iNeuron).left.pent.stim.mean(i) = pent_mean;
                    neurons(iNeuron).left.pent.post.mean(i) = post_pent_mean;
                    neurons(iNeuron).left.pent.stim.end.mean(i) = pent_end_mean;
                    neurons(iNeuron).left.pent.post.start.mean(i) = post_pent_start_mean;
                    neurons(iNeuron).left.pent.on.diff(i) = pent_mean - pre_pent_mean;
                    neurons(iNeuron).left.pent.on.ratio(i) = pent_mean / pre_pent_mean;
                    neurons(iNeuron).left.pent.off.diff(i) = post_pent_mean - pre_pent_mean;
                    neurons(iNeuron).left.pent.off.ratio(i) = post_pent_mean / pre_pent_mean;
                    % NaCl.
                    neurons(iNeuron).left.nacl.pre.trace(i,:) = nacl_pre_trace;
                    neurons(iNeuron).left.nacl.stim.trace(i,:) = nacl_trace;
                    neurons(iNeuron).left.nacl.post.trace(i,:) = nacl_post_trace;
                    % NaCl response.
                    neurons(iNeuron).left.nacl.pre.power(i) = pre_nacl_power;
                    neurons(iNeuron).left.nacl.stim.power(i) = nacl_power;
                    neurons(iNeuron).left.nacl.post.power(i) = post_nacl_power;
                    neurons(iNeuron).left.nacl.pre.mean(i) = pre_nacl_mean;
                    neurons(iNeuron).left.nacl.stim.mean(i) = nacl_mean;
                    neurons(iNeuron).left.nacl.post.mean(i) = post_nacl_mean;
                    neurons(iNeuron).left.nacl.stim.end.mean(i) = nacl_end_mean;
                    neurons(iNeuron).left.nacl.post.start.mean(i) = post_nacl_start_mean;
                    neurons(iNeuron).left.nacl.on.diff(i) = nacl_mean - pre_nacl_mean;
                    neurons(iNeuron).left.nacl.on.ratio(i) = nacl_mean / pre_nacl_mean;
                    neurons(iNeuron).left.nacl.off.diff(i) = post_nacl_mean - pre_nacl_mean;
                    neurons(iNeuron).left.nacl.off.ratio(i) = post_nacl_mean / pre_nacl_mean;
                    
                case 'R' % Right neuron.
                    neurons(iNeuron).right.conf(i) = exps(i).confs(j);
                    neurons(iNeuron).right.scale(i) = exps(i).scales(j);
                    neurons(iNeuron).right.offset(i) = exps(i).offs(j);
                    % Light.
                    neurons(iNeuron).right.light.stim.trace(i,:) = light_trace;
                    neurons(iNeuron).right.light.post.trace(i,:) = light_post_trace;
                    % Light response.
                    neurons(iNeuron).right.light.stim.power(i) = light_power;
                    neurons(iNeuron).right.light.post.power(i) = post_light_power;
                    neurons(iNeuron).right.light.stim.mean(i) = light_mean;
                    neurons(iNeuron).right.light.post.mean(i) = post_light_mean;
                    neurons(iNeuron).right.light.off.diff(i) = post_light_mean - light_mean;
                    neurons(iNeuron).right.light.off.ratio(i) = post_light_mean / light_mean;
                    % Butanone.
                    neurons(iNeuron).right.but.pre.trace(i,:) = but_pre_trace;
                    neurons(iNeuron).right.but.stim.trace(i,:) = but_trace;
                    neurons(iNeuron).right.but.post.trace(i,:) = but_post_trace;
                    % Butanone response.
                    neurons(iNeuron).right.but.pre.power(i) = pre_but_power;
                    neurons(iNeuron).right.but.stim.power(i) = but_power;
                    neurons(iNeuron).right.but.post.power(i) = post_but_power;
                    neurons(iNeuron).right.but.pre.mean(i) = pre_but_mean;
                    neurons(iNeuron).right.but.stim.mean(i) = but_mean;
                    neurons(iNeuron).right.but.post.mean(i) = post_but_mean;
                    neurons(iNeuron).right.but.stim.end.mean(i) = but_end_mean;
                    neurons(iNeuron).right.but.post.start.mean(i) = post_but_start_mean;
                    neurons(iNeuron).right.but.on.diff(i) = but_mean - pre_but_mean;
                    neurons(iNeuron).right.but.on.ratio(i) = but_mean / pre_but_mean;
                    neurons(iNeuron).right.but.off.diff(i) = post_but_mean - pre_but_mean;
                    neurons(iNeuron).right.but.off.ratio(i) = post_but_mean / pre_but_mean;
                    % Pentanedione.
                    neurons(iNeuron).right.pent.pre.trace(i,:) = pent_pre_trace;
                    neurons(iNeuron).right.pent.stim.trace(i,:) = pent_trace;
                    neurons(iNeuron).right.pent.post.trace(i,:) = pent_post_trace;
                    % Pentanedione response.
                    neurons(iNeuron).right.pent.pre.power(i) = pre_pent_power;
                    neurons(iNeuron).right.pent.stim.power(i) = pent_power;
                    neurons(iNeuron).right.pent.post.power(i) = post_pent_power;
                    neurons(iNeuron).right.pent.pre.mean(i) = pre_pent_mean;
                    neurons(iNeuron).right.pent.stim.mean(i) = pent_mean;
                    neurons(iNeuron).right.pent.post.mean(i) = post_pent_mean;
                    neurons(iNeuron).right.pent.stim.end.mean(i) = pent_end_mean;
                    neurons(iNeuron).right.pent.post.start.mean(i) = post_pent_start_mean;
                    neurons(iNeuron).right.pent.on.diff(i) = pent_mean - pre_pent_mean;
                    neurons(iNeuron).right.pent.on.ratio(i) = pent_mean / pre_pent_mean;
                    neurons(iNeuron).right.pent.off.diff(i) = post_pent_mean - pre_pent_mean;
                    neurons(iNeuron).right.pent.off.ratio(i) = post_pent_mean / pre_pent_mean;
                    % NaCl.
                    neurons(iNeuron).right.nacl.pre.trace(i,:) = nacl_pre_trace;
                    neurons(iNeuron).right.nacl.stim.trace(i,:) = nacl_trace;
                    neurons(iNeuron).right.nacl.post.trace(i,:) = nacl_post_trace;
                    % NaCl response.
                    neurons(iNeuron).right.nacl.pre.power(i) = pre_nacl_power;
                    neurons(iNeuron).right.nacl.stim.power(i) = nacl_power;
                    neurons(iNeuron).right.nacl.post.power(i) = post_nacl_power;
                    neurons(iNeuron).right.nacl.pre.mean(i) = pre_nacl_mean;
                    neurons(iNeuron).right.nacl.stim.mean(i) = nacl_mean;
                    neurons(iNeuron).right.nacl.post.mean(i) = post_nacl_mean;
                    neurons(iNeuron).right.nacl.stim.end.mean(i) = nacl_end_mean;
                    neurons(iNeuron).right.nacl.post.start.mean(i) = post_nacl_start_mean;
                    neurons(iNeuron).right.nacl.on.diff(i) = nacl_mean - pre_nacl_mean;
                    neurons(iNeuron).right.nacl.on.ratio(i) = nacl_mean / pre_nacl_mean;
                    neurons(iNeuron).right.nacl.off.diff(i) = post_nacl_mean - pre_nacl_mean;
                    neurons(iNeuron).right.nacl.off.ratio(i) = post_nacl_mean / pre_nacl_mean;
                    
                otherwise % Unknown.
                    error('Unknown left/right neuron: %s', nameLR);
            end
            
        % Store the neuron's stimuli traces.
        else
            neurons(iNeuron).conf(i) = exps(i).confs(j);
            neurons(iNeuron).scale(i) = exps(i).scales(j);
            neurons(iNeuron).offset(i) = exps(i).offs(j);
            % Light.
            neurons(iNeuron).light.stim.trace(i,:) = light_trace;
            neurons(iNeuron).light.post.trace(i,:) = light_post_trace;
            % Light response.
            neurons(iNeuron).light.stim.power(i) = light_power;
            neurons(iNeuron).light.post.power(i) = post_light_power;
            neurons(iNeuron).light.stim.mean(i) = light_mean;
            neurons(iNeuron).light.post.mean(i) = post_light_mean;
            neurons(iNeuron).light.off.diff(i) = post_light_mean - light_mean;
            neurons(iNeuron).light.off.ratio(i) = post_light_mean / light_mean;
            % Butanone.
            neurons(iNeuron).but.pre.trace(i,:) = but_pre_trace;
            neurons(iNeuron).but.stim.trace(i,:) = but_trace;
            neurons(iNeuron).but.post.trace(i,:) = but_post_trace;
            % Butanone response.
            neurons(iNeuron).but.pre.power(i) = pre_but_power;
            neurons(iNeuron).but.stim.power(i) = but_power;
            neurons(iNeuron).but.post.power(i) = post_but_power;
            neurons(iNeuron).but.pre.mean(i) = pre_but_mean;
            neurons(iNeuron).but.stim.mean(i) = but_mean;
            neurons(iNeuron).but.post.mean(i) = post_but_mean;
            neurons(iNeuron).but.stim.end.mean(i) = but_end_mean;
            neurons(iNeuron).but.post.start.mean(i) = post_but_start_mean;
            neurons(iNeuron).but.on.diff(i) = but_mean - pre_but_mean;
            neurons(iNeuron).but.on.ratio(i) = but_mean / pre_but_mean;
            neurons(iNeuron).but.off.diff(i) = post_but_mean - pre_but_mean;
            neurons(iNeuron).but.off.ratio(i) = post_but_mean / pre_but_mean;
            % Pentanedione.
            neurons(iNeuron).pent.pre.trace(i,:) = pent_pre_trace;
            neurons(iNeuron).pent.stim.trace(i,:) = pent_trace;
            neurons(iNeuron).pent.post.trace(i,:) = pent_post_trace;
            % Pentanedione response.
            neurons(iNeuron).pent.pre.power(i) = pre_pent_power;
            neurons(iNeuron).pent.stim.power(i) = pent_power;
            neurons(iNeuron).pent.post.power(i) = post_pent_power;
            neurons(iNeuron).pent.pre.mean(i) = pre_pent_mean;
            neurons(iNeuron).pent.stim.mean(i) = pent_mean;
            neurons(iNeuron).pent.post.mean(i) = post_pent_mean;
            neurons(iNeuron).pent.stim.end.mean(i) = pent_end_mean;
            neurons(iNeuron).pent.post.start.mean(i) = post_pent_start_mean;
            neurons(iNeuron).pent.on.diff(i) = pent_mean - pre_pent_mean;
            neurons(iNeuron).pent.on.ratio(i) = pent_mean / pre_pent_mean;
            neurons(iNeuron).pent.off.diff(i) = post_pent_mean - pre_pent_mean;
            neurons(iNeuron).pent.off.ratio(i) = post_pent_mean / pre_pent_mean;
            % NaCl.
            neurons(iNeuron).nacl.pre.trace(i,:) = nacl_pre_trace;
            neurons(iNeuron).nacl.stim.trace(i,:) = nacl_trace;
            neurons(iNeuron).nacl.post.trace(i,:) = nacl_post_trace;
            % NaCl response.
            neurons(iNeuron).nacl.pre.power(i) = pre_nacl_power;
            neurons(iNeuron).nacl.stim.power(i) = nacl_power;
            neurons(iNeuron).nacl.post.power(i) = post_nacl_power;
            neurons(iNeuron).nacl.pre.mean(i) = pre_nacl_mean;
            neurons(iNeuron).nacl.stim.mean(i) = nacl_mean;
            neurons(iNeuron).nacl.post.mean(i) = post_nacl_mean;
            neurons(iNeuron).nacl.stim.end.mean(i) = nacl_end_mean;
            neurons(iNeuron).nacl.post.start.mean(i) = post_nacl_start_mean;
            neurons(iNeuron).nacl.on.diff(i) = nacl_mean - pre_nacl_mean;
            neurons(iNeuron).nacl.on.ratio(i) = nacl_mean / pre_nacl_mean;
            neurons(iNeuron).nacl.off.diff(i) = post_nacl_mean - pre_nacl_mean;
            neurons(iNeuron).nacl.off.ratio(i) = post_nacl_mean / pre_nacl_mean;
        end
    end
end
end

%% Are we smoothing?
exp_file = 'exp_data.mat';
if is_smooth
    exp_file = 'exp_data_smooth.mat';
end

%% Save the data.
if 1
    save(exp_file, 'exps', 'fps', 'neurons');
end

%% Remove non-responsive worm/stimuli.
if 1
    exps(3).response.nacl = false;
    exps(end).response.pent = false;
    save(exp_file, 'exps', '-append');
end
