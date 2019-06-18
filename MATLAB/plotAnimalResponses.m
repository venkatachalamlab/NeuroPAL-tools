function plotAnimalResponses(neuron_data, names, exp_i, stim_i, fps, varargin)
%PLOTNEURONRESPONSES Plot the responses for a set of neurons.
%
%   Input:
%       exps        - the experiments
%       neuron_data - an array of neurons structs with the info
%       names       - a list of names of neurons to plot
%       exp_i       - the index/neuron of the experiment to plot
%       stim_i      - the index/neuron of the stimulus to plot
%                     1 = butanone
%                     2 = pentanedione
%                     3 = NaCl
%                     4 = light
%       fps         - frames/second
%       is_rescale  - rescale the traces to [0.05,0.95] for visibility?
%          default: false
%       plot_mode   - how should we plot the neural traces?
%          default: 3
%          1 = plot the neural traces overlayed
%          2 = plot the neural traces as mean & SD
%          3 = plot the neural traces as mean & SEM

% Are we rescaling the traces to [0.05,0.95] for visibility?
is_rescale = false;
if ~isempty(varargin)
    is_rescale = (varargin{1});
end

% How are we plotting the neurons?
plot_mode = 3;
if length(varargin) > 1 && ~isempty(varargin{2})
    plot_mode = (varargin{2});
end

% Initialize the neuron names.
neuron_names = {neuron_data.name};
if ischar(names)
    names = {names};
end

% Initialize the experiment list.
if ~iscell(exp_i)
    exp_i = {exp_i};
end

% Initialize the stimulus list.
if ~iscell(stim_i)
    stim_i = {stim_i};
end

% Initialize the experiment data.
neuron = neuron_data(1);
if neuron.isLR
    all_exps = 1:length(neuron.left.but.pre.mean);
else
    all_exps = 1:length(neuron.but.pre.mean);
end

% Initialize the stimulus data.
% Note: we'll throw an error if anything is out of range.
stim_name = ['B', 'P', 'N', 'L'];
all_stims = 1:length(stim_name);
stims_used = [];

% Initialize the trace size data.
% Note: pre/stim/post traces have overlapping data. We drop the transition
% frames to avoid duplicating data in the trace.
light_trace_size = 2 * (round(10 * fps) + 1) - 1; % drop 1 frame
odor_trace_size = 3 * (round(10 * fps) + 1) - 2; % drop 2 frames
nacl_trace_size = 3 * (round(20 * fps) + 1) - 2; % drop 2 frames

% Compute the odor trace indices for NaCl.
odor_trace_i = 1:odor_trace_size;
odor_nacl_i = odor_trace_i + round(10 * fps);

% Compute the light trace indices for odors & NaCl.
light_trace_i = 1:light_trace_size;
light_odor_i = light_trace_i + round(10 * fps);
light_nacl_i = light_trace_i + round(20 * fps);

% Assemble the neural responses.
neurons = cell(length(names), 1);
traces = cell(length(names), 1);
for i = 1:length(names)
    
    % Determine the neuron's name & index.
    nameLR = names{i};
    name = stripNeuronLR(nameLR);
    neuron_i = find(strcmp(neuron_names, name));
    
    % Determine the neuron's side.
    LR = [];
    name_end = length(name) + 1;
    if name_end <= length(nameLR) 
        LR = nameLR(name_end);
    end
    
    % Use all experiments.
    if isempty(exp_i{i})
        exp_i{i} = all_exps;
    end
    exps_i = exp_i{i};
    
    % Use all stimuli.
    if isempty(stim_i{i})
        stim_i{i} = all_stims;
    end
    stims_i = stim_i{i};
    stim_str = nan(1,length(stims_i));
    
    % Store the stims used.
    stims_used = union(stims_used, stims_i);
    
    % Add the neural response.
    if ~neuron_data(neuron_i).isLR

        % Assemble the stimulus responses.
        for j = 1:length(stims_i)
            
            % Determine the stimulus.
            stim_str(j) = stim_name(stims_i(j));
            switch stims_i(j)
                case 1
                    data = neuron_data(neuron_i).but;
                case 2
                    data = neuron_data(neuron_i).pent;
                case 3
                    data = neuron_data(neuron_i).nacl;
                case 4
                    data = neuron_data(neuron_i).light;
            end
            
            % Assemble the response.
            if stims_i(j) <  4
                trace_dataLR = [data.pre.trace(exps_i,:), ...
                    data.stim.trace(exps_i,2:(end-1)), ...
                    data.post.trace(exps_i,:)];
                
            % Assemble the light response.
            else
                trace_dataLR = [data.stim.trace(exps_i,1:(end-1)), ...
                    data.post.trace(exps_i,:)];
            end
            
            % Assemble the light trace data.
            if any(stims_i < 4) && stims_i(j) ==  4
                
                % We're plotting NaCl.
                if any(stims_i == 3)
                    trace_data = nan(length(exps_i), nacl_trace_size);
                    trace_data(:,light_nacl_i) = trace_dataLR;
                    
                % We're plotting odors.
                else
                    trace_data = nan(length(exps_i), odor_trace_size);
                    trace_data(:,light_odor_i) = trace_dataLR;
                end
                
            % Assemble the odor trace data.
            elseif any(stims_i == 3) && stims_i(j) < 3
                trace_data = nan(length(exps_i), nacl_trace_size);
                trace_data(:,odor_nacl_i) = trace_dataLR;
                
            % Assemble the trace data.
            else
                trace_data = trace_dataLR;
            end
            
            % Add the traces.
            traces{i} = cat(1, traces{i}, trace_data);
        end
        
        % Rescale the traces.
        if is_rescale
            traces{i} = traces{i} - min(traces{i}, [], 2);
            traces{i} = traces{i} ./ max(traces{i}, [], 2);
            traces{i} = traces{i} * 0.9 + 0.05;
        end
        
        % Assemble the label.
        exp_str = strrep(num2str(exps_i), '  ', ',');
        neurons{i} = [nameLR ' (' stim_str ', #' exp_str ')'];
        
    % Add the left and/or right neural responses.
    else
        
        % Assemble the stimulus responses.
        for j = 1:length(stims_i)
            
            % Determine the stimulus.
            stim_str(j) = stim_name(stims_i(j));
            switch stims_i(j)
                case 1
                    dataL = neuron_data(neuron_i).left.but;
                    dataR = neuron_data(neuron_i).right.but;
                case 2
                    dataL = neuron_data(neuron_i).left.pent;
                    dataR = neuron_data(neuron_i).right.pent;
                case 3
                    dataL = neuron_data(neuron_i).left.nacl;
                    dataR = neuron_data(neuron_i).right.nacl;
                case 4
                    dataL = neuron_data(neuron_i).left.light;
                    dataR = neuron_data(neuron_i).right.light;
            end
            
            % Assemble the responses.
            if stims_i(j) <  4
                traceL = [dataL.pre.trace(exps_i,:), ...
                    dataL.stim.trace(exps_i,2:(end-1)), ...
                    dataL.post.trace(exps_i,:)];
                traceR = [dataR.pre.trace(exps_i,:), ...
                    dataR.stim.trace(exps_i,2:(end-1)), ...
                    dataR.post.trace(exps_i,:)];
                
            % Assemble the light responses.
            else
                traceL = [dataL.stim.trace(exps_i,1:(end-1)), ...
                    dataL.post.trace(exps_i,:)];
                traceR = [dataR.stim.trace(exps_i,1:(end-1)), ...
                    dataR.post.trace(exps_i,:)];
            end
            
            % Determine the side(s) to assemble.
            trace_length = 2 * length(exps_i);
            if LR == 'L'
                trace_length = length(exps_i);
                traceR = [];
            elseif LR == 'R'
                trace_length = length(exps_i);
                traceL = [];
            end
            
            % Assemble the light trace data.
            if any(stims_i < 4) && stims_i(j) ==  4
                
                % We're plotting NaCl.
                if any(stims_i == 3)
                    trace_data = nan(trace_length, nacl_trace_size);
                    trace_data(:,light_nacl_i) = ...
                        cat(1, traceL, traceR);
                    
                % We're plotting odors.
                else
                    trace_data = nan(trace_length, odor_trace_size);
                    trace_data(:,light_odor_i) = ...
                        cat(1, traceL, traceR);
                end
                
            % Assemble the odor trace data.
            elseif any(stims_i == 3) && stims_i(j) < 3
                trace_data = nan(trace_length, nacl_trace_size);
                trace_data(:,odor_nacl_i) = cat(1, traceL, traceR);
                
            % Assemble the trace data.
            else
                trace_data = cat(1, traceL, traceR);
            end
            
            % Add the traces.
            traces{i} = cat(1, traces{i}, trace_data);
        end
        
        % Rescale the traces.
        if is_rescale
            traces{i} = traces{i} - min(traces{i}, [], 2);
            traces{i} = traces{i} ./ max(traces{i}, [], 2);
            traces{i} = traces{i} * 0.9 + 0.05;
        end
        
        % Assemble the label.
        exp_str = strrep(num2str(exps_i), '  ', ',');
        neurons{i} = [nameLR ' (' stim_str ', #' exp_str ')'];
    end
end

% Determine the maximum trace size.
max_trace = 20;
if any(stims_used == 3)
    max_trace = 60;
elseif any(stims_used < 3)
    max_trace = 30;
end
x_limits = [0,max_trace];

% Setup the figure properties.
figure_size = get(0, 'Screensize');
figure_size(1:2) = round(figure_size(1:2) + figure_size(3:4)./8);
figure_size(3:4) = round(figure_size(3:4) .* 3/4);
trace_background = [0.9, 0.9, 0.9];
trace_background_x = [0, 0, max_trace, max_trace];
plot_label_x = 'Time (s)';
plot_label_y = 'Neural Activity (scaled fluoresence)';

% Setup the font & line properties.
font = 'Arial';
font_size = 10;
plot_line_width = 2;
stim_line_width = 2;

% Create the figure.
fig = figure('NumberTitle', 'off', 'Name', 'Neurons');
set(fig, 'Position', figure_size);
hold on;

% Initialize the stimulus plots.
stim_alpha = 0.3;
stim_trace_x(1,:) = [10,10,20,20]; % butanone
stim_trace_x(2,:) = [10,10,20,20]; % pentanedione
stim_trace_x(3,:) = [20,20,40,40]; % NaCl
stim_trace_x(4,:) = [0,0,10,10]; % light
stim_trace_y = [0, 1, 1, 0];
stim_color(1,:) = [0,1,1]; % cyan
stim_color(2,:) = [1,0,1]; % magenta
stim_color(3,:) = [1,1,0]; % yellow
stim_color(4,:) = [0.5,0.5,1]; % purple

% Initialize the error plots.
error_color = [0,0,0];
error_alpha = 0.25;

% Order the neural traces by name.
neurons = flip(neurons);
traces = flip(traces);
stim_i = flip(stim_i);

% Demarcate the traces.
for i=1:length(neurons)
    if mod(i,2) == 0
        fill(trace_background_x, [i-1, i, i, i-1], trace_background);
    end
end

% Plot the neural traces & their stimuli.
for i=1:length(neurons)
    
    % Plot the stimulus graphic(s).
    stims_i = unique(stim_i{i});
    for j = 1:length(stims_i)
        
        % Offset odor stimuli.
        plot_off = 0;
        if stims_i(j) < 3
             if any(stims_used == 3)
                 plot_off = 10;
             end
             
        % Offset light stimuli.
        elseif stims_i(j) == 4
            if any(stims_used == 3)
                plot_off = 20;
            elseif any(stims_used < 3)
                plot_off = 10;
            end
        end
        
        % Plot the stimulus.
        stims_trace_x = stim_trace_x(stims_i(j),:) + plot_off;
        stims_trace_y = stim_trace_y + i - 1;
        stims_color = stim_color(stims_i(j),:);
        fill(stims_trace_x, stims_trace_y, stims_color, ...
            'FaceAlpha', stim_alpha, 'LineWidth', stim_line_width, ...
            'LineStyle', ':');
    end
    
    % Offset odor plots.
    plot_off = 0;
    if all(stims_i < 3)
        if any(stims_used == 3)
            plot_off = 10;
        end
        
    % Offset light plots.
    elseif all(stims_i == 4)
        if any(stims_used == 3)
            plot_off = 20;
        elseif any(stims_used < 3)
            plot_off = 10;
        end
    end
    
    % Determine the trace time & y-axis offset;
    trace_data = traces{i};
    time_x = (0:(size(trace_data,2) - 1)) / fps + plot_off;
    y_off = i - 1;
    
    % Determine the trace graphic(s).
    if plot_mode > 1
        
        % Compute the mean & error.
        trace_error = nanstd(trace_data, [], 1);
        if plot_mode > 2 % compute the SEM
            trace_error = trace_error / sqrt(size(trace_data,1));
        end
        trace_data = nanmean(trace_data, 1);
        
        % Compute the lower & upper error bounds.
        lower_error = trace_data - trace_error;
        upper_error = trace_data + trace_error;
        
        % Rescale everything uniformly to [0.05,0.95].
        offset = min(lower_error, [], 'all');
        scale = max(upper_error - offset, [], 'all');
        lower_error = (lower_error - offset) / scale;
        lower_error = lower_error * 0.9 + 0.05;
        upper_error = (upper_error - offset) / scale;
        upper_error = upper_error * 0.9 + 0.05;
        trace_data = (trace_data - offset) / scale;
        trace_data = trace_data * 0.9 + 0.05;
        
        % Plot the error
        error_x = [time_x, fliplr(time_x)];
        error_y = [upper_error, fliplr(lower_error)] + y_off;
        fill(error_x, error_y, error_color, 'FaceAlpha', error_alpha, ...
            'LineStyle', 'none');
    
    % Rescale all traces uniformly to [0.05,0.95].
    else
        trace_data = trace_data - min(trace_data, [], 'all');
        trace_data = trace_data / max(trace_data, [], 'all');
        trace_data = trace_data * 0.9 + 0.05;
    end

    % Plot the neural trace.
    plot(time_x, trace_data' + y_off, 'k', 'LineWidth', plot_line_width);
end

% Label the plot.
set(gca, 'FontName', font, 'FontSize', font_size);
xlabel(plot_label_x);
ylabel(plot_label_y);
yticks((1:length(neurons)) - 0.5);
yticklabels(neurons);

% Constrain the axes.
xlim(x_limits);
ylim([0,length(neurons)]);
end
