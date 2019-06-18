function plotExpTraces(exp_name, varargin)
%PLOTEXPTRACES Plot an experiment's neural traces.
%
%   PLOTEXPTRACES(EXP_NAME)
%   PLOTEXPTRACES(EXP_NAME, IS_COMBINE_LR)
%   PLOTEXPTRACES(EXP_NAME, IS_COMBINE_LR, IS_1_FIG)
%   PLOTEXPTRACES(EXP_NAME, IS_COMBINE_LR, IS_1_FIG, PLOT_MODE)
%   PLOTEXPTRACES(EXP_NAME, IS_COMBINE_LR, IS_1_FIG, PLOT_MODE, IS_SMOOTH)
%   PLOTEXPTRACES(EXP_NAME, IS_COMBINE_LR, IS_1_FIG, PLOT_MODE, IS_SMOOTH,
%                 NEURONS)
%   PLOTEXPTRACES(EXP_NAME, IS_COMBINE_LR, IS_1_FIG, PLOT_MODE, IS_SMOOTH,
%                 NEURONS, FRAMES)
%   PLOTEXPTRACES(EXP_NAME, IS_COMBINE_LR, IS_1_FIG, PLOT_MODE, IS_SMOOTH,
%                 NEURONS, FRAMES, IS_INFO)
%   PLOTEXPTRACES(EXP_NAME, IS_COMBINE_LR, IS_1_FIG, PLOT_MODE, IS_SMOOTH,
%                 NEURONS, FRAMES, IS_INFO, IS_SAVE, SAVE_SUFFIX)
%
%   Inputs:
%   exp_name      = the name of the experiment to plot
%   is_combine_LR = combine left/right neural designations?
%       default: false
%   is_1_fig      = show all neurons on one figure?
%       default: false
%   plot_mode     = how should we plot the neural traces?
%       default: 3
%       1 = plot the neural traces individually
%       2 = plot the neural traces as mean & standard deviation
%       3 = plot the neural traces as mean & SEM
%       4 = plot the neural + filtered traces individually
%       5 = plot the neural + filtered traces as mean & standard deviation
%       6 = plot the neural + filtered traces as mean & SEM
%   is_smooth     = smooth the traces?
%                   (causal bandpass filter)
%       default: false
%   neurons       = the names of the neurons to plot
%       default: [] (all neurons)
%   time          = the time (frames) to plot, a 1-2 length vector:
%                   (start_time, end_time) or (start_time)
%       default: [] (all time/frames)
%   is_info       = show the scale & offset on the plots?
%   is_save       = save the figures to files?
%   save_suffix   = a suffix to append to saved files

% Are we plotting a single experiment?
if ischar(exp_name)
    exp_name = {exp_name};
end

% Are we combining left/right neural designations?
is_combine_LR = false;
if ~isempty(varargin)
    is_combine_LR = varargin{1};
end

% Are we showing all neurons on one figure?
is_1_fig = false;
if length(varargin) > 1
    is_1_fig = varargin{2};
end

% How should we plot the neural traces?
plot_mode = 3;
interp_type = 'linear';
plot_mode_label = [];
if length(varargin) > 2
    plot_mode = varargin{3};
end
switch plot_mode
    case {1,4}
        % do nothing
    case {2,5}
        plot_mode_label = ' Mean +/- Standard Deviation';
    case {3,6}
        plot_mode_label = ' Mean +/- SEM';
    otherwise
        error('Unknown plot mode: %d', plot_mode);
end

% Are we smoothing the traces (Savitzky-Golay filter with a ~1.5s window)?
is_smooth = false;
if length(varargin) > 3 && ~isempty(varargin{4})
    is_smooth = varargin{4};
    
    % Determine the passband.
    if islogical(is_smooth) && is_smooth
        %passband = [7, 70, 5];
        %passband = [2, 125, 2];
        passband = [2, 240, 1];
        
    % Set the passband.
    elseif is_smooth
        passband = is_smooth;
        is_smooth = true;
    end
end

% Which neurons are we plotting?
plot_neurons = [];
if length(varargin) > 4
    plot_neurons = varargin{5};
end

% Which frames are we plotting?
plot_times = [];
if length(varargin) > 5
    plot_times = varargin{6};
end

% Are we showing the scale & offset on the plots?
is_plot_info = false;
if length(varargin) > 6 && ~isempty(varargin{7})
    is_plot_info = varargin{7};
end

% Are we saving the figures to files?
is_save = false;
if length(varargin) > 7
    is_save = varargin{8};
end

% Are we appending a suffix to saved figure files?
save_suffix = [];
if length(varargin) > 8
    save_suffix = varargin{9};
end

%% Initialize ganglia, neurons, and stimuli properties.
% Initialize the ganglia & neuron names.
ganglion = [];
if ~is_1_fig
    if is_combine_LR
        load('ganglia.mat');
    else
        load('ganglia_LR.mat');
        ganglion = ganglion_LR;
    end
end

% Initialize the neuron colors.
load('neuron_colors.mat');
neuron_colors = neurons;
neuron_color_list = {neuron_colors.name};
black = [0,0,0];
darken = 0.3*[1,1,1];
trace_error_alpha = 0.3;

% Initialize the stimuli.
load('stimuli.mat');
stim_list = stimulus;
stim_list_keywords = {stim_list.keyword};
stimulus_offset_threshold = 1;
stimulus_alpha = 0.2;

%% Organize the experiments.
% Organize the experiments.
neuron = [];
stimulus = [];
max_fps = 0;
min_trace_time = 0;
max_trace_time = 0;
behaviors = [];
behaviors_offset = [];
behaviors_fps = [];
for i = 1:length(exp_name)
    
    % Open the experiment.
    exp_data = open(exp_name{i});
    
    % Organize the frame times.
    normalize_fps = true; % normalize the fps to a standard rate
    fps = exp_data.fps;
    if isfield(exp_data, 'times')
        exp_times = exp_data.times.times;
        
    % We don't have the time per each frame, interpolate it.
    else
        warning('Cannot find "times" in %s!', exp_name{i});
        
        % Compute the frame times.
        normalize_fps = false;
        num_frames = size(exp_data.gcamp, 2);
        exp_times = linspace(1 / fps, num_frames / fps, num_frames);
        
        % Compute the stimulus times.
        for j = 1:size(exp_data.stimulus,1)
            exp_data.stimulus{j,2} = exp_data.stimulus{j,2} / fps;
            exp_data.stimulus{j,3} = exp_data.stimulus{j,3} / fps;
        end
    end
    
    % Organize the stimuli.
    exp_stim = exp_data.stimulus;
    for j=1:length(exp_stim)
        
        % Which of the stimuli is being applied?
        stim_name = lower(exp_stim{j,1});
        iStim = find(cellfun(@(x) contains(stim_name,x), stim_list_keywords));
        if isempty(iStim)
            error('Unknown stimulus: %s', exp_stim{j,1});
        end
        
        % Organize the stimulus.
        stimulus(i).num(j) = iStim;
        stimulus(i).on(j) = exp_stim{j,2};
        stimulus(i).off(j) = exp_stim{j,3};
        stimulus(i).mid(j) = (stimulus(i).on(j) + stimulus(i).off(j)) / 2;
    end
    
    % Check the previous stimuli.
    stim_offset = 0;
    stimulus(i).offset = 0;
    if i > 1
        % Check the order of the stimuli.
        if any([stimulus(i).num] ~= [stimulus(1).num])
            error('Incompatible order of stimuli between "%s" and "%s".', ...
                exp_name{1}, exp_name{i});
        end
        
        % Determine the offset of the stimuli.
        stim1_mid = [stimulus(1).mid];
        stim_mid = [stimulus(i).mid];
        stim_offset = mean(stim1_mid - stim_mid);
        stim_off_diffs = ([stimulus(i).off] + stim_offset) - [stimulus(1).off];
        stim_on_diffs = ([stimulus(i).on] + stim_offset) - [stimulus(1).on];
        if any(stim_off_diffs > stimulus_offset_threshold) || ...
                any(stim_on_diffs > stimulus_offset_threshold)
            error('Incompatible stimuli intervals between "%s" and "%s".', ...
                exp_name{1}, exp_name{i});
        end
        
        % Set the minimum trace & stimulus offset.
        min_trace_time = min(stim_offset, min_trace_time);
        stimulus(i).offset = stim_offset;
    end
    
    % Fix the neuron names.
    exp_neurons = exp_data.neuron_names;
    if iscategorical(exp_neurons)
        exp_neurons = cellstr(exp_neurons);
    end
    
    % Organize the traces & positions.
    exp_traces = exp_data.gcamp;
    exp_pos = exp_data.positions;

    % Which neurons are we plotting?
    if ~isempty(plot_neurons)
        
        % Organize the plotted neurons.
        if ~iscell(plot_neurons)
            plot_neurons = {plot_neurons};
        end
        
        % Determine which neurons to plot.
        keep_neurons = false(length(exp_neurons),1);
        for j = 1:length(plot_neurons)
            keep = startsWith(exp_neurons, plot_neurons{j});
            keep_neurons(keep) = true;
        end
        
        % Keep the plotted neurons.
        exp_neurons = exp_neurons(keep_neurons);
        exp_traces = exp_traces(keep_neurons,:);
        exp_pos = exp_pos(keep_neurons,:,:);
    end
    
    % Sort the neurons.
    [exp_neurons, sort_i] = sort(exp_neurons);
    exp_traces = exp_traces(sort_i,:);
    exp_pos = exp_pos(sort_i,:,:);
    
    % Normalize the traces & positions to a standard frame rate.
    % Note: the camera's fps is variable & needs to be interpolated to a
    % standard rate.
    %fps = 4;
    if normalize_fps
        [exp_traces, new_times] = normalizeTimedData(exp_traces, exp_times, fps);
        [exp_pos, ~] = normalizeTimedData(exp_pos, exp_times, fps);
        exp_times = new_times;
    end
    
    % Which frames (time) are we plotting?
    if ~isempty(plot_times)

        % Translate the time to frames.
        plot_times(1) = round(plot_times(1) * fps) + 1;
        if length(plot_times) > 1
            plot_times(2) = round(plot_times(2) * fps) + 1;
        else
            plot_times(2) = size(exp_traces,2);
        end
        plot_frames = plot_times(1):plot_times(2);
        
        % Keep the plotted frames.
        exp_traces = exp_traces(:,plot_frames);
        exp_pos = exp_pos(:,:,plot_frames);
        exp_times = exp_times(plot_frames);
        
        % Correct the stimulus times.
        start_off = plot_frames(1) / fps;
        stimulus(i).on = stimulus(i).on - start_off;
        stimulus(i).off = stimulus(i).off - start_off;
        stimulus(i).mid = stimulus(i).mid - start_off;
    end
    
    % Compute the maxes for fps & trace times.
    %spf = nanmedian(diff(exp_times));
    %fps = 1 / spf;
    max_fps = nanmax(fps, max_fps);
    max_trace_length = size(exp_traces,2) - 1;
    max_trace_time = nanmax(max_trace_time, max_trace_length/fps + stim_offset);
    
    % Filter the neural traces.
    exp_ftraces = [];
    if plot_mode > 3
        [exp_ftraces, exp_fscales, exp_foffsets, exp_fbaselines, ~] = ...
            filterTraces(exp_traces, fps, exp_times);
    end
    
    % Clean up the neural traces.
    if is_smooth
        passband_fps = round(passband * fps);
        [exp_traces, exp_scales, exp_offsets] = ...
            cleanTraces(exp_traces, fps, 10, 2, [], 'causal', passband_fps);
        %cleanTraces(exp_traces, fps, [], [], [], 'sgolay', round(1.5*fps));
    else
        [exp_traces, exp_scales, exp_offsets] = ...
            cleanTraces(exp_traces, fps);
            %cleanTraces(exp_traces, fps, 10, 3);
    end
    
    % Scale & offset the filtered neural traces to match the real ones.
    if ~isempty(exp_ftraces)
        fscales = exp_fscales ./ exp_scales;
        foffsets = (exp_foffsets - exp_offsets) ./ exp_scales;
        exp_ftraces = exp_ftraces .* fscales + foffsets;
        exp_fbaselines = exp_fbaselines .* fscales + foffsets;
    end
    
    % Remove neural positions with no trace data.
    [iNeuron, iFrame] = find(isnan(exp_traces));
    for j = 1:length(iNeuron)
        exp_pos(iNeuron(j),:,iFrame(j)) = nan;
    end
    
    % Remove unlabeled neurons.
    unlabeled = startsWith(exp_neurons, 'Q') | ...
        startsWith(exp_neurons, 'STITCH');
    exp_neurons(unlabeled) = [];
    exp_traces(unlabeled,:) = [];
    exp_scales(unlabeled,:) = [];
    exp_offsets(unlabeled,:) = [];
    exp_pos(unlabeled,:,:) = [];
    if ~isempty(exp_ftraces)
        exp_ftraces(unlabeled,:) = [];
        exp_fbaselines(unlabeled,:) = [];
    end
    
    % Warn the user about missing frames.
    missing_frames = find(all(isnan(exp_traces),1));
    if ~isempty(missing_frames)
        missing_str = ['Missing frames: ' num2str(missing_frames(1))];
        for j=2:length(missing_frames)
            missing_str = [missing_str ', ' num2str(missing_frames(j))];
        end
        warning(missing_str);
    end
    
    % Warn the user about missing neurons.
    missing_neurons = find(all(isnan(exp_traces),2));
    if ~isempty(missing_neurons)
        missing_str = ['Missing neurons: ' exp_neurons{missing_neurons(1)}];
        for j=2:length(missing_neurons)
            missing_str = [missing_str ', ' exp_neurons{missing_neurons(j)}];
        end
        warning(missing_str);
    end

    % Organize the neurons.
    neuron_offset = length(neuron);
    for j = 1:length(exp_neurons)
        iNeuron = j + neuron_offset;
        if is_combine_LR
            neuron(iNeuron).name = stripNeuronLR(exp_neurons{j});
        else
            neuron(iNeuron).name = exp_neurons{j};
        end
        neuron(iNeuron).RGB = [];
        neuron(iNeuron).time{1} = exp_times;
        neuron(iNeuron).trace{1} = exp_traces(j,:);
        neuron(iNeuron).scale(1) = exp_scales(j);
        neuron(iNeuron).trace_offset(1) = exp_offsets(j);
        neuron(iNeuron).stim_offset(1) = stim_offset;
        neuron(iNeuron).fps(1) = fps;
        neuron(iNeuron).stimuli = i;
        if ~isempty(exp_ftraces)
            neuron(iNeuron).ftrace{1} = exp_ftraces(j,:);
            neuron(iNeuron).fbaseline(1) = exp_fbaselines(j,:);
        end
    end
    
    % Combine the neurons.
    neuron_names = {neuron.name};
    [~, iUniqueNeurons, iNeurons] = unique(neuron_names);
    if length(iUniqueNeurons) ~= length(neuron_names)
        neuron2 = [];
        for j = 1:length(iUniqueNeurons)
            neuron2(j).name = neuron(iUniqueNeurons(j)).name;
            neuron2(j).RGB = [];
            iNeuron = find(j == iNeurons);
            neuron2(j).time = [neuron(iNeuron).time];
            neuron2(j).trace = [neuron(iNeuron).trace];
            neuron2(j).scale = [neuron(iNeuron).scale];
            neuron2(j).trace_offset = [neuron(iNeuron).trace_offset];
            neuron2(j).stim_offset = [neuron(iNeuron).stim_offset];
            neuron2(j).fps = [neuron(iNeuron).fps];
            neuron2(j).stimuli = [neuron(iNeuron).stimuli];
            if plot_mode > 3
                neuron2(j).ftrace = [neuron(iNeuron).ftrace];
                neuron2(j).fbaseline = [neuron(iNeuron).fbaseline];
            end
        end
        neuron = neuron2;
    end
    
    % Organize the neuron colors.
    for j=1:length(neuron)
        
        % Get the neuron's color.
        neuron(j).RGB = black;
        name = stripNeuronLR(neuron(j).name, true);
        iColor = find(strcmp(neuron_color_list, name));
        if ~isempty(iColor)
            color = neuron_colors(iColor).RGB;
            color = color - darken;
            neuron(j).RGB = max(color, black);
        end
    end
    
    % Compute the behavior.
    behavior = computeBehavior(exp_neurons, exp_traces, exp_pos, fps, false);
    if i > 1 && length(behavior) ~= length(behaviors)
        error('Incompatible behavior computations between "%s" and "%s".', ...
            exp_name{1}, exp_name{i});
    end
    
    % Setup the first experiment's behaviors.
    if i == 1
        behaviors = behavior;
        behaviors_offset(i) = stim_offset;
        behaviors_fps(i) = fps;
        behaviors_time{i} = exp_times;
        for j=1:length(behaviors)
            behaviors(j).data = {behavior(j).data};
        end
        
    % Add the experiment's behaviors.
    else
        for j=1:length(behaviors)
            if ~strcmp(behaviors(j).name, behavior(j).name)
                error('Incompatible behavior computations between "%s" and "%s".', ...
                        exp_name{1}, exp_name{i});
            end
            behaviors(j).data{i} = behavior(j).data;
            behaviors(j).center(i) = behavior(j).center;
            behaviors_offset(i) = stim_offset;
            behaviors_fps(i) = fps;
            behaviors_time{i} = exp_times;
        end
    end
end

%% Initialize the behavioral figure properties.
% Setup the figure properties.
figure_size = get(0, 'Screensize');
figure_size(1:2) = round(figure_size(1:2) + figure_size(3:4)./8);
figure_size(3:4) = round(figure_size(3:4) .* 3/4);
trace_time_end = max_trace_time - min_trace_time;
plot_label_x = 'Time (s)';

%% Plot the behavior.
% Setup the figure.
font = 'Arial';
font_size = 10;
line_width = 1;
fig_name = 'Worm Behavior';
if is_save
    fig = figure('Visible', 'off', 'NumberTitle', 'off', 'Name', fig_name);
    font_size = 6;
    line_width = 0.5;
else
    fig = figure('NumberTitle', 'off', 'Name', fig_name);
    set(gcf, 'Position', figure_size);
end
hold on;

% Plot the behaviors.
num_behaviors = length(behaviors);
for i=1:num_behaviors

    % Plot the behavior.
    subplot(num_behaviors,1,i);
    hold on;
    plotTrace([], behaviors_time, behaviors(i).data, [], [], ...
        behaviors_offset, behaviors_fps, max_fps, ...
        min_trace_time, max_trace_time, ...
        black, trace_error_alpha, line_width, ...
        false, plot_mode, interp_type);
    
    % Label the plot.
    title_text = [behaviors(i).name ': ' behaviors(i).ylabel];
    if i == 1 % add a newline
        title_text = {title_text, ''};
    end
    title(title_text, 'FontName', font, 'FontSize', font_size);
    xlabel(plot_label_x, 'FontName', font, 'FontSize', font_size);
    
    % Compact the plot.
    axis tight;
    xlim([0, trace_time_end]);
    
    % Compactly label the y-axis.
    xlimits = xlim;
    ylimits = ylim;
    mean_y = mean(ylimits);
    line(xlimits, [mean_y, mean_y], 'LineWidth', line_width, ...
        'LineStyle', ':');
    yticks(mean_y);
    
    % Set the tick fonts.
    set(gca, 'FontName', font, 'FontSize', font_size);
    
    % Plot the stimuli.
    add_title = false;
    if i == 1 % title the stimuli
        add_title = true;
    end
    plotStimuli(stimulus, stim_list, ylimits, stimulus_alpha, ...
        min_trace_time, add_title, font, font_size);
end

% Save the figure.
if is_save

    % Name the file.
    behavior_str = 'behavior';
    if length(exp_name) == 1
        filename = strrep(exp_name{1}, '.mat', behavior_str);
    else
        filename = behavior_str;
    end
    
    % Save the file.
    disp(['Printing: ' filename]);
    fig.Renderer = 'Painters';
    orient(fig, 'portrait');
    print(fig, '-fillpage', '-dpdf', filename);
end

%% Initialize the brain activtiy figure properties.
% Setup the figure properties.
figure_size = get(0, 'Screensize');
figure_size(1:2) = round(figure_size(1:2) + figure_size(3:4)./8);
figure_size(3:4) = round(figure_size(3:4) .* 3/4);
trace_background = [0.9, 0.9, 0.9];
trace_background_x = [0, 0, trace_time_end, trace_time_end];
plot_label_x = 'Time (s)';
plot_label_y = ['Neural Activity' plot_mode_label ' (\DeltaF/F_0)'];

% Setup the font & line properties.
font = 'Arial';
font_size = 10;
line_width = 1;
if is_save
    font_size = 6;
    line_width = 0.5;
end

%% Plot the brain activity.
% Plot the traces for each ganglia.
num_plots = 1;
if ~is_1_fig
    num_plots = length(ganglion);
end
for i=1:num_plots
    
    % Organize the gangliar data.
    iNeurons = 1:length(neuron);
    if ~is_1_fig
        [~, ~, iNeurons] = intersect(ganglion(i).neurons, ...
            {neuron.name}, 'stable');
    end
    
    % Plot the neural traces.
    if ~isempty(iNeurons)
        
        % Setup the figure.
        fig_name = 'All Neurons';
        if ~is_1_fig
            fig_name = ganglion(i).name;
        end
        if is_save
            fig = figure('Visible', 'off', 'NumberTitle', 'off', ...
                'Name', fig_name);
        else
            fig = figure('NumberTitle', 'off', 'Name', fig_name);
            set(fig, 'Position', figure_size);
        end
        hold on;
        
        % Order the neural traces by name.
        iNeurons = flip(iNeurons);
        
        % Demarcate the traces.
        for j=1:length(iNeurons)
            if mod(j,2) == 0
                fill(trace_background_x, [j-1, j, j, j-1], ...
                    trace_background);
            end
        end
        
        % Plot the neural traces & create their labels.
        neuron_label = cell(length(iNeurons),1);
        for j=1:length(iNeurons)
            
            % Plot the neural trace.
            plotted_neuron = neuron(iNeurons(j));
            ftrace = [];
            fbaseline = [];
            if plot_mode > 3
                ftrace = plotted_neuron.ftrace;
                fbaseline = plotted_neuron.fbaseline;
            end
            plotTrace(j, plotted_neuron.time, plotted_neuron.trace, ...
                ftrace, fbaseline, plotted_neuron.stim_offset, ...
                plotted_neuron.fps, max_fps, ...
                min_trace_time, max_trace_time, ...
                plotted_neuron.RGB, trace_error_alpha, line_width, ...
                true, plot_mode, interp_type);
            
            % Create the neuron label.
            if is_plot_info
                neuron_label{iNeurons(j)} = ...
                    sprintf('%s (\\lambda=%d, F_{0}=%d)', ...
                    neuron(iNeurons(j)).name, ...
                    round(nanmean(neuron(iNeurons(j)).scale)), ...
                    round(nanmean(neuron(iNeurons(j)).trace_offset)));
            else
                neuron_label{iNeurons(j)} = neuron(iNeurons(j)).name;
            end
        end
        
        % Label the plot.
        if is_save
            title({fig_name, ' '});
        end
        xlabel(plot_label_x);
        ylabel(plot_label_y);
        yticks((1:length(iNeurons)) - 0.5);
        %yticklabels({neuron(iNeurons).name});
        yticklabels(neuron_label(iNeurons));
        
        % Plot the stimuli.
        plotStimuli(stimulus, stim_list, [0 length(iNeurons)], ...
            stimulus_alpha, min_trace_time, true, font, font_size);
        
        % Compact the plot.
        xlim([0, trace_time_end]);
        ylim([0, length(iNeurons)]);
        
        % Set the tick fonts.
        set(gca, 'FontName', font, 'FontSize', font_size);
    
        % Save the figure.
        if is_save
            
            % Name the file.
            ganglia_str = strrep(lower(fig_name), ' ', '_');
            ganglia_str = ['_' ganglia_str save_suffix];
            if length(exp_name) == 1
                filename = strrep(exp_name{1}, '.mat', ganglia_str);
            else
                filename = ganglia_str;
            end
            
            % Save the file.
            disp(['Printing: ' filename]);
            fig.Renderer = 'Painters';
            orient(fig, 'portrait');
            print(fig, '-fillpage', '-dpdf', filename);
        end
    end
end
end

%% Plot behavioral or neural traces.
% Note: the video frames have minimally variable timing.
% Ideally, we should plot using the real time in "times".
function plotTrace(plot_offset, times, traces, ftraces, fbaselines, ...
    time_offsets, fpss, max_fps, min_trace_time, max_trace_time, ...
    color, alpha, line_width, is_rescale, mode, interp_type)

% Initialize the plot info.
fline_style = '-';
fbline_style = ':';
fline_width = 2 * line_width;
fbline_width = 2 * line_width;

% Check the plot mode.
if mode > 3 % plotting filtered traces
    mode = mode - 3;
end

% Check the plot offset.
if isempty(plot_offset)
    plot_offset = 1;
end

% Plot the traces individually.
black = [0,0,0];
num_traces = size(traces,2);
if num_traces == 1 || mode == 1
    line_alpha = 1 / num_traces;
    for i=1:num_traces
        
        % Plot the trace.
        fps = fpss(i);
        trace_y = traces{i};
        trace_x = ((1:length(trace_y)) - 1) ./ fps + time_offsets(i) ...
            - min_trace_time;
        plot(trace_x, trace_y + plot_offset - 1, ...
            'Color', [color line_alpha], 'LineWidth', line_width);
        
        % Plot the filtered trace.
        if ~isempty(ftraces)
            
            % Stay in bounds.
            ftrace_y = ftraces{i};
            fbaseline_y = fbaselines(i);
            min_trace_y = nanmin(trace_y);
            max_trace_y = nanmax(trace_y);
            ftrace_y(ftrace_y < min_trace_y) = min_trace_y;
            ftrace_y(ftrace_y > max_trace_y) = max_trace_y;
            if fbaseline_y < min_trace_y
                fbaseline_y = min_trace_y;
            elseif fbaseline_y > max_trace_y
                fbaseline_y = max_trace_y;
            end
            
            % Plot the filtered trace.
            plot(trace_x, ftrace_y + plot_offset - 1, ...
                'Color', [black line_alpha], 'LineWidth', fline_width, ...
                'LineStyle', fline_style);
            
            % Plot the filtered trace baseline.
            fbaseline_y = zeros(size(trace_x)) + fbaseline_y;
            plot(trace_x, fbaseline_y + plot_offset - 1, ...
                'Color', [black line_alpha], ...
                'LineWidth', fbline_width, 'LineStyle', fbline_style);
        end
    end
    
% Plot the mean trace.
else
    
    % Align the traces.
    frames = round((max_trace_time - min_trace_time) * max_fps);
    trace_x = linspace(0, max_trace_time - min_trace_time, frames);
    combo_traces = nan(num_traces, frames);
    combo_ftraces = [];
    if ~isempty(ftraces)
        combo_ftraces = nan(num_traces, frames);
    end
    for i = 1:num_traces
        
        % Align the traces.
        fps = fpss(i);
        trace_y = traces{i};
        old_trace_x = ((1:length(trace_y)) - 1) ./ fps + time_offsets(i) ...
            - min_trace_time;
        combo_traces(i,:) = interp1(old_trace_x, trace_y, trace_x, ...
            interp_type);
        
        % Align the filtered traces.
        if ~isempty(combo_ftraces)
            combo_ftraces(i,:) = interp1(old_trace_x, ftraces{i}, ...
                trace_x, interp_type);
        end
    end
    
    % Compute mean trace, its error, & the mean baseline.
    mean_trace = nanmean(combo_traces, 1);
    error_trace = nanstd(combo_traces, 0, 1);
    error_nans = isnan(error_trace);
    mean_ftrace = [];
    if ~isempty(combo_ftraces)
        mean_ftrace = nanmean(combo_ftraces, 1);
        mean_fbaseline = nanmean(fbaselines);
    end
    
    % Are we showing the SEM?
    if mode == 3
        error_trace = error_trace ./ sqrt(num_traces);
    end
    
    % Compute the lower & upper traces.
    lower_trace = mean_trace - error_trace;
    upper_trace = mean_trace + error_trace;
    
    % Get rid of nans.
    lower_trace = lower_trace(~error_nans);
    upper_trace = upper_trace(~error_nans);
    lu_trace_x = trace_x(~error_nans);
    
    % Re-scale everything to [0,1].
    if is_rescale
        
        % Set the minimum to 0.
        min_traces = nanmin(nanmin(mean_trace), nanmin(lower_trace));
        if ~isnan(min_traces)
            mean_trace = mean_trace - min_traces;
            lower_trace = lower_trace - min_traces;
            upper_trace = upper_trace - min_traces;
            if ~isempty(mean_ftrace)
                mean_ftrace = mean_ftrace - min_traces;
                mean_fbaseline = mean_fbaseline - min_traces;
            end
        end
        
        % Set the maximum to 1.
        max_traces = nanmax(nanmax(mean_trace), nanmax(upper_trace));
        if ~isnan(max_traces)
            mean_trace = mean_trace ./ max_traces;
            lower_trace = lower_trace ./ max_traces;
            upper_trace = upper_trace ./ max_traces;
            if ~isempty(mean_ftrace)
                mean_ftrace = mean_ftrace ./ max_traces;
                mean_fbaseline = mean_fbaseline ./ max_traces;
            end
        end
    end
    
    % Plot the mean trace.
    fill_x = [lu_trace_x, fliplr(lu_trace_x)];
    fill_y = [upper_trace, fliplr(lower_trace)] + plot_offset - 1;
    fill(fill_x, fill_y, color, 'FaceAlpha', alpha, 'LineStyle', 'none');
    plot(trace_x, mean_trace + plot_offset - 1, 'Color', black, ...
        'LineWidth', line_width);
    
    % Plot the mean filtered trace.
    if ~isempty(mean_ftrace)
        
        % Stay in bounds.
        min_mean_trace = nanmin(mean_trace);
        max_mean_trace = nanmax(mean_trace);
        mean_ftrace(mean_ftrace < min_mean_trace) = min_mean_trace;
        mean_ftrace(mean_ftrace > max_mean_trace) = max_mean_trace;
        if mean_fbaseline < min_mean_trace
            mean_fbaseline = min_mean_trace;
        elseif mean_fbaseline > max_mean_trace
            mean_fbaseline = max_mean_trace;
        end
            
        % Plot the filtered trace.
        plot(trace_x, mean_ftrace + plot_offset - 1, 'Color', black, ...
            'LineWidth', fline_width, 'LineStyle', fline_style);
        
        % Plot the filtered baseline.
        mean_fbaseline = zeros(size(trace_x)) + mean_fbaseline;
        plot(trace_x, mean_fbaseline + plot_offset - 1, 'Color', black, ...
            'LineWidth', fbline_width, 'LineStyle', fbline_style);
    end
end
end

%% Plot stimuli.
function plotStimuli(stimulus, stim_list, ylimits, alpha, ...
    min_trace_time, is_show_name, font, font_size)

% Plot the stimuli.
num_stims = length(stimulus);
for i=1:num_stims
    stim = stimulus(i);
    for j=1:length(stim.num)
        iStim = stim.num(j);
        
        % Plot the stimulus.
        stim_trace_x = [stim.on(j), stim.on(j), stim.off(j), stim.off(j)]  ...
            + stim.offset - min_trace_time;
        stim_trace_y = [ylimits(1), ylimits(2), ylimits(2), ylimits(1)];
        fill(stim_trace_x, stim_trace_y, stim_list(iStim).RGB, ...
            'FaceAlpha', alpha / num_stims, ...
            'LineStyle', ':');
        
        % Label the stimulus.
        if is_show_name && i == num_stims
            stim_mid = stim.mid(j) + stim.offset - min_trace_time;
            text(stim_mid, ylimits(2), stim_list(iStim).name, ...
                'FontName', font, 'FontSize', font_size, ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'bottom');
        end
    end
end
end
