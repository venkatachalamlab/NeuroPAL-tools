function plotNeuronResponses(exps, neurons, fps, stim_i, names, varargin)
%PLOTNEURONRESPONSES Plot the responses for a set of neurons.
%
%   PLOTNEURONRESPONSES(EXPS, NEURONS, FPS, STIM_I, NAMES)
%   PLOTNEURONRESPONSES(EXPS, NEURONS, FPS, STIM_I, NAMES, IS_RESCALE)
%
%   Input:
%       exps       - the experiments
%       neurons    - an array of neurons structs with the info
%       fps        - frames/second
%       stim_i     - the indices of the stimuli to plot
%                    1 = butanone
%                    2 = pentanedione
%                    3 = NaCl
%                    4 = light
%       names      - the names of the neurons to plot
%       is_rescale - rescale to [0.05,0.95]?
%                    default: false

% REMOVED FUNCTIONALITY
%
%       filter    - filter the traces with a causal bandpass?
%                   default: false
%       passband  - the passband for the filter in fps (low, high, crop)
%                   default: [2,10,4]

% Rescale to [0.05,0.95]?
is_rescale = false;
if ~isempty(varargin) && ~isempty(varargin{1})
    is_rescale = varargin{1};
end

% Initialize the stimulus name.
% Note: we'll throw an error if it's out of range.
stim_names = {'Butanone', 'Pentanedione', 'NaCl', 'Light', };

% Initilaize the experiment names.
exp_names = strrep({exps.name}, '_traces.mat', '');
exp_names = strrep(exp_names, '_', ' ');

% Initialize the neuron names.
neuron_names = {neurons.name};
if ischar(names)
    names = {names};
end

% Are we filtering the traces?
%is_filter = false;
%passband = [2,10,4];
%if ~isempty(varargin)
%    is_filter = varargin{1};
%    
%    % Determine the passband.
%    if length(varargin) > 1
%        passband = varargin{1};
%    end
%end
%passband = round(passband * fps);

% Plot the neural responses.
for i = 1:length(names)
    
    % Plot the neuron's response per stimulus.
    for j = 1:length(stim_i)
        
        % Find the neuron.
        k = find(startsWith(neuron_names, names{i}), 1);
        
        % Determine the stimulus.
        stim_name = stim_names{stim_i(j)};
        
        % Plot the neural response.
        if ~neurons(k).isLR
            
            % Assemble the response.
            switch stim_i(j)
                case 1
                    data = neurons(k).but;
                case 2
                    data = neurons(k).pent;
                case 3
                    data = neurons(k).nacl;
                case 4
                    data = neurons(k).light;
            end
            
            % Filter the data.
            %if is_filter
            %    data = filterData(data, passband);
            %end
            
            % Plot the data.
            fig = figure('NumberTitle', 'off', 'Name', ...
                [neurons(k).name ' responding to ' stim_name]);
            pos = get(fig, 'Position');
            width = 1.5 * pos(3);
            height = 2 * pos(4);
            pos(1) = pos(1) - round(width/3);
            pos(3) = width;
            pos(4) = height;
            set(fig, 'Position', pos);
            plotData(data, exp_names, stim_name, fps, neurons(k).name, ...
                is_rescale);
            
        % Plot the left & right neural responses.
        else
            
            % Assemble the response.
            switch stim_i(j)
                case 1
                    dataL = neurons(k).left.but;
                    dataR = neurons(k).right.but;
                case 2
                    dataL = neurons(k).left.pent;
                    dataR = neurons(k).right.pent;
                case 3
                    dataL = neurons(k).left.nacl;
                    dataR = neurons(k).right.nacl;
                case 4
                    dataL = neurons(k).left.light;
                    dataR = neurons(k).right.light;
            end
            
            % Filter the data.
            %if is_filter
            %    dataL = filterData(dataL, passband);
            %    dataR = filterData(dataR, passband);
            %end
            
            % Plot the data.
            fig = figure('NumberTitle', 'off', 'Name', neurons(k).name);
            pos = get(fig, 'Position');
            width = 3 * pos(3);
            height = 2 * pos(4);
            pos(1) = pos(1) - round(width/3);
            pos(3) = width;
            pos(4) = height;
            set(fig, 'Position', pos);
            subplot(1,2,1);
            plotData(dataL, exp_names, stim_name, fps, ...
                [neurons(k).name 'L'], is_rescale);
            subplot(1,2,2);
            plotData(dataR, exp_names, stim_name, fps, ...
                [neurons(k).name 'R'], is_rescale);
            
        end
    end
end
end

% Filter the data.
%function data = filterData(data, passband)
%data.pre.trace = causalBandpassFilter(data.pre.trace, ...
%    passband(1), passband(2), passband(3));
%data.stim.trace = causalBandpassFilter(data.stim.trace, ...
%    passband(1), passband(2), passband(3));
%data.post.trace = causalBandpassFilter(data.post.trace, ...
%    passband(1), passband(2), passband(3));
%end

% Plot the data.
function plotData(data, exp_names, stim_name, fps, name, is_rescale)

% Initialize the plot info.
font_size = 12;
trace_line = 1.5;
stim_line = 1.5;

% Is their any pre-stimulus data.
is_pre = isfield(data, 'pre');

% Assemble the traces.
if is_pre
    traces = [data.pre.trace, data.stim.trace(:,2:(end-1)), data.post.trace];
else  % light response
    traces = [data.stim.trace(:,2:(end-1)), data.post.trace];
end

% Rescale each trace to [0.05,0.95].
if is_rescale
    traces = traces - min(traces, [] ,2);
    traces = traces ./ max(traces, [] ,2);
    traces = traces * 0.9 + 0.05;

% Rescale all traces uniformly to [0.05,0.95].
else
    traces = traces - min(traces(:));
    traces = traces / max(traces(:));
    traces = traces * 0.9 + 0.05;
end

% Measure the effect sizes.
%mean_dP = nanmean(data.on.diff);
if is_pre
    effects = data.on.diff ./ data.pre.mean;
else % light response
    effects = data.off.diff ./ data.stim.mean;
end
mean_effect = nanmean(effects);
std_effect = nanstd(effects);
z_effects = (effects - mean_effect) ./ std_effect;
s_effect = std_effect / mean_effect;
samples = sum(~isnan(z_effects));

% Scale the effect sizes by smoothness.
dtraces = diff(traces,1,2);
std_dtraces = nanstd(dtraces,0,2)';
zs_effects = z_effects ./ std_dtraces;

% Plot the traces for visibility.
for i = 1:size(traces,1)
    traces(i,:) = traces(i,:) + size(traces,1) - i;
end

% Compute the stimulus start & end.
if is_pre
    stim_start = (size(data.pre.trace,2) - 1) / fps;
    stim_end = (size(data.pre.trace,2) + size(data.stim.trace,2) - 2) / fps;
    stim_height = size(data.pre.trace,1);
else % light response
    stim_start = 0;
    stim_end = (size(data.stim.trace,2) - 1) / fps;
    stim_height = size(data.stim.trace,1);
end

% Plot the response.
plot((0:(size(traces,2) - 1)) / fps, traces, 'LineWidth', trace_line);
hold on;
plot([stim_start,stim_start], [0,stim_height], 'k--', 'LineWidth', stim_line);
plot([stim_end,stim_end], [0,stim_height], 'k--', 'LineWidth', stim_line);

% Update the plot labels.
legend_names = cell(size(exp_names));
for i = 1:length(exp_names)
    legend_names{i} = sprintf('%s ~Z=%0.1f', exp_names{i}, zs_effects(i));
end
neuron_name = [name ' (\DeltaS=' num2str(mean_effect,'%0.2f') ...
    ', \sigma/\mu=' num2str(s_effect,'%0.2f') ...
    ', N=' num2str(samples) ')'];

% Label the plot.
legend(legend_names, 'FontSize', font_size);
title([neuron_name ' responding to ' stim_name]);
xlabel('Time (s)');
end
