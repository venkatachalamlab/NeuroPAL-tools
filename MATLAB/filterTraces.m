function [traces, scales, offsets, baselines, removed] = ...
    filterTraces(traces, fps, times, varargin)
%FILTERTRACES Filter neural traces to improve the signal/noise, remove
%             inactive neurons, and scale the traces to [0,1].
%
%   [TRACES, SCALES, OFFSETS, BASELINES, REMOVED] =
%       FILTERTRACES(TRACES, FPS, TIMES)
%   [TRACES, SCALES, OFFSETS, BASELINES, REMOVED] =
%       FILTERTRACES(TRACES, FPS, TIMES, TIME_THRESHOLD)
%   [TRACES, SCALES, OFFSETS, BASELINES, REMOVED] =
%       FILTERTRACES(TRACES, FPS, TIMES, TIME_THRESHOLD, MIN_THRESHOLD)
%   [TRACES, SCALES, OFFSETS, BASELINES, REMOVED] =
%       FILTERTRACES(TRACES, FPS, TIMES, TIME_THRESHOLD, MIN_THRESHOLD,
%                    SNR_THRESHOLD)
%   [TRACES, SCALES, OFFSETS, BASELINES, REMOVED] =
%       FILTERTRACES(TRACES, FPS, TIMES, TIME_THRESHOLD, MIN_THRESHOLD,
%                    SNR_THRESHOLD, SIGMA_THRESHOLD)
%   [TRACES, SCALES, OFFSETS, BASELINES, REMOVED] =
%       FILTERTRACES(TRACES, FPS, TIMES, TIME_THRESHOLD, MIN_THRESHOLD,
%                    SNR_THRESHOLD, SIGMA_THRESHOLD, PASSBAND_FREQ,
%                    STOPBAND_FREQ, FILTER_ORDER)
%
%   Inputs:
%   traces        - the neural traces to clean
%   fps           - the sampling rate in frames/second
%   time_thresh   - the minimum time samples to detect signal, in seconds
%                   default: 90 (1.5 minutes)
%   min_thresh    - the minimum fluoresence intensity to detect signal
%                   default: 8 (min=0, max=255)
%   snr_thresh    - the minimum signal/noise defined as:
%                   std(signal)/std(diff(signal))
%                   default: 1.1
%   sigma_thresh  - remove outliers wherein a single frame is sigma
%                   standard deviations from the mean
%                   default: 10
%   passband_freq - the passband frequency for filtering
%                   default: 1/5
%   stobband_freq - the stopband frequency for filtering
%                   default: 1/2
%   filter_order  - the filter order
%                   default: 12
%   Outputs:
%   traces        - the filtered neural traces
%   scales        - the scales used to normalize to unit traces
%   offsets       - the offsets used to normalize to unit traces
%   baselines     - the baselines (modes) for the filtered neural traces
%   removed       - which neural traces were removed?

%% Note: the video frames have minimally variable timing.
%% Ideally, we should plot using the real time in "times".

% Is there enough data?
time_threshold = 90; % 1.5m of data
if ~isempty(varargin) && ~isempty(varargin{1})
    time_threshold = varargin{1};
end

% Is their enough fluorescence to detect signal?
min_threshold = 8;
if length(varargin) > 1 && ~isempty(varargin{2})
    min_threshold = varargin{2};
end

% Is their enough variance in the fluorescence to detect signal changes?
snr_threshold = 1.1;
if length(varargin) > 2 && ~isempty(varargin{3})
    snr_threshold = varargin{3};
end

% Are we removing extreme outliers?
sigma_threshold = 10;
if length(varargin) > 3 && ~isempty(varargin{4})
    sigma_threshold = varargin{4};
end

% How should we filter the traces?
passband_freq = 1/5;
stopband_freq = 1/2;
filter_order = 12;
if length(varargin) > 6 && ~isempty(varargin{5}) ...
        && ~isempty(varargin{6}) && ~isempty(varargin{7})
    passband_freq = varargin{5};
    stopband_freq = varargin{6};
    filter_order = varargin{7};
end

% Clean up the data by removing the first & last frames (outliers).
traces(:,1) = nan;
traces(:,end) = nan;

% Remove traces with too little data.
nans_threshold = min(ceil(time_threshold * fps), size(traces,2));
num_data = sum(~isnan(traces),2);
remove_nans = num_data < nans_threshold;
traces(remove_nans,:) = nan;

% Dampen extreme outliers: a single frame wherein the signal change is
% larger than sigma_threshold standard deviations from the mean.
if ~isempty(sigma_threshold) && abs(sigma_threshold) > 0
    
    % Compute the thresholds for ouliers.
    extreme_threshold = sigma_threshold * nanstd(traces, 0, 2) + ...
        nanmean(traces, 2);
    diff_traces = diff(traces, 1, 2);
    
    % Find the outliers.
    extreme_max_traces = diff_traces > extreme_threshold;
    extreme_min_traces = diff_traces < -extreme_threshold;
    [iExtreme_neurons, iExtreme_frames] = find(...
        (extreme_max_traces(:,1:(end-1)) & extreme_min_traces(:,2:end)) ...
        | (extreme_min_traces(:,1:(end-1)) & extreme_max_traces(:,2:end)));
    iExtreme_frames = iExtreme_frames + 1;
    
    % Remove the outliers.
    for i = 1:length(iExtreme_neurons)
        traces(iExtreme_neurons(i),iExtreme_frames(i)) = nan;
    end
end

% Remove weak traces.
max_traces = nanmax(traces, [], 2);
remove_mins = max_traces < min_threshold;
traces(remove_mins,:) = nan;

% Remove traces with weak signal/noise.
std_traces = nanstd(traces, 0, 2);
std_diff_traces = nanstd(diff(traces, 1, 2), 0, 2);
snr_ratio = std_traces ./ std_diff_traces;
remove_snrs = snr_ratio < snr_threshold;
traces(remove_snrs,:) = nan;

% Which neurons were removed?
removed = remove_nans | remove_mins | remove_snrs;

% Detrend the traces to remove bleaching.
keepI = find(~removed);
detrend_offsets = nan(size(traces,1),1);
x = 1:size(traces,2);
for i = 1:length(keepI)
    yi = traces(keepI(i),:); % remove nan traces
    yi_data = ~isnan(yi); % which part of the trace is real data?
    
    % Compute individual bleach curves, per neuron, & detrend.
    fi = fit(x(yi_data)', yi(yi_data)', 'exp1');
    if fi.b < 0 % bleach curves must decay
        
        % Detrend the trace.
        fi_yi = fi.a * exp(fi.b * x);
        traces(keepI(i),:) = yi - fi_yi;
        detrend_offsets(keepI(i)) = fi.a;
        
    % The trace was not detrended.
    else
        detrend_offsets(keepI(i)) = 0;
    end
end

% Filter out high-frequency noise.
interp_method = 'pchip'; %'linear';
df = designfilt('lowpassfir', ...
    'PassbandFrequency', passband_freq,...
    'StopbandFrequency',stopband_freq, ...
    'FilterOrder',filter_order, ...
    'SampleRate', fps, 'DesignMethod', 'ls');
edge_frames = ceil(fps / (2 * stopband_freq)); %filter_order/2;
start_removeI = x(1:(edge_frames + 1));
end_removeI = x((end - edge_frames):end);
for i = 1:length(keepI)
    yi = traces(keepI(i),:); % remove nan traces
    yi_data = ~isnan(yi); % which part of the trace is real data?
    
    % Filter out high-frequency noise.
    yi = fillmissing(yi, interp_method); % interpolate missing data
    yi = filtfilt(df, yi); % filter the trace & ensure 0-phase
    
    % Store the filtered trace.
    traces(keepI(i),yi_data) = yi(yi_data);
end
traces(:,start_removeI) = nan;
traces(:,end_removeI) = nan;

% Compute the trace baselines (the mode of their probability distribution).
kernel = 'epanechnikov'; % the theoretically optimal kernel
baselines = nan(size(traces,1),1);
for i = 1:length(keepI)
    yi = traces(keepI(i),:); % remove nan traces

    % Estimate the probability density for the trace and compute its mode.
    [yi_probability, yi_value] = ...
        ksdensity(yi, 'Kernel', kernel); %, 'Support', 'positive');
    [~, yi_probabilityI] = max(yi_probability);
    yi_mode = yi_value(yi_probabilityI);
    
    % Store the trace baseline.
    baselines(keepI(i)) = yi_mode;
end

% Scale the traces.
offsets = nanmin(traces, [], 2);
traces = traces - offsets;
baselines = baselines - offsets;
scales = nanmax(traces, [], 2);
traces = traces ./ scales;
baselines = baselines ./ scales;

% Compute the compounded scales & offsets.
offsets = offsets + detrend_offsets;
end
